"""One serialized generation worker with browser-visible progress."""

from __future__ import annotations

import logging
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy

from src.plasmid_design.errors import DesignError


class JobManager:
    def __init__(self, service):
        self.service = service
        self._lock = threading.Lock()
        self._jobs = {}
        self._active_id = None
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="plasmid-design"
        )

    def submit(self, request: dict) -> dict:
        with self._lock:
            if self._active_id is not None:
                raise DesignError(
                    "已有生成任务正在运行，请等待完成。", code="project_busy"
                )
            job_id = uuid.uuid4().hex
            job = {
                "id": job_id,
                "status": "queued",
                "stage": "queued",
                "message": "任务已提交",
            }
            self._jobs[job_id] = job
            self._active_id = job_id
        try:
            # Validate IDs/revision before accepting; worker checks again before commit.
            preview = self.service.preview(request)
            if not preview["valid"]:
                raise DesignError(
                    "设计检查未通过，请按提示调整组件或限制酶。",
                    code="design_invalid",
                    issues=preview["issues"],
                )
            self._executor.submit(self._run, job_id, deepcopy(request))
        except Exception:
            with self._lock:
                self._active_id = None
                self._jobs.pop(job_id, None)
            raise
        return self.get(job_id)

    def _update(self, job_id, **changes):
        with self._lock:
            self._jobs[job_id].update(changes)

    def _run(self, job_id, request):
        self._update(job_id, status="running", stage="checking", message="正在检查设计")
        try:
            result = self.service.generate(
                request,
                progress=lambda stage, message: self._update(
                    job_id, stage=stage, message=message
                ),
            )
        except DesignError as exc:
            self._update(
                job_id,
                status="failed",
                stage="failed",
                message=str(exc),
                error={"code": exc.code, "message": str(exc)},
            )
        except Exception:
            logging.getLogger(__name__).exception("Plasmid generation failed")
            self._update(
                job_id,
                status="failed",
                stage="failed",
                message="生成失败，旧结果已保留。详细原因见命令行。",
                error={
                    "code": "generation_failed",
                    "message": "生成失败，旧结果已保留。请查看命令行后重试。",
                },
            )
        else:
            self._update(
                job_id,
                status="succeeded",
                stage="complete",
                message="质粒已生成并登记",
                result=result,
            )
        finally:
            with self._lock:
                self._active_id = None
                while len(self._jobs) > 32:
                    self._jobs.pop(next(iter(self._jobs)))

    def get(self, job_id: str) -> dict:
        with self._lock:
            if job_id not in self._jobs:
                raise DesignError(
                    "找不到该生成任务，请刷新页面查看已登记结果。", code="job_not_found"
                )
            return deepcopy(self._jobs[job_id])

    def active(self) -> dict | None:
        with self._lock:
            return (
                deepcopy(self._jobs[self._active_id])
                if self._active_id is not None
                else None
            )

    def close(self):
        self._executor.shutdown(wait=True, cancel_futures=True)
