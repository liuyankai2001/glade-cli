"""Same-origin API and bundled frontend for a project selected by the CLI."""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from filelock import Timeout
from pydantic import BaseModel, ConfigDict, Field, StrictInt, field_validator
from starlette.exceptions import HTTPException
from starlette.middleware.trustedhost import TrustedHostMiddleware

from src.plasmid_design.errors import DesignError
from src.plasmid_design.sequence import (
    DEFAULT_COMPONENT_ORDER,
    normalize_component_order,
)
from src.plasmid_design.service import DesignService
from src.web_api.jobs import JobManager


class DesignRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)
    resistance_id: str = Field(min_length=1, max_length=120)
    replication_id: str = Field(min_length=1, max_length=120)
    t0_id: str | None = Field(default=None, min_length=1, max_length=120)
    t1_id: str | None = Field(default=None, min_length=1, max_length=120)
    expected_revision: StrictInt = Field(ge=0)
    source_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    component_order: list[
        Literal["resistance", "replication", "t0", "t1", "expression"]
    ] = Field(
        default_factory=lambda: list(DEFAULT_COMPONENT_ORDER),
        min_length=5,
        max_length=5,
    )

    @field_validator("component_order", mode="before")
    @classmethod
    def expand_legacy_order(cls, value):
        return normalize_component_order(value)


def error_response(code, message, status, issues=None):
    error = {"code": code, "message": message}
    if issues:
        error["issues"] = issues
    return JSONResponse({"error": error}, status_code=status)


def create_app(config) -> FastAPI:
    service = DesignService(config)
    jobs = JobManager(service)

    @asynccontextmanager
    async def lifespan(app):
        yield
        jobs.close()

    app = FastAPI(
        title="GLADE 质粒设计",
        lifespan=lifespan,
        docs_url=None,
        redoc_url=None,
        openapi_url=None,
    )
    app.state.design_service = service
    app.state.jobs = jobs
    app.add_middleware(TrustedHostMiddleware, allowed_hosts=["127.0.0.1", "localhost"])

    @app.middleware("http")
    async def local_origin(request: Request, call_next):
        if request.method == "POST" and request.headers.get("origin") not in (
            None,
            str(request.base_url).rstrip("/"),
        ):
            return error_response(
                "invalid_origin", "该本地操作需从 GLADE 页面发起。", 403
            )
        response = await call_next(request)
        if request.url.path == "/" or request.url.path.startswith("/api/"):
            response.headers["Cache-Control"] = "no-store"
        return response

    @app.exception_handler(DesignError)
    async def design_error(request, exc):
        if exc.code in {
            "stale_source",
            "stale_result",
            "artifact_invalid",
            "project_busy",
        }:
            status = 409
        elif exc.code in {"artifact_not_found", "job_not_found"}:
            status = 404
        else:
            status = 422
        return error_response(exc.code, str(exc), status, exc.issues)

    @app.exception_handler(RequestValidationError)
    async def request_error(request, exc):
        return error_response(
            "invalid_request", "请求参数不完整或格式有误，请刷新页面后重试。", 422
        )

    @app.exception_handler(Timeout)
    async def lock_error(request, exc):
        return error_response(
            "project_busy", "项目正在被其他命令更新，请稍后重试。", 409
        )

    @app.exception_handler(ValueError)
    async def data_error(request, exc):
        return error_response(
            "invalid_project_data", f"项目或组件库数据无效：{exc}", 422
        )

    @app.exception_handler(HTTPException)
    async def missing_route(request, exc):
        return error_response("not_found", "找不到该页面或文件。", exc.status_code)

    @app.exception_handler(Exception)
    async def server_error(request, exc):
        logging.getLogger(__name__).error("Local web request failed", exc_info=exc)
        return error_response(
            "server_error", "本地服务操作失败，详细原因见命令行。", 500
        )

    @app.get("/api/context")
    def context():
        result = service.context()
        result["active_job"] = jobs.active()
        return result

    @app.get("/api/modules")
    def modules():
        return service.modules()

    @app.post("/api/preview")
    def preview(request: DesignRequest):
        return service.preview(request.model_dump())

    @app.post("/api/generate", status_code=202)
    def generate(request: DesignRequest):
        return jobs.submit(request.model_dump())

    @app.get("/api/jobs/{job_id}")
    def job(job_id: str):
        return jobs.get(job_id)

    @app.get("/api/files/{generation_id}/{file_id}")
    def download(generation_id: str, file_id: str):
        path = service.download(generation_id, file_id)
        return FileResponse(
            path, filename=path.name, media_type="application/octet-stream"
        )

    @app.get("/api/{unknown:path}")
    def unknown_api(unknown: str):
        raise DesignError("找不到该接口。", code="artifact_not_found")

    static = Path(config.root_dir) / "front" / "static"
    if not (static / "index.html").is_file():
        raise DesignError(
            "前端静态文件缺失，请在 front 目录运行 npm ci 和 npm run build。",
            code="frontend_missing",
        )
    app.mount("/", StaticFiles(directory=static, html=True), name="frontend")
    return app
