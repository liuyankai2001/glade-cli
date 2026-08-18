"""Query-scoped progress events for the command-line workflow.

Progress is deliberately separated from business output: console reporters
write to stderr, while the final text or JSON result remains on stdout.
"""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator, Iterator
from contextlib import asynccontextmanager, contextmanager, suppress
from contextvars import ContextVar
from dataclasses import dataclass
import sys
import time
from typing import Literal, Protocol, TextIO


ProgressStatus = Literal[
    "started",
    "completed",
    "info",
    "warning",
    "error",
    "heartbeat",
    "skipped",
]
ProgressMode = Literal["auto", "on", "off"]

_STAGE_LABELS = {
    "startup": "启动",
    "workflow": "工作流",
    "validation": "输入校验",
    "uniprot": "UniProt",
    "kegg": "KEGG",
    "reaction_match": "反应匹配",
    "auxiliary_requirement": "亚基分析",
    "routing": "路由",
    "research_init": "研究初始化",
    "database_retrieval": "数据库检索",
    "database_analysis": "数据库分析",
    "literature_retrieval": "文献检索",
    "literature_analysis": "文献分析",
    "web_retrieval": "网页兜底",
    "web_analysis": "网页分析",
    "host_retrieval": "宿主核验",
    "host_analysis": "宿主分析",
    "final_synthesis": "最终裁决",
    "research": "证据研究",
}
_STATUS_LABELS = {
    "warning": "警告",
    "error": "错误",
    "heartbeat": "等待",
    "skipped": "跳过",
}


@dataclass(frozen=True, slots=True)
class ProgressEvent:
    """One user-facing workflow update."""

    stage: str
    status: ProgressStatus
    message: str
    verbose_only: bool = False


class ProgressReporter(Protocol):
    """Minimal reporter interface injected through a context variable."""

    @property
    def enabled(self) -> bool: ...

    def emit(self, event: ProgressEvent) -> None: ...


class NullProgressReporter:
    """Reporter used for tests, libraries, redirects, and quiet mode."""

    @property
    def enabled(self) -> bool:
        return False

    def emit(self, event: ProgressEvent) -> None:
        _ = event


class MemoryProgressReporter:
    """Collect progress events for programmatic consumers and tests."""

    def __init__(self) -> None:
        self.events: list[ProgressEvent] = []

    @property
    def enabled(self) -> bool:
        return True

    def emit(self, event: ProgressEvent) -> None:
        self.events.append(event)


class ConsoleProgressReporter:
    """Render stable, line-oriented progress to a terminal stream."""

    def __init__(
        self,
        *,
        stream: TextIO | None = None,
        verbose: bool = False,
        clock=time.perf_counter,
    ) -> None:
        self._stream = stream or sys.stderr
        self._verbose = verbose
        self._clock = clock
        self._started = clock()

    @property
    def enabled(self) -> bool:
        return True

    def emit(self, event: ProgressEvent) -> None:
        if event.verbose_only and not self._verbose:
            return
        elapsed = max(0.0, self._clock() - self._started)
        stage = _STAGE_LABELS.get(event.stage, event.stage)
        status = _STATUS_LABELS.get(event.status)
        marker = f" [{status}]" if status else ""
        print(
            f"[{_format_elapsed(elapsed)}] [{stage}]{marker} {event.message}",
            file=self._stream,
            flush=True,
        )


_NULL_REPORTER = NullProgressReporter()
_CURRENT_REPORTER: ContextVar[ProgressReporter] = ContextVar(
    "protein_supply_progress_reporter",
    default=_NULL_REPORTER,
)


def build_console_progress_reporter(
    mode: ProgressMode,
    *,
    verbose: bool = False,
    quiet: bool = False,
    stream: TextIO | None = None,
) -> ProgressReporter:
    """Resolve CLI progress settings without changing final output."""

    effective_stream = stream or sys.stderr
    if quiet or mode == "off":
        return _NULL_REPORTER
    if mode == "auto" and not effective_stream.isatty():
        return _NULL_REPORTER
    return ConsoleProgressReporter(stream=effective_stream, verbose=verbose)


@contextmanager
def progress_reporter_context(
    reporter: ProgressReporter,
) -> Iterator[None]:
    """Install a reporter for the current synchronous/async query context."""

    token = _CURRENT_REPORTER.set(reporter)
    try:
        yield
    finally:
        _CURRENT_REPORTER.reset(token)


def progress_is_enabled() -> bool:
    return _CURRENT_REPORTER.get().enabled


def emit_progress(
    stage: str,
    status: ProgressStatus,
    message: str,
    *,
    verbose_only: bool = False,
) -> None:
    """Emit a best-effort event; progress failures never fail the workflow."""

    try:
        _CURRENT_REPORTER.get().emit(
            ProgressEvent(
                stage=stage,
                status=status,
                message=message,
                verbose_only=verbose_only,
            )
        )
    except Exception:
        return


@asynccontextmanager
async def progress_heartbeat(
    stage: str,
    message: str,
    *,
    interval_seconds: float = 10.0,
    timeout_seconds: float | None = None,
) -> AsyncIterator[None]:
    """Report elapsed waiting time for an operation with no true percentage."""

    if not progress_is_enabled():
        yield
        return

    stopped = asyncio.Event()
    started = time.perf_counter()

    async def heartbeat() -> None:
        while True:
            try:
                await asyncio.wait_for(
                    stopped.wait(),
                    timeout=interval_seconds,
                )
                return
            except TimeoutError:
                elapsed = int(time.perf_counter() - started)
                limit = (
                    f"，上限 {timeout_seconds:g} 秒"
                    if timeout_seconds is not None
                    else ""
                )
                emit_progress(
                    stage,
                    "heartbeat",
                    f"{message}，已等待 {elapsed} 秒{limit}",
                )

    task = asyncio.create_task(heartbeat())
    try:
        yield
    finally:
        stopped.set()
        task.cancel()
        with suppress(asyncio.CancelledError):
            await task


def _format_elapsed(seconds: float) -> str:
    total_seconds = int(seconds)
    minutes, remaining = divmod(total_seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours:02d}:{minutes:02d}:{remaining:02d}"
    return f"{minutes:02d}:{remaining:02d}"
