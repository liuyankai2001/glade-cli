"""Research modes, budgets, and safe tool-call guards."""

from __future__ import annotations

import asyncio
from collections.abc import Awaitable, Callable, Mapping
from dataclasses import dataclass, field
import json
from typing import Any, Literal

from langchain.agents.middleware import (
    ModelCallLimitMiddleware,
    ToolCallLimitMiddleware,
)
from langchain.agents.middleware.types import AgentMiddleware, ToolCallRequest
from langchain_core.messages import ToolMessage
from langgraph.types import Command


ResearchMode = Literal["balanced", "deep"]
ResearcherRole = Literal[
    "supervisor",
    "bio_database",
    "literature",
    "web",
    "host_compatibility",
]


@dataclass(frozen=True, slots=True)
class ResearchBudget:
    """Limits applied to one agent invocation in balanced mode."""

    max_model_calls: int
    max_tool_calls: int
    default_tool_timeout_seconds: float
    per_tool_limits: Mapping[str, int] = field(default_factory=dict)
    per_tool_timeouts: Mapping[str, float] = field(default_factory=dict)
    grouped_tool_limits: tuple[tuple[frozenset[str], int], ...] = ()


@dataclass(frozen=True, slots=True)
class ResearchPolicy:
    """Complete runtime policy for one public research mode."""

    mode: ResearchMode
    workflow_timeout_seconds: float
    supervisor_max_tokens: int
    worker_max_tokens: int
    model_timeout_seconds: float
    web_search_mode: Literal["request", "auto", "playwright"]
    playwright_navigation_timeout_ms: int
    budgets: Mapping[ResearcherRole, ResearchBudget] | None

    def budget_for(self, role: ResearcherRole) -> ResearchBudget | None:
        """Return the role budget, or no budget in deep mode."""

        if self.budgets is None:
            return None
        return self.budgets[role]


BALANCED_RESEARCH_POLICY = ResearchPolicy(
    mode="balanced",
    workflow_timeout_seconds=300,
    supervisor_max_tokens=4096,
    worker_max_tokens=3200,
    model_timeout_seconds=60,
    web_search_mode="request",
    playwright_navigation_timeout_ms=8000,
    budgets={
        "supervisor": ResearchBudget(
            max_model_calls=8,
            max_tool_calls=5,
            default_tool_timeout_seconds=100,
            per_tool_limits={"task": 5},
        ),
        "bio_database": ResearchBudget(
            max_model_calls=8,
            max_tool_calls=8,
            default_tool_timeout_seconds=15,
            per_tool_limits={
                "intact_search_interactions": 1,
                "intact_get_interactions": 1,
            },
            per_tool_timeouts={
                "intact_search_interactions": 12,
                "intact_get_interactions": 12,
            },
            grouped_tool_limits=(
                (
                    frozenset(
                        {
                            "intact_search_interactions",
                            "intact_get_interactions",
                        }
                    ),
                    1,
                ),
            ),
        ),
        "literature": ResearchBudget(
            max_model_calls=8,
            max_tool_calls=8,
            default_tool_timeout_seconds=20,
            per_tool_limits={
                "SemanticScholar_search_papers": 1,
                "SemanticScholar_get_pdf_snippets": 1,
            },
            per_tool_timeouts={
                "SemanticScholar_search_papers": 15,
                "SemanticScholar_get_pdf_snippets": 15,
            },
            grouped_tool_limits=(
                (
                    frozenset(
                        {
                            "SemanticScholar_search_papers",
                            "SemanticScholar_get_pdf_snippets",
                        }
                    ),
                    1,
                ),
            ),
        ),
        "web": ResearchBudget(
            max_model_calls=6,
            max_tool_calls=4,
            default_tool_timeout_seconds=15,
            per_tool_limits={"search": 2, "fetchWebContent": 2},
        ),
        "host_compatibility": ResearchBudget(
            max_model_calls=8,
            max_tool_calls=8,
            default_tool_timeout_seconds=30,
        ),
    },
)


DEEP_RESEARCH_POLICY = ResearchPolicy(
    mode="deep",
    workflow_timeout_seconds=900,
    supervisor_max_tokens=4096,
    worker_max_tokens=4096,
    model_timeout_seconds=90,
    # Deep mode broadens sources and query depth, but ordinary primary-source
    # pages do not require a browser engine. Request mode avoids paying a
    # Playwright startup cost on every query-scoped research session.
    web_search_mode="request",
    playwright_navigation_timeout_ms=8000,
    budgets=None,
)


def get_research_policy(mode: ResearchMode) -> ResearchPolicy:
    """Resolve one validated public research mode."""

    if mode == "balanced":
        return BALANCED_RESEARCH_POLICY
    if mode == "deep":
        return DEEP_RESEARCH_POLICY
    raise ValueError(f"unsupported research mode: {mode}")


class ToolExecutionGuardMiddleware(AgentMiddleware):
    """Deduplicate tool calls and convert timeouts into auditable failures."""

    def __init__(
        self,
        *,
        default_timeout_seconds: float,
        per_tool_timeouts: Mapping[str, float] | None = None,
        grouped_tool_limits: tuple[tuple[frozenset[str], int], ...] = (),
    ) -> None:
        super().__init__()
        self.default_timeout_seconds = default_timeout_seconds
        self.per_tool_timeouts = dict(per_tool_timeouts or {})
        self.grouped_tool_limits = grouped_tool_limits
        self._group_counts = [0 for _ in grouped_tool_limits]
        self._completed: dict[str, ToolMessage] = {}
        self._inflight: dict[str, asyncio.Future[ToolMessage]] = {}
        self._lock = asyncio.Lock()

    async def awrap_tool_call(
        self,
        request: ToolCallRequest,
        handler: Callable[
            [ToolCallRequest],
            Awaitable[ToolMessage | Command[Any]],
        ],
    ) -> ToolMessage | Command[Any]:
        tool_name = request.tool_call.get("name") or "unknown_tool"
        cache_key = _tool_cache_key(tool_name, request.tool_call.get("args"))

        async with self._lock:
            cached = self._completed.get(cache_key)
            if cached is not None:
                return _clone_tool_message(cached, request)

            inflight = self._inflight.get(cache_key)
            if inflight is None:
                exceeded_group = self._exceeded_group(tool_name)
                if exceeded_group is not None:
                    return ToolMessage(
                        content=(
                            f"Tool group budget was exhausted before {tool_name}. "
                            "Use the evidence already collected, record the source "
                            "as unavailable, and do not treat this as negative "
                            "biological evidence."
                        ),
                        tool_call_id=request.tool_call["id"],
                        name=tool_name,
                        status="error",
                    )
                inflight = asyncio.get_running_loop().create_future()
                self._inflight[cache_key] = inflight
                self._increment_groups(tool_name)
                owns_call = True
            else:
                owns_call = False

        if not owns_call:
            result = await inflight
            return _clone_tool_message(result, request)

        timeout_seconds = self.per_tool_timeouts.get(
            tool_name,
            self.default_timeout_seconds,
        )
        try:
            async with asyncio.timeout(timeout_seconds):
                result = await handler(request)
            if not isinstance(result, ToolMessage):
                return result
        except TimeoutError:
            result = ToolMessage(
                content=(
                    f"Tool {tool_name} timed out after {timeout_seconds:g} seconds. "
                    "Record this as an unavailable source; it is not negative "
                    "biological evidence and must not be retried in this run."
                ),
                tool_call_id=request.tool_call["id"],
                name=tool_name,
                status="error",
            )
        except Exception as exc:
            async with self._lock:
                waiter = self._inflight.pop(cache_key)
                if not waiter.done():
                    waiter.set_exception(exc)
                    waiter.exception()
            raise

        async with self._lock:
            self._completed[cache_key] = result
            waiter = self._inflight.pop(cache_key)
            if not waiter.done():
                waiter.set_result(result)
        return result

    def _exceeded_group(self, tool_name: str) -> int | None:
        for index, (tool_names, limit) in enumerate(self.grouped_tool_limits):
            if tool_name in tool_names and self._group_counts[index] >= limit:
                return index
        return None

    def _increment_groups(self, tool_name: str) -> None:
        for index, (tool_names, _) in enumerate(self.grouped_tool_limits):
            if tool_name in tool_names:
                self._group_counts[index] += 1


def build_budget_middleware(
    budget: ResearchBudget | None,
) -> list[AgentMiddleware]:
    """Build the middleware stack for one balanced-mode researcher."""

    if budget is None:
        return []

    middleware: list[AgentMiddleware] = [
        ToolExecutionGuardMiddleware(
            default_timeout_seconds=budget.default_tool_timeout_seconds,
            per_tool_timeouts=budget.per_tool_timeouts,
            grouped_tool_limits=budget.grouped_tool_limits,
        ),
        ToolCallLimitMiddleware(
            run_limit=budget.max_tool_calls,
            exit_behavior="continue",
        ),
    ]
    middleware.extend(
        ToolCallLimitMiddleware(
            tool_name=tool_name,
            run_limit=limit,
            exit_behavior="continue",
        )
        for tool_name, limit in budget.per_tool_limits.items()
    )
    middleware.append(
        ModelCallLimitMiddleware(
            run_limit=budget.max_model_calls,
            exit_behavior="end",
        )
    )
    return middleware


def _tool_cache_key(tool_name: str, arguments: Any) -> str:
    try:
        serialized = json.dumps(
            arguments,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        )
    except (TypeError, ValueError):
        serialized = repr(arguments)
    return f"{tool_name}:{serialized}"


def _clone_tool_message(
    message: ToolMessage,
    request: ToolCallRequest,
) -> ToolMessage:
    return ToolMessage(
        content=message.content,
        additional_kwargs=message.additional_kwargs,
        response_metadata=message.response_metadata,
        name=request.tool_call.get("name") or message.name,
        tool_call_id=request.tool_call["id"],
        artifact=message.artifact,
        status=message.status,
    )
