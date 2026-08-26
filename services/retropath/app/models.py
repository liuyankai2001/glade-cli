from __future__ import annotations

from pydantic import BaseModel, Field, model_validator


class JobParameters(BaseModel):
    max_steps: int = Field(default=3, ge=1, le=10)
    topx: int = Field(default=100, ge=1, le=1000)
    dmin: int = Field(default=2, ge=0, le=16)
    dmax: int = Field(default=16, ge=2, le=16)
    mwmax_source: int = Field(default=1000, ge=1, le=5000)
    msc_timeout: int = Field(default=10, ge=1, le=60)

    @model_validator(mode="after")
    def validate_diameters(self) -> "JobParameters":
        if self.dmin > self.dmax:
            raise ValueError("dmin must be less than or equal to dmax")
        return self


TERMINAL_STATUSES = {
    "succeeded",
    "no_solution",
    "source_in_sink",
    "failed",
    "timed_out",
}


def public_job(job: dict[str, object]) -> dict[str, object]:
    return {
        "job_id": job["job_id"],
        "status": job["status"],
        "created_at": job["created_at"],
        "started_at": job.get("started_at"),
        "finished_at": job.get("finished_at"),
        "return_code": job.get("return_code"),
        "failure_code": job.get("failure_code"),
        "error": job.get("error"),
        "parameters": job["parameters"],
    }
