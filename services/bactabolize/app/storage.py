import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


def utcnow():
    return datetime.now(timezone.utc).isoformat()


class JobStorage:
    def __init__(self, jobs_dir):
        self.jobs_dir = Path(jobs_dir).resolve()
        self.database_path = self.jobs_dir / 'jobs.sqlite3'

    def connect(self):
        connection = sqlite3.connect(str(self.database_path), timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    def initialize(self):
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.execute('PRAGMA journal_mode=WAL')
            db.execute('CREATE TABLE IF NOT EXISTS jobs (job_id TEXT PRIMARY KEY, status TEXT, parameters TEXT, created_at TEXT, started_at TEXT, finished_at TEXT, return_code INTEGER, failure_code TEXT, error TEXT, job_dir TEXT)')
            db.execute("UPDATE jobs SET status='failed', finished_at=?, failure_code='service_restarted', error='service restarted during execution' WHERE status='running'", (utcnow(),))

    def get_job(self, job_id):
        with self.connect() as db:
            row = db.execute('SELECT * FROM jobs WHERE job_id=?',(job_id,)).fetchone()
        if row is None:
            return None
        result = dict(row)
        result['parameters'] = json.loads(result['parameters'])
        return result

    def active_count(self):
        with self.connect() as db:
            return db.execute("SELECT COUNT(*) FROM jobs WHERE status IN ('queued','running')").fetchone()[0]

    def queued_job_ids(self):
        with self.connect() as db:
            return [r[0] for r in db.execute("SELECT job_id FROM jobs WHERE status='queued' ORDER BY created_at")]

    def create_job(self, job_id, parameters, job_dir):
        with self.connect() as db:
            db.execute("INSERT INTO jobs (job_id,status,parameters,created_at,job_dir) VALUES (?,'queued',?,?,?)", (job_id,json.dumps(parameters),utcnow(),str(job_dir)))
        return self.get_job(job_id)

    def claim_job(self, job_id):
        with self.connect() as db:
            changed = db.execute("UPDATE jobs SET status='running',started_at=? WHERE job_id=? AND status='queued'",(utcnow(),job_id)).rowcount
        return self.get_job(job_id) if changed else None

    def finish_job(self, job_id, status, return_code=None, failure_code=None, error=None):
        with self.connect() as db:
            db.execute('UPDATE jobs SET status=?, finished_at=?,return_code=?,failure_code=?,error=? WHERE job_id=?',(status,utcnow(),return_code,failure_code,error,job_id))
