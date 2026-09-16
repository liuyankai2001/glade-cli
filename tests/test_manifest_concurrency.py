import json
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from src.write_manifest.store import (
    manifest_update_lock,
    read_design_manifest,
    update_design_manifest,
)


class ManifestConcurrencyTests(unittest.TestCase):
    def test_reader_waits_for_writer_commit_and_reads_new_revision(self):
        # Windows cannot replace a JSON file while another thread has it open.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            update_design_manifest(
                path, target_compound_id="C00031", sections={"value": 1}
            )
            started, finished = threading.Event(), threading.Event()
            received = []

            def read():
                started.set()
                received.append(read_design_manifest(path))
                finished.set()

            with manifest_update_lock(path):
                worker = threading.Thread(target=read)
                worker.start()
                self.assertTrue(started.wait(1))
                completed_before_commit = finished.wait(0.05)
                update_design_manifest(
                    path, target_compound_id="C00031", sections={"value": 2}
                )
            worker.join(timeout=2)
            self.assertFalse(
                completed_before_commit,
                "reader must not keep the file open across a commit",
            )
            self.assertTrue(finished.is_set())
            self.assertEqual(received[0]["value"], 2)

    def test_only_one_process_revision_writer_can_commit(self):
        # Simultaneous writes with revision=0 must not both report success.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            barrier = threading.Barrier(2)

            def write(name):
                barrier.wait()
                try:
                    return update_design_manifest(
                        path,
                        target_compound_id="C00031",
                        sections={name: {}},
                        expected_revision=0,
                    )
                except Exception as exc:
                    return exc

            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(write, ["first", "second"]))
            self.assertEqual(sum(isinstance(r, dict) for r in results), 1)
            self.assertTrue(all(isinstance(r, (dict, ValueError)) for r in results))
            data = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(data["revision"], 1)
            self.assertEqual(sum(k in data for k in ["first", "second"]), 1)


if __name__ == "__main__":
    unittest.main()
