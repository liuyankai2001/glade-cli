import argparse
import importlib
import importlib.util
import unittest
import contextlib
import io
import socket
import tempfile
from pathlib import Path
from unittest.mock import patch


class WebCliTests(unittest.TestCase):
    def parser(self):
        spec = importlib.util.find_spec("src.cli.commands.web")
        self.assertIsNotNone(spec, "single-command web launcher must be registered")
        parser = argparse.ArgumentParser()
        importlib.import_module("src.cli.commands.web").register(
            parser.add_subparsers(dest="command", required=True)
        )
        return parser

    def test_single_filename_and_launch_defaults(self):
        args = self.parser().parse_args(["web", "-i", "项目.json"])
        self.assertEqual(args.input, "项目.json")
        self.assertEqual(args.port, 8000)
        self.assertFalse(args.no_browser)
        self.assertTrue(callable(args.func))

    def test_paths_and_invalid_ports_are_rejected_before_project_load(self):
        parser = self.parser()
        for filename in ["../other.json", "parts\\a.json", "C:\\other.json"]:
            with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                parser.parse_args(["web", "-i", filename])
        for port in ["0", "65536"]:
            with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                parser.parse_args(["web", "-i", "a.json", "--port", port])

    def test_ctrl_c_stops_launcher_without_a_traceback(self):
        from plasmid_fixtures import make_project

        with tempfile.TemporaryDirectory() as tmp:
            config = make_project(Path(tmp))
            with socket.socket() as listener:
                listener.bind(("127.0.0.1", 0))
                config.port = listener.getsockname()[1]
            config.no_browser = True
            self.parser()
            run_web = importlib.import_module("src.cli.commands.web").run_web
            escaped = False
            output = io.StringIO()
            with (
                patch("uvicorn.Server.run", side_effect=KeyboardInterrupt),
                contextlib.redirect_stdout(output),
            ):
                try:
                    run_web(config)
                except KeyboardInterrupt:
                    escaped = True
            self.assertFalse(escaped, "Ctrl+C must be a normal launcher exit")


if __name__ == "__main__":
    unittest.main()
