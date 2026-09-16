"""Start the bundled local plasmid design workbench with one command."""

from __future__ import annotations

from argparse import ArgumentTypeError


def _filename(value: str) -> str:
    if (
        not value
        or value in (".", "..")
        or any(c in value for c in ("/", "\\", ":"))
        or not value.lower().endswith(".json")
    ):
        raise ArgumentTypeError("--input 只需填写 inputs 目录下的 JSON 文件名。")
    return value


def _port(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise ArgumentTypeError("端口必须为 1–65535 的整数。") from exc
    if not 1 <= number <= 65535:
        raise ArgumentTypeError("端口必须为 1–65535 的整数。")
    return number


def register(subparsers):
    parser = subparsers.add_parser("web", help="启动本地质粒设计网页")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=_filename,
        help="inputs 目录下的项目 JSON 文件名",
    )
    parser.add_argument("--port", type=_port, default=8000, help="本地端口，默认 8000")
    parser.add_argument(
        "--no-browser", action="store_true", help="启动服务但不自动打开浏览器"
    )
    parser.set_defaults(func=run_web)


def run_web(config):
    import socket
    import webbrowser

    import uvicorn

    from src.pathway_analyze.target_id import validate_target_compound_id
    from src.web_api.app import create_app

    validate_target_compound_id(config.target_name)
    port = int(getattr(config, "port", 8000))
    url = f"http://127.0.0.1:{port}"
    app = create_app(config)

    class LocalServer(uvicorn.Server):
        async def startup(self, sockets=None):
            await super().startup(sockets=sockets)
            if self.started:
                print(f"GLADE 质粒设计：{url}\n按 Ctrl+C 停止服务。", flush=True)
                if not getattr(config, "no_browser", False):
                    try:
                        webbrowser.open(url)
                    except (webbrowser.Error, OSError):
                        print(f"浏览器未自动打开，请访问 {url}", flush=True)

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        try:
            listener.bind(("127.0.0.1", port))
        except OSError as exc:
            raise SystemExit(f"无法使用端口 {port}，请换一个 --port 后重试。") from exc
        server = LocalServer(
            uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning")
        )
        try:
            server.run(sockets=[listener])
        except KeyboardInterrupt:
            print("GLADE 网页服务已停止。", flush=True)
