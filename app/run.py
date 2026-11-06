"""Launch the Terranova dashboard.

    python -m app.run              # http://127.0.0.1:8000
    python -m app.run --port 9000
    python -m app.run --reload     # development

Everything served is precomputed. If the pipeline has not been run, the app
still starts and the UI reports which artefacts are missing.
"""

from __future__ import annotations

import argparse
import webbrowser
from threading import Timer

import uvicorn


def main() -> None:
    parser = argparse.ArgumentParser(description="Terranova dashboard")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true",
                        help="auto-reload on code changes (development)")
    parser.add_argument("--no-browser", action="store_true",
                        help="do not open a browser window on start")
    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}"
    print(f"\n  Terranova dashboard -> {url}\n  (ctrl-c to stop)\n")

    if not args.no_browser and not args.reload:
        Timer(1.2, lambda: webbrowser.open(url)).start()

    uvicorn.run(
        "app.backend.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level="info",
    )


if __name__ == "__main__":
    main()
