"""Run the Predictive Surf Atlas locally: `python -m surfmap gui`.

Serves web/ (the same page that is published as a shareable link) and adds three actions
the page shows only here: download new data, rebuild the PDF, and refresh the app's own data
at the local resolution. Each runs as a background `python -m surfmap ...` job whose log the
page polls.
"""

import http.server
import json
import subprocess
import sys
import threading
import webbrowser

from . import config

WEB = config.ROOT / "web"
HEAD = (b'<!doctype html><html lang="en"><head><meta charset="utf-8">'
        b'<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"></head><body>')


class _Jobs:
    def __init__(self, res):
        self.res, self.lock, self.lines, self.busy, self.reload = res, threading.Lock(), [], False, False

    def steps(self, name):
        py = [sys.executable, "-m", "surfmap"]
        return {
            "fetch": [py + ["fetch", "--source", s] for s in ("ww3", "points", "copernicus")],
            "build": [py + ["build"]],
            "export": [py + ["export", "--res", str(self.res)]],
        }[name]

    def start(self, name):
        with self.lock:
            if self.busy:
                return False
            self.busy, self.reload, self.lines = True, False, [f"$ {name}"]
        threading.Thread(target=self._run, args=(name,), daemon=True).start()
        return True

    def _run(self, name):
        ok = True
        for cmd in self.steps(name):
            self._say("$ " + " ".join(cmd[2:]))
            proc = subprocess.Popen(cmd, cwd=config.ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in proc.stdout:
                if "warn" not in line.lower():
                    self._say(line.rstrip())
            ok = proc.wait() == 0 and ok
        with self.lock:
            self.busy = False
            self.reload = ok and name == "export"
            self.lines.append("Done." if ok else "Finished with errors; see above.")

    def _say(self, line):
        with self.lock:
            self.lines = (self.lines + [line])[-200:]

    def state(self):
        with self.lock:
            return dict(busy=self.busy, reload=self.reload, log="\n".join(self.lines[-60:]))


def serve(port=8765, res=1, open_browser=True):
    if not (WEB / "data" / "meta.json").exists():
        print(f"No app data yet; exporting at {res} degrees first.")
        from . import export
        export.export(res)
    jobs = _Jobs(res)

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k):
            super().__init__(*a, directory=str(WEB), **k)

        def _json(self, obj):
            body = json.dumps(obj).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path in ("/", "/index.html"):  # published pages get this skeleton from the host
                body = HEAD + (WEB / "index.html").read_bytes()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            elif self.path == "/api/ping":
                self._json(dict(local=True))
            elif self.path == "/api/log":
                self._json(jobs.state())
            else:
                super().do_GET()

        def do_POST(self):
            name = self.path.rsplit("/", 1)[-1]
            if self.path.startswith("/api/run/") and name in ("fetch", "build", "export"):
                self._json(dict(started=jobs.start(name)))
            else:
                self.send_error(404)

        def log_message(self, *args):
            pass

    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}/"
    print(f"Predictive Surf Atlas at {url}  (Ctrl-C to stop)")
    if open_browser:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
