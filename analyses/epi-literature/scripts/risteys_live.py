"""Fetch Risteys LiveView-rendered sections (relationships table etc.) for an endpoint.
Usage: python3 risteys_live.py ENDPOINT outdir"""
import sys, re, json, requests, websocket, time, os, html
ep = sys.argv[1]; out = sys.argv[2]
host = sys.argv[3] if len(sys.argv) > 3 else "risteys.finngen.fi"
url = f"https://{host}/endpoints/{ep}"
s = requests.Session()
r = s.get(url, timeout=120)
page = r.text
open(f"{out}/{ep}.page.html", "w").write(page)
tok = re.search(r'csrf-token" content="([^"]+)"', page).group(1)
views = re.findall(r'data-phx-session="([^"]+)" data-phx-static="([^"]+)" id="([^"]+)"', page)
print("views", len(views), list(s.cookies.keys()))
cookie = "; ".join(f"{k}={v}" for k, v in s.cookies.items())
prox = os.environ.get("HTTPS_PROXY", "").replace("http://", "")
ph, pp = prox.split(":")
results = {}
def render(node, comps):
    # render phoenix liveview rendered tree to string
    if isinstance(node, (str, int, float)):
        return str(node)
    if isinstance(node, int):
        return render(comps.get(str(node), {}), comps)
    if node is None: return ""
    if "s" in node and isinstance(node["s"], list):
        st = node["s"]
        if "d" in node:  # comprehension
            outp = []
            for row in node["d"]:
                parts = [st[0]]
                for i, v in enumerate(row):
                    parts.append(render(v, comps)); parts.append(st[i + 1])
                outp.append("".join(parts))
            return "".join(outp)
        parts = [st[0]]
        for i in range(len(st) - 1):
            v = node.get(str(i), "")
            if isinstance(v, int) and not isinstance(v, bool):
                v = comps.get(str(v), {})
            parts.append(render(v, comps)); parts.append(st[i + 1])
        return "".join(parts)
    if "d" in node:
        return "".join(render(x, comps) for row in node["d"] for x in row)
    return ""
def merge(a, b):
    if not isinstance(a, dict) or not isinstance(b, dict): return b
    for k, v in b.items():
        if k in a and isinstance(a[k], dict) and isinstance(v, dict) and "s" not in v:
            merge(a[k], v)
        else:
            a[k] = v
    return a
for sess, stat, vid in views:
    ws = websocket.create_connection(
        f"wss://{host}/live/websocket?_csrf_token={tok}&vsn=2.0.0",
        header=[f"Cookie: {cookie}"], origin=f"https://{host}",
        http_proxy_host=ph, http_proxy_port=int(pp), proxy_type="http",
        sslopt={"ca_certs": "/root/.ccr/ca-bundle.crt"}, timeout=90)
    ws.send(json.dumps(["1", "1", f"lv:{vid}", "phx_join", {"url": url, "params": {"_csrf_token": tok, "_mounts": 0}, "session": sess, "static": stat}]))
    tree = None; comps = {}
    t0 = time.time()
    while time.time() - t0 < 60:
        try:
            msg = json.loads(ws.recv())
        except Exception as e:
            break
        jr, ref, topic, event, payload = msg
        if event == "phx_reply":
            resp = payload.get("response", {})
            if "rendered" in resp:
                tree = resp["rendered"]; comps = tree.pop("c", {}) if isinstance(tree, dict) else {}
            elif payload.get("status") != "ok":
                print("reply", payload); break
        elif event == "diff":
            c = payload.pop("c", {})
            merge(comps, c)
            if tree is not None: merge(tree, payload)
        # wait a few seconds for async loads
        if tree is not None and time.time() - t0 > 15: break
    ws.close()
    if tree is not None:
        results[vid] = render(tree, comps)
        print(vid, len(results[vid]))
with open(f"{out}/{ep}.live.html", "w") as f:
    for k, v in results.items():
        f.write(f"<!-- VIEW {k} -->\n{v}\n")
