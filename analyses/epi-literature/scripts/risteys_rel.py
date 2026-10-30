"""Fetch Risteys (FinnGen R13/FinRegistry) relationships table for an endpoint via Phoenix LiveView websocket.
Output: CSV with case overlap, FinRegistry survival HR (endpoint -> other? see Risteys docs), FinnGen rg, hits."""
import re, json, requests, websocket, os, sys, html, csv
ep = sys.argv[1]; outdir = sys.argv[2]; host = "risteys.finngen.fi"; url = f"https://{host}/endpoints/{ep}"
s = requests.Session(); page = s.get(url, timeout=120).text
open(f"{outdir}/{ep}.page.html", "w").write(page)
tok = re.search(r'csrf-token" content="([^"]+)"', page).group(1)
views = re.findall(r'data-phx-session="([^"]+)" data-phx-static="([^"]+)" id="([^"]+)"', page)
cookie = "; ".join(f"{k}={v}" for k, v in s.cookies.items())
port = int(os.environ["HTTPS_PROXY"].rsplit(":", 1)[1])
ws = websocket.create_connection(f"wss://{host}/live/websocket?_csrf_token={tok}&vsn=2.0.0", header=[f"Cookie: {cookie}"], origin=f"https://{host}",
    http_proxy_host="127.0.0.1", http_proxy_port=port, proxy_type="http", sslopt={"ca_certs": "/root/.ccr/ca-bundle.crt"}, timeout=60)
# the relationships view is the one whose join reply contains form-relationships-table
payload = None
for sess, stat, vid in views:
    ws.send(json.dumps(["1", "1", f"lv:{vid}", "phx_join", {"url": url, "params": {"_csrf_token": tok, "_mounts": 0}, "session": sess, "static": stat}]))
    m = json.loads(ws.recv())
    if "form-relationships-table" in json.dumps(m):
        while True:
            m2 = json.loads(ws.recv())
            if m2[3] == "diff" and "8" in m2[4]:
                payload = m2[4]; break
        break
ws.close()
json.dump(payload, open(f"{outdir}/{ep}.rel_raw.json", "w"))
tmpl = payload["8"].get("p", {})
def rend(x):
    if isinstance(x, str): return x
    st = tmpl[str(x["s"])] if isinstance(x.get("s"), int) else x["s"]
    out = st[0]
    for i in range(len(st) - 1): out += rend(x.get(str(i), "")) + st[i + 1]
    return out
rows = []
for r in payload["8"]["d"]:
    name = re.search(r'/endpoints/([^"]+)', r[0]).group(1)
    cells = [re.sub(r"<[^>]+>", " ", rend(c)) for c in r[3:]]
    cells = [html.unescape(re.sub(r"\s+", " ", c)).strip() for c in cells]
    rows.append([name, html.unescape(r[2])] + cells)
with open(f"{outdir}/{ep}.relationships.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["endpoint", "longname", "cases_FR", "cases_FG", "HR_FR", "rg_FG", "hits_FG", "coloc_FG"]); w.writerows(rows)
print(ep, len(rows)); print(json.dumps({k: v for k, v in tmpl.items()})[:1500])
