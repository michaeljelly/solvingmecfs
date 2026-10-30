#!/usr/bin/env python3
"""Generate the site's figures as inline SVG fragments (colours from CSS variables). usage: python3 -I gen_figs.py <svg_out_dir>"""
import math, os, sys
OUT = sys.argv[1]
def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

STYLE_COMMON = """
.{c} text{{font-family:var(--f-sans);fill:var(--ink)}}
.{c} .nm{{font-size:14px;font-weight:600}}
.{c} .gl{{font-size:12.5px;fill:var(--muted)}}
.{c} .an{{font-size:12.5px;font-style:italic;fill:var(--muted)}}
.{c} .anw{{font-size:12.5px;font-style:italic;fill:var(--warm)}}
.{c} .hw{{font-size:13px;font-weight:600;fill:var(--warm)}}
.{c} .ln{{fill:none;stroke:var(--faint);stroke-width:1.2}}
.{c} .lw{{fill:none;stroke:var(--warm);stroke-width:2}}
.{c} .ah{{fill:var(--faint)}} .{c} .aw{{fill:var(--warm)}}
.{c} .dot{{fill:var(--bg);stroke:var(--ink);stroke-width:1.2}}
.{c} .dotw{{fill:var(--bg);stroke:var(--warm);stroke-width:1.6}}
.{c} .dn{{font-size:11px;font-weight:600;text-anchor:middle;dominant-baseline:central}}
.{c} .dnw{{fill:var(--warm)}}
.{c} a{{cursor:pointer;text-decoration:none}} .{c} a text{{text-decoration:none}} .{c} a:hover .nm{{fill:var(--link)}} .{c} a:focus-visible .dot{{stroke:var(--link);stroke-width:2.5}}
"""

def arrowhead(x, y, ang, cls="ah", s=6):
    a = math.radians(ang)
    p = [(x, y), (x - s * math.cos(a) + s * .55 * math.sin(a), y - s * math.sin(a) - s * .55 * math.cos(a)),
         (x - s * math.cos(a) - s * .55 * math.sin(a), y - s * math.sin(a) + s * .55 * math.cos(a))]
    return f'<path class="{cls}" d="M{p[0][0]:.1f} {p[0][1]:.1f}L{p[1][0]:.1f} {p[1][1]:.1f}L{p[2][0]:.1f} {p[2][1]:.1f}Z"/>'

import json
GRAPH = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "graph.json"), encoding="utf8"))
NODES = {s["n"]: (s["name"], s["gloss"]) for s in GRAPH["steps"]}
WEAK = {s["n"] for s in GRAPH["steps"] if s["found"] in ("early", "lc-only", "inferred", "contested")}
UNTESTED = {(e["from"], e["to"]) for e in GRAPH["edges"] if e["status"] == "untested"}
def dash(a, b): return ' style="stroke-dasharray:2 4"' if (a, b) in UNTESTED else ""
CRASH = {5, 6, 7}

def node(n, x, y, anchor, dy_name, dy_gloss, tx):
    nm, gl = NODES[n]
    w = n in CRASH and n != 5
    dc, nc = ("dotw", "dn dnw") if n in CRASH else ("dot", "dn")
    ds = ' style="stroke-dasharray:3 2.2"' if n in WEAK else ""
    return (f'<a href="/model#step-{n}"><circle class="{dc}"{ds} cx="{x:.1f}" cy="{y:.1f}" r="10"/><text class="{nc}" x="{x:.1f}" y="{y:.1f}">{n}</text>'
            f'<text class="nm" x="{tx:.1f}" y="{y+dy_name:.1f}" text-anchor="{anchor}">{esc(nm)}</text>'
            f'<text class="gl" x="{tx:.1f}" y="{y+dy_gloss:.1f}" text-anchor="{anchor}">{esc(gl)}</text></a>')

def model_wide():
    c = "fm-w"; W, H = 960, 500; cx, cy, rx, ry = 520, 255, 190, 185
    angs = {2: 180, 3: 128.6, 4: 77.1, 5: 25.7, 6: -25.7, 7: -77.1, 8: -128.6}
    P = {n: (cx + rx * math.cos(math.radians(a)), cy - ry * math.sin(math.radians(a))) for n, a in angs.items()}
    out = []
    order = [2, 3, 4, 5, 6, 7, 8, 2]
    for i in range(7):
        a0, a1 = angs[order[i]], angs[order[i + 1]] if order[i + 1] != 2 else -180
        warm = order[i] in (5, 6)
        # arc from a0 to a1 (decreasing angle = clockwise on screen), trimmed by 14 deg
        s, e = a0 - 9, a1 + 9
        pts = []
        for k in range(31):
            t = s + (e - s) * k / 30
            pts.append((cx + rx * math.cos(math.radians(t)), cy - ry * math.sin(math.radians(t))))
        d = "M" + " L".join(f"{x:.1f} {y:.1f}" for x, y in pts)
        out.append(f'<path class="{"lw" if warm else "ln"}"{dash(order[i], order[i + 1])} d="{d}"/>')
        (x1, y1), (x2, y2) = pts[-2], pts[-1]
        out.append(arrowhead(x2, y2, math.degrees(math.atan2(y2 - y1, x2 - x1)), "aw" if warm else "ah", 7 if warm else 6))
    # trigger and node 1
    out.append('<text class="an" x="16" y="168">viral infection (EBV, SARS-CoV-2, others)</text>')
    out.append('<path class="ln" style="stroke-dasharray:3 3" d="M40 178 V214"/>' + arrowhead(40, 222, 90))
    out.append(node(1, 40, 240, "start", 0, 17, 58))
    out.append(f'<path class="ln"{dash(1, 2)} d="M52 266 Q 200 300 312 262"/>' + arrowhead(316, 260, -20))
    # nodes 2..8
    x, y = P[2]; out.append(node(2, x, y, "start", 5, 22, x + 18))
    x, y = P[3]; out.append(node(3, x, y, "end", 0, 17, x - 18))
    x, y = P[4]; out.append(node(4, x, y, "middle", -36, -20, x))
    x, y = P[5]; out.append(node(5, x, y, "start", 0, 17, x + 18))
    x, y = P[6]; out.append(node(6, x, y, "start", 0, 17, x + 18))
    x, y = P[7]; out.append(node(7, x, y, "middle", 32, 48, x))
    x, y = P[8]; out.append(node(8, x, y, "end", 0, 17, x - 18))
    # exertion input to 5
    x5, y5 = P[5]
    out.append(f'<text class="anw" x="{x5+18:.1f}" y="{y5-58:.1f}">exertion, standing or sensory load</text>')
    out.append(f'<path class="lw" style="stroke-width:1.2" d="M{x5+8:.1f} {y5-52:.1f} L{x5+2:.1f} {y5-16:.1f}"/>' + arrowhead(x5 + 1.4, y5 - 13, 99, "aw"))
    # crash annotation inside the ring
    out.append(f'<text class="hw" x="{cx+150:.1f}" y="{cy+72:.1f}" text-anchor="end">post-exertional malaise</text>')
    out.append(f'<text class="anw" x="{cx+150:.1f}" y="{cy+89:.1f}" text-anchor="end">5 → 6 → 7, peaks 24–72 h</text>')
    body = "".join(out)
    return (f'<svg class="{c}" viewBox="0 0 {W} {H}" role="group" aria-labelledby="{c}-t" xmlns="http://www.w3.org/2000/svg"><title id="{c}-t">A working model of ME/CFS: eight abnormalities arranged as a loop</title>'
            f'<style>{STYLE_COMMON.format(c=c)}</style>{body}</svg>')

def model_narrow():
    c = "fm-n"; W = 340; y0 = 70; step = 58; H = y0 + step * 7 + 50
    out = ['<text class="an" x="10" y="20">viral infection (EBV, SARS-CoV-2, others)</text>',
           '<path class="ln" style="stroke-dasharray:3 3" d="M30 28 V50"/>' + arrowhead(30, 57, 90)]
    for i, n in enumerate(range(1, 9)):
        y = y0 + i * step
        out.append(node(n, 30, y, "start", 0, 17, 50))
        if n < 8:
            warm = n in (5, 6)
            out.append(f'<path class="{"lw" if warm else "ln"}"{dash(n, n + 1)} d="M30 {y+28} V{y+step-15}"/>' + arrowhead(30, y + step - 12, 90, "aw" if warm else "ah"))
    y2, y8 = y0 + step, y0 + 7 * step
    out.append(f'<path class="ln"{dash(8, 2)} d="M300 {y8} H318 V{y2} H306"/>' + arrowhead(302, y2, 180))
    y4 = y0 + 3 * step
    out.append(f'<text class="anw" x="50" y="{y4+40}">← exertion, standing or sensory load</text>')
    y7 = y0 + 6 * step
    out.append(f'<path class="lw" style="stroke-width:1.5" d="M8 {y0+4*step-8} V{y7+8}"/>')
    out.append(f'<text class="anw" x="50" y="{y8+44}">5–7: post-exertional malaise</text>')
    return (f'<svg class="{c}" viewBox="0 0 {W} {H+10}" role="group" aria-labelledby="{c}-t" xmlns="http://www.w3.org/2000/svg"><title id="{c}-t">A working model of ME/CFS: eight abnormalities, the last feeding back to the second</title>'
            f'<style>{STYLE_COMMON.format(c=c)}</style>{"".join(out)}</svg>')

# ---------- Figure 2: time course of post-exertional malaise
def tx(h, x0, x1):  # smooth compressive time axis
    return x0 + (x1 - x0) * math.log(1 + h / 24) / math.log(1 + 240 / 24)

def smooth(points):
    d = f"M{points[0][0]:.1f} {points[0][1]:.1f}"
    for i in range(1, len(points)):
        x0, y0 = points[i - 1]; x1, y1 = points[i]; mx = (x0 + x1) / 2
        d += f" C{mx:.1f} {y0:.1f} {mx:.1f} {y1:.1f} {x1:.1f} {y1:.1f}"
    return d

def crash(wide=True):
    c = "fc-w" if wide else "fc-n"
    if wide: W, x0, x1, top, base = 900, 40, 880, 40, 220
    else: W, x0, x1, top, base = 340, 12, 330, 40, 200
    amp = base - top
    Y = lambda v: base - v * amp
    ctrl = [(0, 0), (2, .12), (8, .03), (16, 0), (240, 0)]
    me = [(0, 0), (2, .14), (8, .07), (32, .86), (56, .92), (84, .76), (160, .24), (240, .1)]
    st = STYLE_COMMON.format(c=c) + f"""
.{c} .ax{{stroke:var(--ink);stroke-width:1}} .{c} .tk{{font-size:11.5px;fill:var(--muted);text-anchor:middle;font-variant-numeric:lining-nums}}
.{c} .bl{{stroke:var(--faint);stroke-width:1;stroke-dasharray:2 3}}
.{c} .cc{{fill:none;stroke:var(--cool);stroke-width:2}} .{c} .mc{{fill:none;stroke:var(--warm);stroke-width:2.6}}
.{c} .band{{fill:var(--early);opacity:.13}} .{c} .bt{{font-size:12.5px;font-weight:600;fill:var(--early)}}
.{c} .rug{{stroke:var(--ink);stroke-width:2}} .{c} .rl{{font-size:12px;fill:var(--ink-2)}} .{c} .rh{{font-size:12px;font-weight:600;fill:var(--ink)}}
.{c} .cl{{font-size:13px;font-weight:600}}
"""
    o = []
    bx0, bx1 = tx(48, x0, x1), tx(72, x0, x1)
    o.append(f'<rect class="band" x="{bx0:.1f}" y="{top-12}" width="{bx1-bx0:.1f}" height="{base-top+12}"/>')
    o.append(f'<text class="bt" x="{(bx0+bx1)/2:.1f}" y="{top-18}" text-anchor="middle">rarely sampled</text>' if wide else '')
    o.append(f'<line class="bl" x1="{x0}" y1="{base}" x2="{x1}" y2="{base}"/>')
    o.append(f'<text class="an" x="{x1}" y="{base-6}" text-anchor="end">usual level</text>')
    o.append(f'<line class="ax" x1="{x0}" y1="{base+14}" x2="{x1}" y2="{base+14}"/>')
    ticks = [(0, "0"), (12, "12 h"), (24, "24 h"), (48, "48 h"), (72, "72 h"), (96, "4 d"), (168, "7 d"), (240, "10 d")] if wide else [(0, "0"), (24, "24 h"), (72, "72 h"), (168, "7 d")]
    for h, lab in ticks:
        x = tx(h, x0, x1); o.append(f'<line class="ax" x1="{x:.1f}" y1="{base+14}" x2="{x:.1f}" y2="{base+19}"/><text class="tk" x="{x:.1f}" y="{base+33}">{lab}</text>')
    o.append(f'<text class="an" x="{x0}" y="{top-18}">symptom severity (schematic)</text>' if wide else f'<text class="an" x="{x0}" y="{top-24}">symptom severity, schematic</text>')
    o.append(f'<text class="tk" x="{x0}" y="{base+33}" style="text-anchor:start;opacity:0">.</text>')
    hs = [i * 0.5 for i in range(0, 481)]
    bump = lambda h: 0.13 * (h / 2) * math.exp(1 - h / 2)
    rise = lambda h: 1 / (1 + math.exp(-(h - 20) / 5))
    mef = lambda h: bump(h) + rise(h) * (0.08 + 0.86 / (1 + math.exp((h - 100) / 28)))
    poly = lambda f: "M" + " L".join(f"{tx(h,x0,x1):.1f} {Y(f(h)):.1f}" for h in hs)
    o.append(f'<path class="cc" d="{poly(bump)}"/>')
    o.append(f'<path class="mc" d="{poly(mef)}"/>')
    px = tx(26, x0, x1)
    o.append(f'<text class="cl" x="{px-8:.1f}" y="{Y(.8)-4:.1f}" text-anchor="end" style="fill:var(--warm)">ME/CFS</text>')
    o.append(f'<text class="cl" x="{tx(16,x0,x1)+4:.1f}" y="{base-10}" style="fill:var(--cool)">healthy</text>' if wide else f'<text class="cl" x="{tx(30,x0,x1):.1f}" y="{base-8}" style="fill:var(--cool)">healthy</text>')
    o.append(f'<text class="an" x="{tx(104,x0,x1)+10:.1f}" y="{Y(mef(104))-14:.1f}">recovery over days to weeks</text>' if wide else "")
    # rug of sampling times
    ry = base + 64
    rows = [("Two-day exercise studies, incl. plasma proteome (Germain 2025)" if wide else "Two-day exercise tests, plasma proteome", [0, 24]), ("Muscle biopsy (Appelman 2024, long COVID)", [24]), ("Natural crashes, severe patients", [])]
    o.append(f'<text class="rh" x="{x0}" y="{ry-14}">When studies have sampled</text>')
    for i, (lab, hs) in enumerate(rows):
        y = ry + 10 + i * (30 if wide else 40)
        for h in hs:
            x = tx(h, x0, x1); o.append(f'<line class="rug" x1="{x:.1f}" y1="{y-7}" x2="{x:.1f}" y2="{y+7}"/>')
        if not hs: o.append(f'<text class="bt" x="{x0}" y="{y+5}" style="font-weight:400;font-style:italic">almost never</text>')
        lx = tx(30, x0, x1) if wide else x0
        ly = y + 4 if wide else y + 22
        if not wide and not hs: lx, ly = x0 + 92, y + 5
        o.append(f'<text class="rl" x="{lx:.1f}" y="{ly}">{esc(lab)}</text>')
    o.append(f'<rect class="band" x="{bx0:.1f}" y="{ry}" width="{bx1-bx0:.1f}" height="{(3*(30 if wide else 40))}"/>')
    H = ry + 3 * (30 if wide else 40) + 14
    if not wide:
        o.append(f'<text class="bt" x="{(bx0+bx1)/2:.1f}" y="{H+4}" text-anchor="middle">rarely sampled</text>'); H += 10
    return (f'<svg class="{c}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="{c}-t" xmlns="http://www.w3.org/2000/svg"><title id="{c}-t">Schematic time course of symptoms after one exertion, in a healthy person and in ME/CFS, with the times at which studies have collected samples</title>'
            f'<style>{st}</style>{"".join(o)}</svg>')

for name, svg in [("fig-model-wide.svg", model_wide()), ("fig-model-narrow.svg", model_narrow()), ("fig-crash-wide.svg", crash(True)), ("fig-crash-narrow.svg", crash(False))]:
    open(os.path.join(OUT, name), "w", encoding="utf8").write(svg); print(name, len(svg))
