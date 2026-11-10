# Roofline of the smartsed GPU time loop at 5 m (A40). Usage: roofline.py kern_sum.csv out.svg
# Bytes = minimum DRAM traffic per call (hand-counted from the source); flops = rough hand count.
import csv, sys, math
N, NNZ = 1060123, 5292909; LNNZ = (NNZ + N)//2
BW, PEAK = 551.7, 492.3            # measured: run/perf/stream.cu, run/perf/peak.cu
HS = 8                             # rain gauges (inner loop of computePrecipitation)
# label, name prefix in nsys, bytes per call, flops per call
K = [
 ("IC0 factorisation (csric0)",   "void csric0_lower",                 20*LNNZ+4*N,       10*N),
 ("sediment update",              "computeKernel_sedUpdate",           44*N,              4*N),
 ("SpSV triangular solve",        "void cusparse::spsm_v2_kernel",     12*LNNZ+20*N,      2*LNNZ),
 ("sediment flux x",              "computeKernel_sedFluxHorizontal",   40*N,              3*N),
 ("sediment flux y",              "computeKernel_sedFluxVertical",     40*N,              3*N),
 ("SpMV",                         "void cusparse::csrmv_v3_kernel",    12*NNZ+20*N,       2*NNZ),
 ("axpy",                         "void axpy_kernel_val",              24*N,              2*N),
 ("matrix assembly (gather)",     "buildMatrix_cell_gather",           148*N,             40*N),
 ("precipitation",                "computePrecipitation",              (60+8*HS)*N,       (6+10*HS)*N),
 ("dot",                          "void dot_kernel",                   16*N,              2*N),
 ("friction",                     "computeKernel_friction",            44*N,              20*N),
 ("back-trace x (bilinear)",      "bilinearInterpolationHorizontal(",  32*N,              30*N),
]
rows = list(csv.DictReader(open(sys.argv[1])))
tot = sum(float(r["Total Time (ns)"]) for r in rows)
pts = []
for lab, pre, by, fl in K:
    r = next(r for r in rows if r["Name"].startswith(pre))
    t = float(r["Avg (ns)"])*1e-9; share = 100*float(r["Total Time (ns)"])/tot
    ai = fl/by; p = fl/t/1e9; roof = min(PEAK, BW*ai)
    pts.append((lab, share, t*1e6, by/t/1e9, ai, p, 100*p/roof))
print(f"{'kernel':30s} share%  avg_us  eff_GB/s  AI    GFLOP/s  %roof")
for i, q in enumerate(pts, 1): print(f"{i:2d} {q[0]:27s} {q[1]:5.1f} {q[2]:8.1f} {q[3]:8.0f} {q[4]:6.3f} {q[5]:7.1f} {q[6]:6.0f}")
print("covered share %.1f" % sum(q[1] for q in pts))
# ---- SVG (log-log) ----
W, H, L, T, PW, PH = 1040, 560, 70, 50, 520, 440
x0, x1, y0, y1 = math.log10(0.02), math.log10(4), math.log10(0.5), math.log10(1000)
X = lambda v: L + PW*(math.log10(v)-x0)/(x1-x0); Y = lambda v: T + PH*(1-(math.log10(v)-y0)/(y1-y0))
INK, MUTE, GRID, BLUE = "#1f2328", "#6b7280", "#e5e7eb", "#2563eb"
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" font-family="Helvetica,Arial,sans-serif" font-size="12">',
     f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
     f'<text x="{L}" y="24" font-size="15" font-weight="bold" fill="{INK}">Roofline of the SMART-SED time loop, 5 m grid (1.06 M cells), NVIDIA A40, double precision</text>']
for v in (0.02,0.05,0.1,0.2,0.5,1,2,4):
    o.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{T+PH}" stroke="{GRID}"/><text x="{X(v):.1f}" y="{T+PH+16}" text-anchor="middle" fill="{MUTE}">{v:g}</text>')
for v in (1,10,100,1000):
    o.append(f'<line x1="{L}" x2="{L+PW}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{GRID}"/><text x="{L-8}" y="{Y(v)+4:.1f}" text-anchor="end" fill="{MUTE}">{v:g}</text>')
o.append(f'<text x="{L+PW/2}" y="{T+PH+38}" text-anchor="middle" fill="{INK}">arithmetic intensity [flop / byte]</text>')
o.append(f'<text transform="translate(18,{T+PH/2}) rotate(-90)" text-anchor="middle" fill="{INK}">performance [GFLOP/s]</text>')
ridge = PEAK/BW
o.append(f'<polyline fill="none" stroke="{INK}" stroke-width="2" points="{X(0.02):.1f},{Y(BW*0.02):.1f} {X(ridge):.1f},{Y(PEAK):.1f} {X(4):.1f},{Y(PEAK):.1f}"/>')
o.append(f'<text x="{X(0.03):.1f}" y="{Y(BW*0.03)-10:.1f}" fill="{INK}" transform="rotate(-37 {X(0.03):.1f} {Y(BW*0.03)-10:.1f})">memory roof: {BW:.0f} GB/s (measured)</text>')
o.append(f'<text x="{X(4)-4:.1f}" y="{Y(PEAK)-8:.1f}" text-anchor="end" fill="{INK}">compute roof: {PEAK:.0f} GFLOP/s (measured)</text>')
for i, q in sorted(enumerate(pts, 1), key=lambda z: -z[1][1]):
    r = 4 + 2.6*math.sqrt(q[1]); cx, cy = X(q[4]), Y(q[5])
    o.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{BLUE}" fill-opacity="0.55" stroke="#ffffff" stroke-width="2"><title>{q[0]}: {q[1]:.1f}% of GPU time, {q[6]:.0f}% of roof</title></circle>')
    o.append(f'<text x="{cx:.1f}" y="{cy+4:.1f}" text-anchor="middle" font-size="11" font-weight="bold" fill="{INK}">{i}</text>')
lx = L + PW + 30
o.append(f'<text x="{lx}" y="{T+4}" font-weight="bold" fill="{INK}">#  kernel</text><text x="{lx+270}" y="{T+4}" text-anchor="end" font-weight="bold" fill="{INK}">GPU time</text><text x="{lx+370}" y="{T+4}" text-anchor="end" font-weight="bold" fill="{INK}">% of roof</text>')
for i, q in enumerate(pts, 1):
    y = T + 6 + 24*i
    o.append(f'<text x="{lx}" y="{y}" fill="{INK}">{i}</text><text x="{lx+22}" y="{y}" fill="{INK}">{q[0]}</text><text x="{lx+270}" y="{y}" text-anchor="end" fill="{INK}">{q[1]:.1f} %</text><text x="{lx+370}" y="{y}" text-anchor="end" fill="{INK}">{q[6]:.0f} %</text>')
y = T + 6 + 24*(len(pts)+1) + 10
for k, s in enumerate(["Circle area = share of GPU kernel time (190 loop iterations).",
                       "Time: Nsight Systems. Bytes and flops: counted by hand from",
                       "the source (minimum traffic), so positions are approximate.",
                       "Hardware counters (ncu) are not available to users on Zeus."]):
    o.append(f'<text x="{lx}" y="{y+16*k}" fill="{MUTE}" font-size="11">{s}</text>')
o.append('</svg>'); open(sys.argv[2], "w").write("\n".join(o))
