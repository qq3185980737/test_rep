import subprocess, re, os, time, sys
EXE = r"D:\test_rep\mopac_cpp_public\x64\Release\mopac.exe"
RUN = os.path.dirname(os.path.abspath(__file__))
os.chdir(RUN)
# Far approach: H3 starts ~2.6 A from O2, must cross the barrier by thermal energy
GEOM = """O  5.00000000  5.00000000  5.00000000
O  5.00000000  5.00000000  6.20000000
H  7.50000000  5.00000000  6.80000000
"""
TMP = [1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]

def run_one(T):
    name = f"_far_{T}"
    mop = f"""PM7 MOLDYN TEMP={T} BOX=12 MDSTEPS=400 DT=0.2
Mixture: far-approach scan {T}K
Constraint:
{GEOM}"""
    open(name+".mop","w").write(mop)
    if os.path.exists(name+".out"): os.remove(name+".out")
    t0=time.time()
    try:
        rc = subprocess.run([EXE, name+".mop"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300).returncode
    except subprocess.TimeoutExpired:
        rc = "TIMEOUT"
    dt = time.time()-t0
    out = name+".out"
    if not os.path.exists(out) or rc!=0:
        return dict(T=T, rc=rc, time=round(dt,1), E0=None, minHO=None, lastHO=None, pct=None, maxOO=None, note="no-out")
    lines = open(out, encoding='utf-8', errors='ignore').read().splitlines()
    E=[]
    for l in lines:
        if "Total Energy (Potential + Kinetic)" in l:
            try: E.append(float(l.split("=")[1].strip()))
            except: pass
    frames=[]
    for i,l in enumerate(lines):
        if re.search(r"Atom\s+X\(A\)", l):
            blk=[]
            for j in range(i+1, min(i+6,len(lines))):
                m=re.match(r'\s*\d+\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)', lines[j])
                if m: blk.append([float(m.group(1)),float(m.group(2)),float(m.group(3))])
            if len(blk)>=3: frames.append(blk)
    def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    minHO=1e9; bonded=0; maxOO=0
    for f in frames:
        ho=d(f[2],f[1]); minHO=min(minHO,ho)
        if ho<1.1: bonded+=1
        maxOO=max(maxOO, d(f[0],f[1]))
    lastHO=d(frames[-1][2],frames[-1][1])
    return dict(T=T, rc=rc, time=round(dt,1), E0=E[0] if E else None,
                minHO=round(minHO,3), lastHO=round(lastHO,3),
                pct=round(100.0*bonded/len(frames),1), maxOO=round(maxOO,3))

rows=[run_one(T) for T in TMP]
for r in rows: print(r)
print("\n=== SUMMARY (far approach, 400 steps) ===")
print(f"{'T(K)':>6} {'rc':>5} {'t(s)':>5} {'E0(H)':>10} {'minHO':>6} {'lastHO':>6} {'bond%':>6} {'maxOO':>6}")
for r in rows:
    E0s=f"{r['E0']:.4f}" if r['E0'] is not None else "-"
    print(f"{r['T']:>6} {str(r['rc']):>5} {r['time']:>5} {E0s:>10} {str(r['minHO']):>6} {str(r['lastHO']):>6} {str(r['pct']):>6} {str(r['maxOO']):>6}")
