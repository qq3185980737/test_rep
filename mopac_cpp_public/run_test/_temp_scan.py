import subprocess, re, os, sys, time
EXE = r"D:\test_rep\mopac_cpp_public\x64\Release\mopac.exe"
RUN = os.path.dirname(os.path.abspath(__file__))
os.chdir(RUN)
GEOM = """O  5.00000000  5.00000000  5.00000000
O  5.00000000  5.00000000  6.20000000
H  6.20000000  5.00000000  6.60000000
"""
TMP = [2500,3000,3250,3500,3750,4000,4500,5000]

def run_one(T):
    name = f"_scan_{T}"
    mop = f"""PM7 MOLDYN TEMP={T} BOX=10 MDSTEPS=400 DT=0.2
Mixture: temperature scan {T}K, H+O2->OH+O
Constraint:
{GEOM}"""
    open(name+".mop","w").write(mop)
    for ext in (".out",):
        if os.path.exists(name+ext): os.remove(name+ext)
    t0=time.time()
    try:
        r = subprocess.run([EXE, name+".mop"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
        rc = r.returncode
    except subprocess.TimeoutExpired:
        rc = "TIMEOUT"
    dt = time.time()-t0
    out = name+".out"
    if not os.path.exists(out) or rc != 0:
        return dict(T=T, rc=rc, time=round(dt,1), n=0, E0=None, Elast=None,
                    minHO=None, lastHO=None, bonded=None, maxOO=None, note="no-out")
    lines = open(out, encoding='utf-8', errors='ignore').read().splitlines()
    # energy
    E=[]
    for l in lines:
        if "Total Energy (Potential + Kinetic)" in l:
            try: E.append(float(l.split("=")[1].strip()))
            except: pass
    # frames
    frames=[]
    for i,l in enumerate(lines):
        if re.search(r"Atom\s+X\(A\)", l):
            blk=[]
            for j in range(i+1, min(i+6,len(lines))):
                m=re.match(r'\s*\d+\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)', lines[j])
                if m: blk.append([float(m.group(1)),float(m.group(2)),float(m.group(3))])
            if len(blk)>=3: frames.append(blk)
    def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    if not frames:
        return dict(T=T, rc=rc, time=round(dt,1), n=len(E), E0=E[0] if E else None,
                    Elast=E[-1] if E else None, note="no-frame")
    minHO=1e9; bonded=0; maxOO=0
    for f in frames:
        ho=d(f[2],f[1])
        minHO=min(minHO,ho)
        if ho<1.1: bonded+=1
        maxOO=max(maxOO, d(f[0],f[1]))
    last=frames[-1]
    lastHO=d(last[2],last[1])
    return dict(T=T, rc=rc, time=round(dt,1), n=len(E),
                E0=E[0] if E else None, Elast=E[-1] if E else None,
                minHO=round(minHO,3), lastHO=round(lastHO,3),
                bonded=bonded, tot=len(frames), pct=round(100.0*bonded/len(frames),1),
                maxOO=round(maxOO,3))

rows=[]
for T in TMP:
    rows.append(run_one(T))
    print(f"done T={T}: {rows[-1]}")
    sys.stdout.flush()

print("\n=== SUMMARY ===")
print(f"{'T(K)':>6} {'rc':>6} {'t(s)':>6} {'E0(H)':>10} {'dE(mH)':>7} {'minHO':>6} {'lastHO':>6} {'bond%':>6} {'maxOO':>6}")
for r in rows:
    dE = (r['Elast']-r['E0'])*1000 if r['E0'] is not None and r['Elast'] is not None and r['n'] else None
    dEs = f"{dE:.1f}" if dE is not None else "-"
    E0s = f"{r['E0']:.4f}" if r['E0'] is not None else "-"
    print(f"{r['T']:>6} {str(r['rc']):>6} {r['time']:>6} {E0s:>10} {dEs:>7} {str(r['minHO']):>6} {str(r['lastHO']):>6} {r.get('pct','-'):>6} {str(r['maxOO']):>6}")
