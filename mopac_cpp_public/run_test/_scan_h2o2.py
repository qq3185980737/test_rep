import subprocess, re, os, time
EXE = r"D:\test_rep\mopac_cpp_public\x64\Release\mopac.exe"
RUN = os.path.dirname(os.path.abspath(__file__))
os.chdir(RUN)
# True H2 + O2: H2 molecule attacks O2 non-collinearly from +x side
GEOM = """O  5.00000000  5.00000000  5.00000000
O  5.00000000  5.00000000  6.20000000
H  7.20000000  5.00000000  6.50000000
H  7.90000000  5.00000000  6.90000000
"""
TMP = [2000,2500,3000,3500,4000,4500,5000,6000,7000]

def run_one(T):
    name = f"_h2o2_{T}"
    mop = f"""PM7 MOLDYN TEMP={T} BOX=12 MDSTEPS=500 DT=0.2
Mixture: H2 + O2 scan {T}K
Constraint:
{GEOM}"""
    open(name+".mop","w").write(mop)
    if os.path.exists(name+".out"): os.remove(name+".out")
    try:
        rc = subprocess.run([EXE, name+".mop"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300).returncode
    except subprocess.TimeoutExpired:
        rc = "TIMEOUT"
    if not os.path.exists(name+".out") or rc!=0:
        return dict(T=T, rc=rc, E0=None, minHO=None, bondpct=None, HHmax=None, OOmax=None, last=None)
    lines = open(name+".out", encoding='utf-8', errors='ignore').read().splitlines()
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
            if len(blk)>=4: frames.append(blk)
    def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    minHO=1e9; bond=0; HHmax=0; OOmax=0
    for f in frames:
        o1,o2,h3,h4=f[0],f[1],f[2],f[3]
        for h in (h3,h4):
            for o in (o1,o2):
                hh=d(h,o)
                minHO=min(minHO,hh)
                if hh<1.1: bond+=1
        HHmax=max(HHmax, d(h3,h4))
        OOmax=max(OOmax, d(o1,o2))
    last=frames[-1]
    lHO=min(d(last[2],last[0]),d(last[2],last[1]),d(last[3],last[0]),d(last[3],last[1]))
    lHH=d(last[2],last[3]); lOO=d(last[0],last[1])
    return dict(T=T, rc=rc, E0=E[0] if E else None, minHO=round(minHO,3),
                bondpct=round(100.0*bond/len(frames),1), HHmax=round(HHmax,3),
                OOmax=round(OOmax,3), last=f"HO={lHO:.2f} HH={lHH:.2f} OO={lOO:.2f}")

rows=[run_one(T) for T in TMP]
out=[f"{'T':>6} {'rc':>6} {'E0':>10} {'minHO':>6} {'bond%':>6} {'HHmax':>6} {'OOmax':>6} {'last(HO/HH/OO)':>22}"]
for r in rows:
    E0=f"{r['E0']:.4f}" if r['E0'] is not None else "-"
    out.append(f"{r['T']:>6} {str(r['rc']):>6} {E0:>10} {str(r['minHO']):>6} {str(r['bondpct']):>6} {str(r['HHmax']):>6} {str(r['OOmax']):>6} {str(r['last']):>22}")
open("_h2o2_summary.txt","w").write("\n".join(out))
print("\n".join(out))
