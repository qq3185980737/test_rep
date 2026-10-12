import re, os, glob
RUN = os.path.dirname(os.path.abspath(__file__))
os.chdir(RUN)
def parse(fn):
    lines=open(fn,encoding='utf-8',errors='ignore').read().splitlines()
    E=[]
    for l in lines:
        if "Total Energy (Potential + Kinetic)" in l:
            try: E.append(float(l.split("=")[1].strip()))
            except: pass
    frames=[]
    for i,l in enumerate(lines):
        if re.search(r"Atom\s+X\(A\)",l):
            blk=[]
            for j in range(i+1,min(i+6,len(lines))):
                m=re.match(r'\s*\d+\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)',lines[j])
                if m: blk.append([float(m.group(1)),float(m.group(2)),float(m.group(3))])
            if len(blk)>=3: frames.append(blk)
    def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    minHO=1e9;bonded=0;maxOO=0
    for f in frames:
        ho=d(f[2],f[1]); minHO=min(minHO,ho)
        if ho<1.1: bonded+=1
        maxOO=max(maxOO,d(f[0],f[1]))
    lastHO=d(frames[-1][2],frames[-1][1])
    dE=(E[-1]-E[0])*1000 if E else None
    return len(E), (E[0] if E else None), dE, minHO, lastHO, round(100.0*bonded/len(frames),1), maxOO
files=sorted(glob.glob("_far_*.out"), key=lambda f:int(re.search(r'_far_(\d+)\.out',f).group(1)))
print(f"{'T(K)':>6} {'n':>4} {'E0(H)':>10} {'dE(mH)':>7} {'minHO':>6} {'lastHO':>6} {'bond%':>6} {'maxOO':>6}")
for f in files:
    T=int(re.search(r'_far_(\d+)\.out',f).group(1))
    n,E0,dE,minHO,lastHO,pct,maxOO=parse(f)
    print(f"{T:>6} {n:>4} {E0:>10.4f} {dE:>7.1f} {minHO:>6.3f} {lastHO:>6.3f} {pct:>6.1f} {maxOO:>6.3f}")
