import re, sys
fn = sys.argv[1] if len(sys.argv)>1 else "_md_h_o2_hi.out"
lines = open(fn, encoding='utf-8', errors='ignore').read().splitlines()
frames=[]
for i,l in enumerate(lines):
    if re.search(r"Atom\s+X\(A\)", l):
        blk=[]
        for j in range(i+1, min(i+6,len(lines))):
            m=re.match(r'\s*(\d+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)', lines[j])
            if m: blk.append([float(m.group(2)),float(m.group(3)),float(m.group(4))])
        if len(blk)>=3: frames.append(blk)
def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
ho=[d(f[2],f[1]) for f in frames]
oo=[d(f[0],f[1]) for f in frames]
print("H3O2=" + ",".join(f"{x:.3f}" for x in ho))
print("O1O2=" + ",".join(f"{x:.3f}" for x in oo))
