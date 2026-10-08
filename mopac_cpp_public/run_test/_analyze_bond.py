import re, sys
fn = sys.argv[1] if len(sys.argv)>1 else "_md_h_o2_hi.out"
lines = open(fn, encoding='utf-8', errors='ignore').read().splitlines()
frames=[]; 
for i,l in enumerate(lines):
    if re.search(r"Atom\s+X\(A\)", l):
        blk=[]
        for j in range(i+1, min(i+6,len(lines))):
            m=re.match(r'\s*(\d+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)', lines[j])
            if m: blk.append([float(m.group(2)),float(m.group(3)),float(m.group(4))])
        if len(blk)>=3: frames.append(blk)
def d(a,b): return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
minHO=1e9; minOO=1e9; maxOO=0; minHO_at=0
for k,f in enumerate(frames):
    h3,o1,o2=f[2],f[0],f[1]
    ho=d(h3,o2); oo=d(o1,o2)
    if ho<minHO: minHO=ho; minHO_at=k
    minOO=min(minOO,oo); maxOO=max(maxOO,oo)
print(f"frames={len(frames)} min H3-O2={minHO:.3f} A (frame {minHO_at})  O1-O2 range min={minOO:.3f} max={maxOO:.3f}")
last=frames[-1]
h3,o1,o2=last[2],last[0],last[1]
print(f"LAST: O1={o1} O2={o2} H3={h3}")
print(f"LAST H3-O2={d(h3,o2):.3f} A  O1-O2={d(o1,o2):.3f} A")
# how long H3 stayed near O (bonded, <1.1 A)
near=sum(1 for k,f in enumerate(frames) if d(f[2],f[1])<1.1)
print(f"frames with H3-O2<1.1A (bonded): {near}/{len(frames)}")
