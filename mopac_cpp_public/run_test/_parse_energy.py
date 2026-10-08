import re
def parse(path):
    E=[]
    cur=None
    for ln in open(path, encoding="utf-8", errors="replace"):
        m=re.search(r"MD Step\s+(\d+)", ln)
        if m: cur=int(m.group(1))
        m=re.search(r"Total Energy.*=\s*(-?[\d.]+)", ln)
        if m and cur is not None: E.append((cur,float(m.group(1))))
    return E
for path in [r"D:\test_rep\mopac_cpp_public\run_test\_md_h2o2_6atom.out",
             r"D:\test_rep\mopac_cpp_public\run_test\_md_h2o2_50mol.out"]:
    E=parse(path)
    e0=E[0][1]
    n=len(E); stride=max(1,n//60)
    sel=[(s,(e-e0)*1000.0) for i,(s,e) in enumerate(E) if i%stride==0]
    print(path.split("\\")[-1], "E0=%.6f"%e0, "n=%d"%n)
    print("dE_series=", sel)
