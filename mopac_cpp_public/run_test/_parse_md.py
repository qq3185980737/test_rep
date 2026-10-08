import re
def parse(path):
    steps=[]; T=[]; E=[]; Ek=[]
    cur=None
    for ln in open(path, encoding="utf-8", errors="replace"):
        m=re.search(r"MD Step\s+(\d+)", ln)
        if m: cur=int(m.group(1))
        m=re.search(r"Total Energy.*=\s*(-?[\d.]+)", ln)
        if m and cur is not None: E.append((cur,float(m.group(1))))
        m=re.search(r"Instant Temperature.*=\s*([\d.]+)", ln)
        if m and cur is not None: T.append((cur,float(m.group(1))))
        m=re.search(r"Total Kinetic Energy.*=\s*([\d.]+)", ln)
        if m and cur is not None: Ek.append((cur,float(m.group(1))))
    return T,E,Ek
for path in [r"D:\test_rep\mopac_cpp_public\run_test\_md_h2o2_6atom.out",
             r"D:\test_rep\mopac_cpp_public\run_test\_md_h2o2_50mol.out"]:
    T,E,Ek=parse(path)
    print("FILE", path.split("\\")[-1])
    print("steps:", len(T))
    # downsample to <=60 points for plotting
    n=len(T); stride=max(1, n//60)
    sel_T=[t for i,t in enumerate(T) if i%stride==0]
    sel_E=[e for i,e in enumerate(E) if i%stride==0]
    print("T_series=[", ",".join("%d:%.1f"%(s,t) for s,t in sel_T), "]")
    print("E_min=%.4f E_max=%.4f"%(min(e[1] for e in E), max(e[1] for e in E)))
    print("E_first=%.4f E_last=%.4f"%(E[0][1],E[-1][1]))
    print("T_first=%.1f T_last=%.1f"%(T[0][1],T[-1][1]))
