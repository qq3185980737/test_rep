import random, math
random.seed(7)
# 25 H2 + 25 O2 = 50 molecules, 75 atoms, in a 20 A box
out = []
out.append("PM7 TEMP=2000 BOX=20 MDSTEPS=60")
out.append("Mixture: 25 H2 + 25 O2 (random thermal)")
out.append("Constraint: ")
cells = [(i, j, k) for i in range(5) for j in range(5) for k in range(2)]
mols = []
for n in range(50):
    ci, cj, ck = cells[n]
    cx = 2 + ci * 4.0
    cy = 2 + cj * 4.0
    cz = 2 + ck * 8.0
    u = random.random(); v = random.random()
    th = 2 * math.pi * u; ph = math.acos(2 * v - 1)
    dx = math.sin(ph) * math.cos(th); dy = math.sin(ph) * math.sin(th); dz = math.cos(ph)
    bond = 0.74 if n < 25 else 1.21
    el = "H" if n < 25 else "O"
    mols.append((el,
                 cx - 0.5 * bond * dx, cy - 0.5 * bond * dy, cz - 0.5 * bond * dz,
                 cx + 0.5 * bond * dx, cy + 0.5 * bond * dy, cz + 0.5 * bond * dz))
for (el, x1, y1, z1, x2, y2, z2) in mols:
    out.append("%s %12.6f %12.6f %12.6f" % (el, x1, y1, z1))
    out.append("%s %12.6f %12.6f %12.6f" % (el, x2, y2, z2))
open(r"D:\test_rep\mopac_cpp_public\run_test\_md_h2o2_50mol.mop", "w").write("\n".join(out) + "\n")
print("atoms:", 2 * len(mols))
print("molecules:", len(mols))
print("total lines:", len(out))
