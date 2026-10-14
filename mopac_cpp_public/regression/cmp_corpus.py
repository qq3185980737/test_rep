# -*- coding: utf-8 -*-
"""Compare FINAL HEAT values between the recorded C++ outputs (out_cpp/) and
the Fortran .arc summaries (out_ftn/) of the released regression corpus.

Run from anywhere: the directories are resolved relative to this script.
Expected result on the archived corpus: 327 paired inputs, 320 bitwise
identical, 7 EF-optimization value mismatches (documented 1-ULP build
sensitivity), 22 inputs without a Fortran .arc counterpart.
"""
import re
import os
import glob
import json

_HERE = os.path.dirname(os.path.abspath(__file__))
CXX_DIR = os.path.join(_HERE, "out_cpp")
FTN_DIR = os.path.join(_HERE, "out_ftn")

HEAT_RE_OUT = re.compile(r"final heat of formation\s*=\s*([-+\d.]+)", re.IGNORECASE)
HEAT_RE_ARC = re.compile(r"heat of formation\s*=\s*([-+\d.]+)", re.IGNORECASE)


def last_heat(text, regex):
    hits = regex.findall(text)
    return hits[-1] if hits else None


def collect(dirpath, regex, ext):
    out = {}
    for p in glob.glob(os.path.join(dirpath, "*" + ext)):
        base = os.path.splitext(os.path.basename(p))[0]
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                h = last_heat(fh.read(), regex)
        except Exception:
            h = None
        out[base] = h
    return out


cxx = collect(CXX_DIR, HEAT_RE_OUT, ".out")
ftn = collect(FTN_DIR, HEAT_RE_ARC, ".arc")

names = sorted(set(cxx) & set(ftn))
missing_cxx = sorted(set(ftn) - set(cxx))
missing_ftn = sorted(set(cxx) - set(ftn))

identical = []
diff = []
for n in names:
    a, b = cxx[n], ftn[n]
    if a is None or b is None:
        diff.append((n, a, b, "missing_heat"))
    elif a == b:
        identical.append(n)
    else:
        diff.append((n, a, b, "value_mismatch"))

report = {
    "total_paired": len(names),
    "bitwise_identical": len(identical),
    "diff": len(diff),
    "missing_cxx_output": len(missing_cxx),
    "missing_ftn_output": len(missing_ftn),
    "diff_list": diff[:40],
    "missing_cxx_head": missing_cxx[:20],
    "missing_ftn_head": missing_ftn[:20],
}

print(json.dumps(report, ensure_ascii=False, indent=1))
