"""Fraction of LHE events inside invariant-mass windows of the resonances.

Used to compare the 2015 stand-alone results with MadEvent, which keeps only
events with |m(l nu) - m_W| < bwcutoff * Gamma_W for the W's declared as decay
chains in the process card (cut_bw in SubProcesses/myamp.f).

usage: python scripts/lhe_masswindow.py <events.lhe[.gz]> [bwcutoff]
"""
import gzip
import math
import sys

MW, GW = 80.4, 2.0476
MT, GT = 173.0, 1.4911


def events(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        ev = None
        for line in f:
            s = line.strip()
            if s.startswith("<event"):
                ev = []
            elif s.startswith("</event"):
                yield ev[1:]
                ev = None
            elif ev is not None:
                ev.append(s.split())


def mass(ps):
    e = sum(p[3] for p in ps)
    px = sum(p[0] for p in ps)
    py = sum(p[1] for p in ps)
    pz = sum(p[2] for p in ps)
    m2 = e * e - px * px - py * py - pz * pz
    return math.sqrt(max(m2, 0.0))


def main():
    path = sys.argv[1]
    cut = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
    n = n_w = n_t = n_wt = 0
    for ev in events(path):
        parts = [(int(r[0]), int(r[1]), int(r[2]),
                  [float(x) for x in r[6:10]]) for r in ev if len(r) >= 10]
        fin = [(pid, p) for pid, st, m1, p in parts if st == 1]
        lp = [p for pid, p in fin if pid in (-11, -13)]
        lm = [p for pid, p in fin if pid in (11, 13)]
        nu = [p for pid, p in fin if pid in (12, 14)]
        nub = [p for pid, p in fin if pid in (-12, -14)]
        b = [p for pid, p in fin if pid == 5]
        bb = [p for pid, p in fin if pid == -5]
        if not (lp and lm and nu and nub and b and bb):
            continue
        n += 1
        mwp = mass([lp[0], nu[0]])
        mwm = mass([lm[0], nub[0]])
        mtp = mass([lp[0], nu[0], b[0]])
        mtm = mass([lm[0], nub[0], bb[0]])
        inw = abs(mwp - MW) < cut * GW and abs(mwm - MW) < cut * GW
        int_ = abs(mtp - MT) < cut * GT and abs(mtm - MT) < cut * GT
        n_w += inw
        n_t += int_
        n_wt += inw and int_

    def frac(k):
        f = k / n
        return f, math.sqrt(f * (1 - f) / n)

    print(f"events: {n}   window: {cut} x Gamma")
    for name, k in (("both W  in window", n_w), ("both t  in window", n_t),
                    ("W and t in window", n_wt)):
        f, e = frac(k)
        print(f"  {name}: {f:.4f} +- {e:.4f}")


if __name__ == "__main__":
    main()
