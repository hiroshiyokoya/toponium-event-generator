"""Phase-space points for the matrix-element comparison (issue #9).

Takes gg-initiated events from a 2015 blvblv LHE file and writes the momenta
in the MadEvent (2010) order:  g g mu+ vm mu- vm~ b b~
one event per block of 8 lines "E px py pz".

usage: python make_points.py <events.lhe[.gz]> <out.dat> [n] [idprup]
"""
import gzip
import sys

ORDER = (-13, 14, 13, -14, 5, -5)


def main():
    src, out = sys.argv[1], sys.argv[2]
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    idprup = int(sys.argv[4]) if len(sys.argv) > 4 else None
    op = gzip.open if src.endswith(".gz") else open
    n = 0
    with op(src, "rt") as f, open(out, "w") as g:
        ev = None
        for line in f:
            t = line.strip()
            if t.startswith("<event"):
                ev = []
                continue
            if t.startswith("</event"):
                rows = [r.split() for r in ev[1:]]
                ins = [r for r in rows if r[1] == "-1"]
                fin = {int(r[0]): r for r in rows if r[1] == "1"}
                okp = idprup is None or int(ev[0].split()[1]) == idprup
                if okp and all(r[0] == "21" for r in ins) and all(k in fin for k in ORDER):
                    for r in ins + [fin[k] for k in ORDER]:
                        px, py, pz, e = (float(x) for x in r[6:10])
                        g.write(f"{e:.10e} {px:.10e} {py:.10e} {pz:.10e}\n")
                    n += 1
                    if n >= nmax:
                        break
                ev = None
                continue
            if ev is not None:
                ev.append(t)
    print(f"wrote {n} points to {out}")


if __name__ == "__main__":
    main()
