"""Compare Green-function tables row by row (docs/REVIEW.md).

usage:
  python compare_tables.py free <table> <ie_step>            # |Gamma - 1| (alpha_s = 0)
  python compare_tables.py diff <table_a> <table_b> <ie_step> # |a - b| on computed rows

Rows (energy index IE) not computed in a partial run are skipped by taking
every <ie_step>-th row.
"""
import sys


def read(path):
    with open(path) as f:
        head = [f.readline().split() for _ in range(4)]
        ne, np_ = int(head[1][0]), int(head[1][1])
        emin, emax = float(head[2][0]), float(head[2][1])
        rows, cur = [], []
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == "---":
                rows.append(cur)
                cur = []
            else:
                cur.append(complex(float(t[0]), float(t[1])))
    assert len(rows) == ne + 1 and all(len(r) == np_ + 1 for r in rows)
    es = [emin + i * (emax - emin) / ne for i in range(ne + 1)]
    return head, es, rows


def main():
    mode = sys.argv[1]
    if mode == "free":
        _, es, rows = read(sys.argv[2])
        step = int(sys.argv[3])
        for ie in range(0, len(rows), step):
            dev = max(abs(g - 1) for g in rows[ie])
            print(f"E = {es[ie]:8.2f}  max_p |Gamma-1| = {dev:.2e}")
    else:
        ha, es, ra = read(sys.argv[2])
        hb, _, rb = read(sys.argv[3])
        step = int(sys.argv[4])
        print("headers:", " ".join(ha[0]), "|", " ".join(hb[0]))
        for ie in range(0, len(ra), step):
            d = max(abs(a - b) for a, b in zip(ra[ie], rb[ie]))
            m = max(abs(a) for a in ra[ie])
            print(f"E = {es[ie]:8.2f}  max_p |a-b| = {d:.2e}  (max |a| = {m:.3f})")


if __name__ == "__main__":
    main()
