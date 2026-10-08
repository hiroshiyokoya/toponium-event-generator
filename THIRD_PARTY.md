# Third-party code and data

The MIT licence in `LICENSE` covers the code written by the authors of this
repository: the Green-function programs (`legacy/paper-2010/GrnEP_pro/`), the
Green-function modifications and tables (`legacy/paper-2010/*/GrnFnc_*/`), and
the scripts and docs.

The rest of `legacy/paper-2010/pp*_Grn_v42/` is a MadGraph/MadEvent v4 process
directory. It contains code and data by other authors, shipped unchanged so that
the 2010 results can be reproduced. That material is **not** covered by the MIT
licence. Status as of 2026-10-08:

| Path (in each `pp*_Grn_v42/`) | Origin | Licence status |
|---|---|---|
| `Source/`, `SubProcesses/`, `bin/`, `Cards/`, `HTML/`, `lib/` (except below) | MadGraph/MadEvent v4.4.42 (MG/ME team: F. Maltoni, T. Stelzer et al.) | No licence file in the v4 process directory. The successor MG5_aMC is released under an adapted University of Illinois/NCSA licence (redistribution allowed, notices kept, no endorsement). We keep all original headers. **To be confirmed** with the MadGraph team for v4. |
| `Source/DHELAS/` | HELAS (H. Murayama, I. Watanabe, K. Hagiwara, KEK-91-11), extended by the MG team | No explicit licence. Distributed with every MadGraph release. |
| `Source/CERNLIB/` (`abend.f`, `dlsqp2.f`, `lenocc.f`, `mtlprt.f`, `mtlset.f`, `radmul.f`) | CERN Program Library, as bundled by MadGraph | CERNLIB is distributed under the GNU GPL. These files stay under the GPL. |
| `Source/PDF/` (`Ctq*.f`, `cteq3.f`, `Partonx5.f`, `jeppe02.f`, ...), `lib/Pdfdata/` (`cteq*.tbl`, `mrs*.dat`, `mrst2002nlo.dat`) | CTEQ and MRST parton distribution codes and grids, as bundled by MadGraph | Publicly distributed by the PDF groups, also via LHAPDF. No explicit licence. Please cite the PDF papers. |
| `Source/MadWeight_File/Python/*.py` (13 files: `accep.py`, `Cards.py`, `clean.py`, `cluster.py`, ...) | MadWeight (matrix-element reweighting tool bundled with MG/ME v4; not used by this generator) | Header says "license: GNU" (GNU GPL, version not stated). These files stay under the GPL. |
| `Source/MadWeight_File/Python/progressbar.py` | Python progress-bar module bundled with MadWeight | GNU LGPL, version 2.1 or later (stated in the file). |

## Not shipped: BASES/SPRING

The stand-alone version (`standalone/`) uses BASES/SPRING V5.1 (S. Kawabata,
Comput. Phys. Commun. 88 (1995) 309). It is distributed by the CPC Program
Library on Mendeley Data under the **CPC licence**. That licence allows
academic and non-profit use, requires citing the paper, and does not allow
passing the code to third parties. It is therefore not included in this
repository. Users download it themselves and put the zip into `third_party/`.
`scripts/extract_bases51.sh` extracts it at build time (see README).

`standalone/bases51/bases51_compat.f` was written for this repository and is
under the MIT licence. It re-implements `BSSETD` and `BSSETP` (KEK V5.0 API),
`XHSAVE2` (the author's add-on) and the CERNLIB `DATIME`, `UCOPY`, `TIMEX`,
`TIMEST`. It contains no BASES/SPRING code.

`standalone/common/intvegas.f` is a VEGAS variant (G. P. Lepage), used only
when `IBSS=0`. Its origin and licence are still to be checked.

If you are an author of any of the above and want a change, please open an
issue.
