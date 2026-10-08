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

Not shipped here (to be handled in #5): BASES/SPRING (KEK), used by the 2015
stand-alone version. Its redistribution terms have not been checked yet.

If you are an author of any of the above and want a change, please open an
issue.
