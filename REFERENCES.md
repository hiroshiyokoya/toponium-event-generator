# References

## This work
- Y. Sumino and H. Yokoya, *Bound-state effects on kinematical distributions of
  top quarks at hadron colliders*, JHEP 09 (2010) 034,
  [arXiv:1007.0075](https://arxiv.org/abs/1007.0075).

## MadGraph / MadEvent
The 2010 code (`legacy/paper-2010/`) is a **MadGraph/MadEvent v4** process
directory (`MGMEVersion.txt`: 4.4.42, `TemplateVersion.txt`: 2.4.21) modified
to include the Green-function (bound-state) corrections. It is *not* based on
the current MadGraph5_aMC@NLO. Please cite the papers of the MG/ME version that
was used:

- J. Alwall, P. Demin, S. de Visscher, R. Frederix, M. Herquet, F. Maltoni,
  T. Plehn, D. L. Rainwater and T. Stelzer, *MadGraph/MadEvent v4: The New Web
  Generation*, JHEP 09 (2007) 028,
  [arXiv:0706.2334](https://arxiv.org/abs/0706.2334).
- F. Maltoni and T. Stelzer, *MadEvent: Automatic Event Generation with
  MadGraph*, JHEP 02 (2003) 027,
  [arXiv:hep-ph/0208156](https://arxiv.org/abs/hep-ph/0208156).
- T. Stelzer and W. F. Long, *Automatic Generation of Tree Level Helicity
  Amplitudes*, Comput. Phys. Commun. 81 (1994) 357,
  [arXiv:hep-ph/9401258](https://arxiv.org/abs/hep-ph/9401258).
- HELAS (helicity amplitude subroutines, `Source/DHELAS`): H. Murayama,
  I. Watanabe and K. Hagiwara, *HELAS: HELicity Amplitude Subroutines for
  Feynman Diagram Evaluations*, KEK-91-11 (1992).

For the current version (not used here):

- J. Alwall, R. Frederix, S. Frixione, V. Hirschi, F. Maltoni, O. Mattelaer,
  H.-S. Shao, T. Stelzer, P. Torrielli and M. Zaro, *The automated computation
  of tree-level and next-to-leading order differential cross sections, and
  their matching to parton shower simulations*, JHEP 07 (2014) 079,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301).

The 2015 extended code (all final states; to be added) uses MadGraph II/HELAS
matrix elements (`me_*.f`) inside a stand-alone Fortran program.

For the current status of MG5_aMC and MadGraph7, see
[docs/madgraph-status.md](docs/madgraph-status.md).

## Other components
- PDFs: CTEQ6L1 (built into the MadEvent v4 `Source/PDF`): J. Pumplin,
  D. R. Stump, J. Huston, H. L. Lai, P. Nadolsky and W. K. Tung, *New
  generation of parton distributions with uncertainties from global QCD
  analysis*, JHEP 07 (2002) 012,
  [arXiv:hep-ph/0201195](https://arxiv.org/abs/hep-ph/0201195).
- Licences of bundled third-party code and data: see
  [THIRD_PARTY.md](THIRD_PARTY.md).
- Companion reference cited in the original `GrnFnc_*/memo.txt`:
  "arXiv:1006.7014". *This ID does not exist* (arXiv returns 404 and INSPIRE
  has no record of it). It is probably a pre-submission number of
  arXiv:1007.0075 itself, but this is not confirmed.
