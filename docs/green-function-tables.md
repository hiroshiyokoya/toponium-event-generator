# Green-function tables: how they are made and how the generators use them

The bound-state (toponium) effects enter the event generators through one complex
number per phase-space point, the **correction factor**

$$
\Gamma^{(c)}(E,p) \;=\; \frac{G^{(c)}(E,p)}{G_0(E,p)}, \qquad
G_0(E,p) = -\frac{1}{E + i\Gamma_t - p^2/m_t},
$$

where $G^{(c)}$ is the momentum-space Green function of the $t\bar t$ pair in colour
state $c$ (singlet or octet), and $G_0$ is the free one. This is Eq. (3.1) of
arXiv:1007.0075. Computing $G^{(c)}$ takes time, so it is computed once and tabulated as
a function of the energy $E$ and of the momentum $p = |\vec p_t|$ in the $t\bar t$ rest
frame. The generators then read the table and interpolate.

## 1. The tables

| File | Colour | Region |
|---|---|---|
| `grnep1thre.tbl` | singlet | threshold |
| `grnep1high.tbl` | singlet | high energy |
| `grnep8thre.tbl` | octet | threshold |
| `grnep8high.tbl` | octet | high energy |

The repository ships one set, made for $m_t = 173$ GeV, $\Gamma_t = 1.4911$ GeV and
$\alpha_s(\mu_B) = 0.1534$ at $\mu_B = 20$ GeV, with the one-loop fixed-order QCD
potential. It is in `legacy/paper-2010/ppblvblv_Grn_v42/GrnFnc_blvblv/`, and the identical
set is in `ppblvbjj_Grn_v42/GrnFnc_blvbjj/` and `ppbwbw_Grn_v42/GrnFnc_bwbw/`. Both versions of
the generator use it.

### Format

Plain text, written with list-directed `WRITE`:

```
  173.  1.4911  0.1534        <- m_t, Gamma_t, alpha_s(mu_B)
 800 300                      <- NE, NP  (number of intervals)
 -50.  30.                    <- E_min, E_max [GeV]
  0.  300.                    <- threshold table: p_min, p_max [GeV]
                                 high-energy table: RR (see below)
 Re Gamma  Im Gamma           <- (NE+1) blocks of (NP+1) lines,
 ...                             E index outer, p index inner,
 --- 0 ---                       each block closed by a separator line
 ...
```

### Grids

| | Threshold (`*thre.tbl`) | High energy (`*high.tbl`) |
|---|---|---|
| $E$ range | −50 … 30 GeV | 30 … 830 GeV |
| $E$ grid | 801 points, uniform (0.1 GeV) | 401 points, uniform (2 GeV) |
| $p$ range | 0 … 300 GeV | around the on-shell momentum $p_0 = \sqrt{m_t E}$ |
| $p$ grid | 301 points, uniform (1 GeV) | 102 points, geometric: dense near $p_0$, ratio `RR` = 0.9 on each side |

The grid constants are in `tblgrnep.inc` (`NE1, NP1, EMIN1, …, NE2, NP2, EMIN2, EMAX2, RR`).
They must be the same in the program that writes the table and in the one that reads it.

## 2. Making tables (`legacy/paper-2010/GrnEP_pro/`)

### What the programs do

`GrnEPthre.f` (threshold grid) and `GrnEPhigh.f` (high-energy grid) do the following for
each colour state (`ICLR` = 1 singlet, 2 octet) and each $E$ on the grid:

1. Solve the S-wave Schrödinger equation for the Green function in coordinate space,
   $\left[-\nabla^2/m_t + V_{\rm QCD}(r) - (E + i\Gamma_t)\right] G(\vec r) = \delta^3(\vec r)$.
   This uses a Runge–Kutta method in $\log r$ (`SOLVE`, `REFINE`).
2. Fourier-transform it to momentum space with Simpson's rule (`FOURIER`), for all $p$ on
   the grid, and multiply by $-(E + i\Gamma_t - p^2/m_t)$ to get $\Gamma(E,p)$.
3. Write the table (`SAVETABLE`).

As a by-product, $\mathrm{Im}\,G(E, r=0)$ is written to `fort.31` … `fort.34`, and the
potential to `fort.1`/`fort.2`.

### Inputs

They are set at the top of the main programs of **both** `GrnEPthre.f` and `GrnEPhigh.f`.
Keep them identical in the two files.

| Quantity | Where | Default |
|---|---|---|
| $m_t$ | `MU = 173D0` | 173 GeV |
| $\Gamma_t$ | `GAMT = 1.4911D0` (or `RGAMT(...)`, the LO width from $G_F, m_W, m_b, V_{tb}$, `gamt.f`) | 1.4911 GeV |
| $\mu_B$, $\alpha_s(\mu_B)$ | `MUB = 20D0`, `ASB = ALPHAS(MUB, MZ, ASZ, NF, NORD)` (QCD-PEGASUS, $\alpha_s(M_Z) = 0.118$, NLO running) | 20 GeV, 0.15329 |
| QCD potential | `Vqcd.f` (copy one of the files below to it) | `Vqcd_1LFO.f` |

The shipped tables carry $\alpha_s(\mu_B) = 0.1534$ in their header, not the 0.15329 from
QCD-PEGASUS. They were probably made with the alternative line in the main program,
`ASB = 0.1534D0`, which is commented out.

| Potential | File | Remarks |
|---|---|---|
| one-loop, fixed order | `Vqcd_1LFO.f` | default; singlet ($C = -C_F$) and octet ($C = 1/6$) |
| two-loop, fixed order | `Vqcd_2LFO.f` | M. Peter, Nucl. Phys. B 501 (1997) 471, and Kniehl et al., Phys. Lett. B 607 (2005) 96 |
| two-loop, RG-improved | `Vqcd_2LRGI.f` | Sumino et al., Phys. Rev. D 47 (1993) 56 (`POTQCD.f`, `EI.f`). Singlet only. Needs the initialisation `CALL SETVPR(...)`, which is commented out in the main programs; without it the program crashes |

### Running (Docker)

Work in a copy, because the programs write their tables and `fort.*` logs into the current
directory:

```bash
docker run --rm -it -v "$PWD:/work" toponium-eg:dev
# inside the container
cp -r /work/legacy/paper-2010/GrnEP_pro /tmp/grn && cd /tmp/grn
cp Vqcd_1LFO.f Vqcd.f              # choose the potential
# edit MU, GAMT, MUB/ASB at the top of GrnEPthre.f and GrnEPhigh.f
make LD=gfortran FFLAGS="-O3 -ffixed-line-length-132 -std=legacy -w"
./GrnEPthre.exe > thre.log         # writes grnep1thre.tbl, grnep8thre.tbl
./GrnEPhigh.exe > high.log         # writes grnep1high.tbl, grnep8high.tbl
```

The original `Makefile` names `g77`, so `LD=gfortran` is needed. All three potentials
compile with gfortran 13. With gfortran 13 the threshold program does about 150 energy
points per minute (about 11 minutes for both colours), and the high-energy program about
30 (about 26 minutes). The four tables thus take about 40 minutes. The 2010 memo spoke of
several hours per program.

`plotgrnep.f` (`make plotgrnep.exe`) is a small example that reads a table. It writes
$|\Gamma|^2$ and $|G|^2$ along $p$ at fixed $E$ to `fort.7` and `fort.8`.

## 3. How the generators use the tables

### Kinematics

For each event, `GRNMOD` computes

$$
E = \frac{\hat s}{4 m_t} - m_t, \qquad
p = \frac{\sqrt{\hat s}}{2}\,\lambda^{1/2}\!\left(1, \frac{s_t}{\hat s}, \frac{s_{\bar t}}{\hat s}\right),
$$

from the partonic $\hat s$ and the invariant masses $s_t = (p_b + p_{W^+})^2$ and
$s_{\bar t} = (p_{\bar b} + p_{W^-})^2$. $E \simeq \sqrt{\hat s} - 2m_t$ near threshold,
and $p$ is the top momentum in the $t\bar t$ rest frame.

It then reads $\Gamma^{(c)}(E,p)$ by bilinear interpolation on the grid (`READGRN` →
`HOKAN2D`). The tables are read once per colour, on the first call.

### Where the factor enters

`GGTT` (2010) / `GGTTRES` (2015) multiply the resonant $gg \to t\bar t$ (or $q\bar q \to t\bar t$)
amplitudes. The switches are `IGRN`, `INR` and `ITR`; see
[versions-2010-2015.md](versions-2010-2015.md#switches).

- `IGRN=1` (overall correction): the amplitude is multiplied by
  `COR*Γ + (1-COR)*INR`, where `COR` = 1 for `ITR=0`. For `ITR=1`,
  `COR` = $\prod_{i=t,\bar t}\left[m_t\Gamma_t / (\sqrt{s_i}\,\Gamma_t(s_i))\right]^{1/2}$, with
  the LO width $\Gamma_t(s)$ at virtuality $s$ (`CORR`).
- `IGRN=2` (total-angular-momentum prescription): the t- and u-channel amplitudes are split.
  Only the part computed with the modified top propagator (`FVIXXZ`) is multiplied by
  $\Gamma$, and the s-channel amplitude is left unchanged. The two prescriptions differ at
  subleading order. The paper uses `IGRN=1` (its Sec. 2.1).
- The gg channel uses the singlet or octet table according to the colour decomposition.
  $q\bar q$ uses the octet table.

### Where the tables are read from

| Version | Location | Notes |
|---|---|---|
| MadEvent (2010) | `GrnFnc_*/grnep*.tbl`, linked into each `SubProcesses/P*` directory by `link.sh` | read from `./` or `../`, because MadEvent runs jobs in `G*` subdirectories. Grid constants are hard-coded in `SETUPGRN`, in the matrix-element files |
| Stand-alone | the run directory | `scripts/run_standalone.sh` copies the tables from `$TABLES` (default: the `GrnFnc_blvblv` set). Grid constants come from `standalone/common/tblgrnep.inc` |

To use another set, point `TABLES` to its directory:

```bash
TABLES=/work/mytables scripts/run_standalone.sh blvblv test IGRN=1 ...
```

For the MadEvent version, replace the four files in `GrnFnc_*/` of your working copy
before `link.sh`. `scripts/run_paper2010.sh` copies the process directory first, so the
repository stays unchanged if you edit that copy.

### Keep the parameters consistent

The table is only valid for its $m_t$, $\Gamma_t$ and $\alpha_s(\mu_B)$, which are in its
first line. Neither version stops if they differ from the generator's parameters:

- **MadEvent (2010)** compares them with `TMASS` and prints "Parameters inconsistent" to
  unit 9, but carries on.
- **Stand-alone** reads them into `COMMON/GrnPAR/` (used for the $p$ grid) and does not
  compare.

So when you change $m_t$ or $\Gamma_t$ in `param_card.dat` / `param_card_sm.dat`, make new
tables with the same values. The 2010 memo says the same.

### Outside the table

| Region | MadEvent (2010) | Stand-alone |
|---|---|---|
| $E < -50$ GeV | $\Gamma = 1$ (no correction) | extrapolated from the first interval |
| $E > 830$ GeV ($\sqrt{\hat s} \gtrsim 833$ GeV) | $\Gamma = 1$ | extrapolated linearly from the last interval |
| $p > p_{\max}$ | clamped to $p_{\max}$ | clamped to $p_{\max}$ |
| $p < p_{\min}$ (high-energy grid only, where $p_{\min} > 0$) | extrapolated from the first interval | extrapolated from the first interval |

Inside the range the two versions use identical code. The contributions from outside are
small (see #9 and `scripts/check/grn_extrap.f`).

## References

- Y. Sumino and H. Yokoya, JHEP 09 (2010) 034, arXiv:1007.0075, Sec. 2.1 (Green function,
  prescriptions) and Sec. 3 (tabulation).
- Potentials: see the table above. $\alpha_s$: A. Vogt, Comput. Phys. Commun. 170 (2005) 65,
  hep-ph/0408244 (QCD-PEGASUS).
