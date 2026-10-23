C     ==============================================================
C     GrnEPhigh.f -- Green-function correction table, high-energy grid
C     ==============================================================
C     Computes the correction factor
C        Gamma(E,p) = G(E,p)/G0(E,p),  G0 = -1/(E + i*Gamma_t - p^2/m_t)
C     for the colour-singlet (ICLR=1) and colour-octet (ICLR=2) ttbar
C     system (arXiv:1007.0075, Eq. (3.1)). G(E,p) is the S-wave
C     momentum-space Green function, E the energy measured from the
C     threshold 2*m_t, p = |p_t| in the ttbar rest frame.
C     Units: GeV and GeV^-1.
C
C     Method, for each colour and each E on the grid:
C      1. SOLVE, REFINE: the radial equation for u(r) = r*G(r)
C            u'' + m_t*(E + i*Gamma_t - V(r))*u = 0
C         (reduced mass m_t/2) is integrated with 4th-order
C         Runge-Kutta in x = log10(r). The solution is normalised to
C         u -> 1 for r -> 0 and decays for r -> infinity.
C      2. FOURIER: sine transform of u to momentum space (Simpson
C         rule in log r plus an analytic piece for r < R_UV), then
C         Gamma = -(E + i*Gamma_t - p^2/m_t) * Gtilde(E,p).
C      3. SAVETABLE: writes the table of the colour state.
C
C     Inputs, set in the main program below (keep them identical in
C     GrnEPthre.f and GrnEPhigh.f):
C       MU     top-quark mass m_t [GeV]
C       GAMT   top-quark width Gamma_t [GeV] (or RGAMT, gamt.f)
C       MUB    scale mu_B of the potential [GeV]
C       ASB    alpha_s(mu_B), from QCD-PEGASUS (alphas_pegasus.f)
C       Vqcd.f the QCD potential VQCD(r): copy Vqcd_1LFO.f (default),
C              Vqcd_2LFO.f or Vqcd_2LRGI.f to Vqcd.f before make
C     Grid constants: tblgrnep.inc. Step numbers: incgrnep.inc.
C
C     Outputs (in the current directory):
C       grnep1high.tbl, grnep8high.tbl   tables (singlet, octet)
C       fort.1, fort.2    r, V(r) of singlet, octet (check)
C       fort.21           convergence log of SOLVE (the "R->oo
C                         converge" lines go to fort.11)
C       fort.22           E, r0, (u, du/dx) at r0 = 10**XIN
C       fort.23           E, R_UV, r0, R_IR
C       fort.33, fort.34  E, Im B for singlet, octet, i.e.
C                         Im{G}(E,r=0) up to m_t/(4 pi) (see SOLVE)
C
C     High-energy vs threshold program: this file uses a uniform E
C     grid, E = EMIN2..EMAX2 = 30..830 GeV (NE2 = 400 intervals),
C     and for each E a p grid of NP2+1 = 102 points around the
C     on-shell momentum p0 = sqrt(m_t*E), geometric with ratio RR on
C     each side (see PP). It uses NLOOPH (SOLVE) and NLOOP2 (REFINE,
C     FOURIER) steps, ten times more than GrnEPthre.f, which covers
C     E = -50..30 GeV and p = 0..300 GeV on uniform grids.
C
C     Build and run (docs/green-function-tables.md), in a copy of
C     this directory, since tables and fort.* go to the cwd:
C       make LD=gfortran FFLAGS="-O3 -ffixed-line-length-132
C            -std=legacy -w"
C       ./GrnEPhigh.exe > high.log
C     ==============================================================
      PROGRAM MAIN
      IMPLICIT NONE
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION MUB,ASB
      COMMON/GrnPOT/MUB,ASB
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      DOUBLE PRECISION XIN,XUV,XIR
      COMMON/GrnX/XIN,XUV,XIR
      DOUBLE PRECISION X,R,V,VQCD
      EXTERNAL VQCD
C...  External program
      INTEGER NF,NORD
      DOUBLE PRECISION ALPHAS
      DOUBLE PRECISION RGAMT
      DOUBLE PRECISION GF,MW,MB,VTB
      PARAMETER (GF=1.16639D-5,MW=80.4D0,MB=5D0,VTB=1D0)
*     PARAMETER (GF=1.17801D-5,MW=80D0,MB=4.7D0,VTB=0.998D0)
      DOUBLE PRECISION MZ,ASZ
      PARAMETER (MZ=91.2D0,ASZ=0.118D0)
*     PARAMETER (MZ=91.2D0,ASZ=0.12D0)
      DOUBLE PRECISION ALFST,ALFSP,RQCD
C...  Top Mass, Width
      MU = 173D0
      GAMT = 1.4911D0
*     MU = 176.17D0
*     GAMT = 1.602D0
*     GAMT = RGAMT (MU,GF,MW,MB,VTB)
*     WRITE (6,*) MU,GAMT
C     alpha_s(mu_B) at mu_B = MUB, run from alpha_s(M_Z) = ASZ with
C     NF = 5 flavours at NLO (QCD-PEGASUS); this gives 0.15329.
C     Note (docs/REVIEW.md P11): the shipped tables carry 0.1534 in
C     their header; they were made with the commented-out
C     ASB = 0.1534D0 below.
C...  CALL alphas_pegasus.f
      NF = 5                    ! num of flavor
      NORD = 1                  ! 1:NLO,2:NNLO
      MUB = 20D0
      ASB = ALPHAS (MUB,MZ,ASZ,NF,NORD)
      WRITE (6,*) MUB,ASB
C...  or give them directly
*     MUB = 20D0
*     ASB = 0.1534D0
C...  Potential and Color setting
*     ICLR = 1                  ! Color of ttbar pair
C     Loop over the colour states; VQCD (Vqcd.f) reads ICLR.
      DO 100 ICLR = 1, 2        ! 1:Singlet,2:Octet
      WRITE (6,*) ICLR,MU,GAMT,ASB
C     Note (docs/REVIEW.md P12): Vqcd_2LRGI.f needs this SETVPR call
C     (POTQCD.f). It is commented out, so with the 2LRGI potential
C     the program crashes (already in the potential check below).
C     SETVPR takes alpha_s(M_Z) = ASZ, not MUB, ASB.
*     CALL SETVPR (ASZ,MU,MB,MZ,GAMT,ALFST,ALFSP,RQCD)
C...  Setting Initial Conditions
C     In x = log10(r [GeV^-1]): XIN is the start point r0 = 0.1 of
C     the integration, XUV and XIR the first trial end points towards
C     r -> 0 (1e-7) and r -> infinity (100) (see SOLVE).
      XIN = -1D0
      XUV = -7D0
      XIR =  2D0
C...  Check QCD Potential
C     Write r, V(r) at 101 points to fort.ICLR (fort.1, fort.2).
      DO X = XUV, XIR, (XIR-XUV)/100.
         R = 10D0**X
         V = VQCD(R)
         WRITE (ICLR,*) R,V
      ENDDO
C.....
      CALL CALCGRNEP
C.....
 100  CONTINUE
      STOP
      END
C
C
C     --------------------------------------------------------------
C     CALCGRNEP: loop over the energy grid E = EE(IE), IE = 0..NE2,
C     for the colour ICLR. For each E, SOLVE (which calls REFINE)
C     solves the radial equation and FOURIER transforms the solution
C     to all p of the grid, filling GRNEP(IE,IP). SAVETABLE then
C     writes the table.
C     --------------------------------------------------------------
      SUBROUTINE CALCGRNEP
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INCLUDE 'incgrnep.inc'
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER IE, EINT,PINT
      DOUBLE PRECISION EE,PP
      EXTERNAL EE,PP,EINT,PINT
      DOUBLE PRECISION DIMB
      COMMON/GrnIMG/DIMB
C.....
C     Values (u, du/dx) at x0 = XIN of the two independent solutions
C     GA, GB used in SOLVE: (1, 1) and (1, -1).
      GA0(1) = DCMPLX( 1D0,0D0)
      GA0(2) = DCMPLX( 1D0,0D0)
      GB0(1) = DCMPLX( 1D0,0D0)
      GB0(2) = DCMPLX(-1D0,0D0)
C...  Energy Loop
      DO 50 IE = 0, NE2
         E = EE (IE)
C     p0 = sqrt(m_t*E), the on-shell momentum: centre of the p grid
C     (PP, PINT) through COMMON/GrnOS/.
         POS = DSQRT(MU*E)
*     WRITE (6,*) IE,E
C...
         CALL SOLVE (IE)
         CALL FOURIER (IE)
         WRITE (6,*) IE, E, " done, ", DIMB
C...
 50   CONTINUE
 51   CONTINUE
C...  WRITE TABLE TO THE FILE
      CALL SAVETABLE
C...
 100  CONTINUE
      RETURN
      END
C
C     --------------------------------------------------------------
C     SOLVE: solves u'' + m_t*(E + i*Gamma_t - V(r))*u = 0 for the
C     current E and colour, with GG(1) = u, GG(2) = du/dx in
C     x = log10(r). Two solutions GA, GB start at x0 = XIN with the
C     values GA0, GB0 and are integrated with RUNGE_KUTTA, first
C     towards r -> 0 (step 1), then towards r -> infinity (step 2).
C     From their values at a radius r it forms
C       G0 = A1*GA + A2*GB   with G0 = 0, dG0/dr = 1 at r
C       G1 = A5*GA + A6*G0   with G1 = 1, Im dG1/dr = 0 at r
C     taken at r = R_UV (proxy for r -> 0), and at r = R_IR
C       u  = A7*GA + A8*G0   with u = 0 at R_IR (and u = 1 at R_UV)
C       B  = -G1/G0          so that u = G1 + B*G0.
C     With G(r) = m_t*u/(4 pi r) (see FOURIER), Im G(E,0) is
C     m_t/(4 pi) * Im B. Output: GP0 = (u, du/dx) at x0 (COMMON
C     /GrnSOL/), R_UV and R_IR (COMMON/GrnREG/); then REFINE.
C     --------------------------------------------------------------
      SUBROUTINE SOLVE (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE=(1D0,0D0))
      INTEGER I,J
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      INTEGER IE
      DOUBLE PRECISION EPS,ERR
      DOUBLE PRECISION H,X,R
      DOUBLE PRECISION R0,X1,R1
      DOUBLE PRECISION XIN,XUV,XIR
      COMMON/GrnX/XIN,XUV,XIR
      DOUBLE COMPLEX GA(N),GAN(N),DGA
      DOUBLE COMPLEX GB(N),GBN(N),DGB
      DOUBLE COMPLEX GP(N),GPN(N)
      DOUBLE PRECISION GPMN
      DOUBLE COMPLEX A1,A2,A5,A6,A7,A8,A1P,A2P,A7P,A8P,BP
      DOUBLE COMPLEX B,B1,B2,GP0(N)
      COMMON/GrnSOL/B1,B2,GP0
      DOUBLE COMPLEX G0N,G1N
      DOUBLE PRECISION DIMB
      DOUBLE PRECISION DEL,D0
      INTEGER I1,I2,I7,I8
      COMMON/GrnIMG/DIMB
C     EPS: accuracy of u = 1 at R_UV (step 3). ERR = EPS/NLOOPH:
C     tolerance on the relative change of a coefficient in one step.
C...  Define Accuracy
      EPS = 1D-4
      ERR = EPS/DBLE(NLOOPH)
C     --- Step 1: r -> 0. Integrate GA, GB from x0 = XIN to X1 = XUV
C     (H < 0). After each step compute A1, A2 (regular solution G0),
C     A5 = A7 = 1/GA and A6. A flag I1, I2, I7 is set once the
C     relative change of A1, A2, A7 in one step is below ERR (it is
C     not reset afterwards). When all three are set, R_UV = r and
C     the coefficients of that step are kept.
C...  Take R->0 Limit
C...  Set Final Point of Evaluation
 90   CONTINUE
      X0 = XIN
      X1 = XUV
 91   CONTINUE
C...  Initialize
      X = X0
      R0 = 10D0**X0
      DO I=1,N
         GA(I) = GA0(I)
         GB(I) = GB0(I)
      ENDDO
C...  Reflesh Flag
      I1=0
      I2=0
      I7=0
      I8=0
C...  Start to solve
      H = (X1-X0)/NLOOPH        ! Step of Evaluations
      DO I=1,NLOOPH,1
         CALL RUNGE_KUTTA (X,H,GA,GAN)
         CALL RUNGE_KUTTA (X,H,GB,GBN)
         X = X+H
         R = 10D0**X
         DGA = 1D0/L10/R * GAN(2)
         DGB = 1D0/L10/R * GBN(2)
         A1 = CONE/(DGA - DGB*GAN(1)/GBN(1))
         A2 = -A1*GAN(1)/GBN(1)
         A5 = CONE/GAN(1)
         A6 = DCMPLX(0D0,-DIMAG(A5*DGA))
         A7 = A5
C     Note (docs/REVIEW.md P6): A1P, A2P, A7P are not initialised
C     before the first step of the first energy.
         IF ( I1.EQ.0 .AND. ABS((A1-A1P)/A1P).LE.ERR ) I1 = 1
         IF ( I2.EQ.0 .AND. ABS((A2-A2P)/A2P).LE.ERR ) I2 = 1
         IF ( I7.EQ.0 .AND. ABS((A7-A7P)/A7P).LE.ERR ) I7 = 1
         IF ( I1*I2*I7 .EQ. 1 ) THEN
            RUV = R
            WRITE (21,*) IE, R, "; R->0 converge"
            GOTO 100
         ENDIF
         DO J=1,N
            GA(J) = GAN(J)
            GB(J) = GBN(J)
         ENDDO
         A1P = A1
         A2P = A2
         A7P = A7
      ENDDO
C     Not converged: move the end point one decade closer to r = 0
C     and restart from x0.
C     Note (docs/REVIEW.md P13): the number of retries is unbounded.
      WRITE (21,*) IE, X1, "; R->0 NOT CONVERGE"
      X1 = X1 - 1D0
      GOTO 91
 100  CONTINUE
C     --- Step 2: r -> infinity. Integrate GA, GB again from x0 to
C     X1 = XIR. With A1, A2, A5, A6, A7 from step 1 form G0, G1,
C     A8 = -A7*GA/G0 (so that u = A7*GA + A8*G0 vanishes at r) and
C     B = -G1/G0.
C
C     Convergence: DEL, the relative change of A8 in one step, has to
C     drop below ERR (D0 = ERR after a reset) and then not increase;
C     I8 doubles at each such step and is reset to 1 otherwise.
C     I8 > 1000, i.e. 10 such steps in a row, gives R_IR = r.
C     I8 starts at 0 here (1 in GrnEPthre.f), so the counting starts
C     only after the first reset.
C...  Take R->Infinity Limit
C...  Set Final Point of Evaluation
      X1 = XIR
 101  CONTINUE
C...  Initialize
      X = X0
      DO I=1,N
         GA(I) = GA0(I)
         GB(I) = GB0(I)
      ENDDO
C...  Start to solve
      H = (X1-X0)/NLOOPH        ! Step of Evaluations
      DO I=1,NLOOPH,1
         CALL RUNGE_KUTTA (X,H,GA,GAN)
         CALL RUNGE_KUTTA (X,H,GB,GBN)
         X = X+H
         R = 10D0**X
         DGA = 1D0/L10/R * GAN(2)
         DGB = 1D0/L10/R * GBN(2)
         G0N = A1*GAN(1) + A2*GBN(1)
         G1N = A5*GAN(1) + A6*G0N
         A8 = - A7*GAN(1)/G0N
         B = -G1N/G0N
C     Note (docs/REVIEW.md P6): A8P and D0 are not initialised
C     before the first step of the first energy.
         DEL = ABS((A8-A8P)/A8P)
         IF ( DEL.LE.D0 ) THEN
            I8 = I8*2
            D0 = DEL
         ELSE
            I8 = 1
            D0 = ERR
         ENDIF
         IF ( I8.GT.1000 ) THEN
            RIR = R
            WRITE (11,*) IE, X, "; R->oo converge"
            GOTO 110
         ENDIF
         DO J=1,N
            GA(J) = GAN(J)
            GB(J) = GBN(J)
         ENDDO
         A8P = A8
         BP = B
      ENDDO
C     Not converged: extend the end point by 0.2 in log10(r) and
C     restart from x0. NLOOPT is increased, but the loops of this
C     program use NLOOPH, so only the end point (and the step size)
C     changes.
C     Note (docs/REVIEW.md P13): the number of retries is unbounded.
      WRITE (21,*) IE, X1, "; R->INFINITY, NOT CONVERGE"
      NLOOPT = NLOOPT + 5000
      X1 = X1 + .2
      GOTO 101
 110  CONTINUE
C     Physical solution u = A7*GA + A8*G0 = B1*GA + B2*GB; GP0 holds
C     its (u, du/dx) at x0. Im B is written to fort.33/34.
      WRITE (21,*) IE, E, ": Runge-Kutta Solved ", NLOOPH
      WRITE (21,*) "R_0 =", R0
      B1 = A7 + A1*A8
      B2 =      A2*A8
      GP0(1) = B1*GA0(1) + B2*GB0(1)
      GP0(2) = B1*GA0(2) + B2*GB0(2)
      DIMB = DIMAG(B)
      WRITE (22,*) E, R0, GP0
      WRITE (32+ICLR,*) E, DIMB ! PRINT Im{G}
C     --- Step 3: integrate u itself (GP0 at x0) to refine the end
C     points. R_UV is moved out to the first r, coming from x0, where
C     u = (1, 0) within EPS; if this never happens R_UV is kept.
C...  Check End Points
C...  Re-calc R_UV where the function sufficiently reach to (1.,0.)
C...  Set Final Point of Evaluation
      R1 = RUV
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  R -> R_IR
      H = (X1-X0)/NLOOP2        ! Step of Evaluations
      DO I=1,NLOOP2,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
         IF ( ABS(DREAL(GPN(1))-1D0).LT.EPS .AND.
     -        ABS(DIMAG(GPN(1))-0D0).LT.EPS ) THEN
            WRITE (21,*) "R_UV", RUV, " ->", R
            RUV = R
            GOTO 10
         ENDIF
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
 10   CONTINUE
C     R_IR: integrate u out to 1.5*R_IR. R_IR becomes the first r
C     with |u| < 1e-4, or else the r where |u| is smallest.
C     Note (docs/REVIEW.md P5): the cut |u| < 1e-4 ignores the weight
C     r*u of the p -> 0 integrand in FOURIER; the dropped tail
C     r > R_IR gives errors of a few % at small p for E > 0.
C...  Re-calc R_IR as the minimum point of ABS(GP)
C...  Set Final Point of Evaluation
      R1 = RIR * 1.5
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
      GPMN = ABS(GP0(1))
C...  R -> R_IR
      H = (X1-X0)/NLOOP2        ! Step of Evaluations
      DO I=1,NLOOP2,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
         IF ( ABS(GPN(1)) .LT. GPMN ) THEN
            R0 = R
            GPMN = DBLE(ABS(GPN(1)))
         ENDIF
         IF ( ABS(GPN(1)) .LT. 1D-4 ) THEN
            R0 = R
            GOTO 95
         ENDIF
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
 95   CONTINUE
      WRITE (21,*) "R_IR:", RIR, " ->", R0
      RIR = R0
      R0 = 10D0**X0
C     fort.23: E, R_UV, r0 = 10**XIN, R_IR.
      WRITE (23,*) E, RUV, R0, RIR
      CALL REFINE (IE)
      RETURN
      END
C
C     --------------------------------------------------------------
C     REFINE: tabulates u(r) for FOURIER. Starting from GP0 at
C     r0 = 10**XIN, it integrates the physical solution with NLOOP2
C     equal steps in log10(r):
C       1: from r0 down to R_UV -> RP1(i), REGP1(i), IMGP1(i)
C       2: from r0 up to R_IR   -> RP2(i), REGP2(i), IMGP2(i)
C     i = 0..NLOOP2; i = 0 is r0 in both (COMMON/GrnTBL/).
C     --------------------------------------------------------------
      SUBROUTINE REFINE (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      INTEGER I,J
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER IE
      DOUBLE PRECISION H,X,R
      DOUBLE PRECISION R0,X1,R1
      DOUBLE COMPLEX GP(N),GPN(N)
      DOUBLE COMPLEX B1,B2,GP0(N)
      COMMON/GrnSOL/B1,B2,GP0
C...  1: R->0 direction
C...  Set Final Point of Evaluation
      R1 = RUV
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  Store starting point data
      R0 = 10D0**X0
      REGP1 (0) = DREAL(GP0(1))
      IMGP1 (0) = DIMAG(GP0(1))
      RP1 (0) = R0
      REGP2 (0) = DREAL(GP0(1))
      IMGP2 (0) = DIMAG(GP0(1))
      RP2 (0) = R0
*      WRITE (24,*) R0,DREAL(GP0(1)),DIMAG(GP0(1)),ABS(GP0(1))
C...  R -> R_UV
      H = (X1-X0)/NLOOP2        ! Step of Evaluations
      DO I=1,NLOOP2,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
C...
         REGP1 (I) = DREAL(GPN(1))
         IMGP1 (I) = DIMAG(GPN(1))
         RP1 (I) = R
*         WRITE (24,*) R,DREAL(GPN(1)),DIMAG(GPN(1)),ABS(GPN(1))
C...
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
C...  2: R->oo direction
C...  Set Final Point of Evaluation
      R1 = RIR
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  R -> R_IR
      H = (X1-X0)/NLOOP2        ! Step of Evaluations
      DO I=1,NLOOP2,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
C...
         REGP2 (I) = DREAL(GPN(1))
         IMGP2 (I) = DIMAG(GPN(1))
         RP2 (I) = R
*         WRITE (25,*) R, DREAL(GPN(1)),DIMAG(GPN(1)),ABS(GPN(1))
C...
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
C...
      RETURN
      END
C
C     --------------------------------------------------------------
C     FOURIER: momentum-space Green function for all p of the grid.
C     With G(r) = m_t*u(r)/(4 pi r), the S-wave transform is
C       Gtilde(p) = (m_t/p) * int_0^inf dr sin(p r) u(r)
C                 = int d(ln r) r*u(r)*S1(p,r),  S1 = m_t*sin(p r)/p.
C     The r range is split into
C       [R_UV, r0], [r0, R_IR]: Simpson rule in ln r on the tables
C                               of REFINE (steps DL1, DL2)
C       [0, R_UV]: u = 1, analytic, m_t*(1-cos(p R_UV))/p^2 (SIGUV)
C     and r > R_IR is dropped. Then
C       Gamma(E,p) = -(E + i*Gamma_t - p^2/m_t) * Gtilde = G/G0
C     is stored in GRNEP(IE,IP).
C     --------------------------------------------------------------
      SUBROUTINE FOURIER (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      INTEGER I,IE,IP
      DOUBLE COMPLEX CIMG
      PARAMETER ( CIMG=(0D0,1D0))
      DOUBLE PRECISION P,PP
      EXTERNAL PP
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION DL1,DL2,SIGRE,SIGIM,SIGUV
      DOUBLE PRECISION SIGRE1,SIGIM1,SIGRE2,SIGIM2
      DOUBLE PRECISION S1
      EXTERNAL S1
      DOUBLE PRECISION SMP(2)
      DATA SMP/2D0,4D0/         ! /Even, Odd/
      DOUBLE COMPLEX GAMMA
      DOUBLE COMPLEX GRNEP(0:NE2,0:NP2)
      COMMON/GrnMOM/GRNEP
      DL1 = DABS(DLOG(RP1(1))-DLOG(RP1(0)))
      DL2 = DABS(DLOG(RP2(1))-DLOG(RP2(0)))
      DO 100 IP = 0,NP2,1
         P = PP (IP)
C...  Daikei Koushiki
*         SIGRE1 = 0D0
*         SIGRE2 = 0D0
*         SIGIM1 = 0D0
*         SIGIM2 = 0D0
*         DO I = 0,NLOOP2-1,1
*            SIGRE1 = SIGRE1 + DL1 * ( RP1(I)*REGP1(I)*S1(P,RP1(I),MU)
*     -           + RP1(I+1)*REGP1(I+1)*S1(P,RP1(I+1),MU) ) / 2D0
*            SIGRE2 = SIGRE2 + DL2 * ( RP2(I)*REGP2(I)*S1(P,RP2(I),MU)
*     -           + RP2(I+1)*REGP2(I+1)*S1(P,RP2(I+1),MU) ) / 2D0
*            SIGIM1 = SIGIM1 + DL1 * ( RP1(I)*IMGP1(I)*S1(P,RP1(I),MU)
*     -           +  RP1(I+1)*IMGP1(I+1)*S1(P,RP1(I+1),MU) ) / 2D0
*            SIGIM2 = SIGIM2 + DL2 * ( RP2(I)*IMGP2(I)*S1(P,RP2(I),MU)
*     -           + RP2(I+1)*IMGP2(I+1)*S1(P,RP2(I+1),MU) ) / 2D0
*         ENDDO
*         WRITE (6,*) SIGRE1,SIGRE2,SIGIM1,SIGIM2
C     Simpson weights in i: (1, 4, 2, 4, ..., 2, 4, 1) * DL/3;
C     SMP(MOD(I,2)+1) = 4 for odd i, 2 for even i (NLOOP2 is even).
C...  Simpthon Integral Formula
         SIGRE1 = DL1/3D0 * ( RP1(0)*REGP1(0)*S1(P,RP1(0),MU)
     -        + RP1(NLOOP2)*REGP1(NLOOP2)*S1(P,RP1(NLOOP2),MU) )
         SIGRE2 = DL2/3D0 * ( RP2(0)*REGP2(0)*S1(P,RP2(0),MU)
     -        + RP2(NLOOP2)*REGP2(NLOOP2)*S1(P,RP2(NLOOP2),MU) )
         SIGIM1 = DL1/3D0 * ( RP1(0)*IMGP1(0)*S1(P,RP1(0),MU)
     -        + RP1(NLOOP2)*IMGP1(NLOOP2)*S1(P,RP1(NLOOP2),MU) )
         SIGIM2 = DL2/3D0 * ( RP2(0)*IMGP2(0)*S1(P,RP2(0),MU)
     -        + RP2(NLOOP2)*IMGP2(NLOOP2)*S1(P,RP2(NLOOP2),MU) )
C...
         DO I = 1,NLOOP2-1, 1
            SIGRE1 = SIGRE1 + DL1/3D0 *
     -           SMP(MOD(I,2)+1)*RP1(I)*REGP1(I)*S1(P,RP1(I),MU)
            SIGRE2 = SIGRE2 + DL2/3D0 *
     -           SMP(MOD(I,2)+1)*RP2(I)*REGP2(I)*S1(P,RP2(I),MU)
            SIGIM1 = SIGIM1 + DL1/3D0 *
     -           SMP(MOD(I,2)+1)*RP1(I)*IMGP1(I)*S1(P,RP1(I),MU)
            SIGIM2 = SIGIM2 + DL2/3D0 *
     -           SMP(MOD(I,2)+1)*RP2(I)*IMGP2(I)*S1(P,RP2(I),MU)
         ENDDO
*         WRITE (6,*) SIGRE1,SIGRE2,SIGIM1,SIGIM2
C...
C     Piece r < R_UV with u = 1 (real part only).
C     Note (docs/REVIEW.md P9): at p = 0 the factor P makes SIGUV = 0;
C     the intended value is MU*RUV**2/2 (negligible, R_UV ~ 1e-6).
         IF (P.EQ.0D0) THEN
            SIGUV = P*MU*RUV**2/2D0
         ELSE
            SIGUV = MU/P**2*(1D0-DCOS(P*RUV))
         ENDIF
         SIGRE = SIGRE1 + SIGRE2 + SIGUV
         SIGIM = SIGIM1 + SIGIM2
*         WRITE (6,*) P,SIGRE1, SIGRE2, SIGUV
C...
C     Gamma = G/G0 with G0 = -1/(E + i*Gamma_t - p^2/m_t).
         GAMMA = - (E + CIMG*GAMT - P**2/MU) * DCMPLX(SIGRE,SIGIM)
         GRNEP(IE,IP) = GAMMA
*         WRITE (6,*) PP(IP), DREAL(GAMMA),DIMAG(GAMMA), IE
 100  CONTINUE
      RETURN
      END
C
C     --------------------------------------------------------------
C     SAVETABLE: writes GRNEP to grnep1high.tbl (ICLR = 1) or
C     grnep8high.tbl (ICLR = 2), list-directed:
C       line 1: MU, GAMT, ASB    (m_t, Gamma_t, alpha_s(mu_B))
C       line 2: NE2, NP2         (numbers of intervals)
C       line 3: EMIN2, EMAX2
C       line 4: RR               (ratio of the p grid)
C     then for IE = 0..NE2: NP2+1 lines "Re Gamma, Im Gamma"
C     (IP = 0..NP2), closed by a line "--- IE ---".
C     The readers (readgrnep.f) must use the same tblgrnep.inc.
C     --------------------------------------------------------------
      SUBROUTINE SAVETABLE
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION MUB,ASB
      COMMON/GrnPOT/MUB,ASB
      INTEGER IE,IP,IOF
      PARAMETER (IOF=28)
      DOUBLE COMPLEX GRNEP(0:NE2,0:NP2)
      COMMON/GrnMOM/GRNEP
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      INCLUDE 'filedata.inc'
      IF (ICLR.EQ.1) THEN
         WRITE (6,*) "   Open file : ", H1
         OPEN (IOF,FILE=H1)
      ELSEIF (ICLR.EQ.2) THEN
         WRITE (6,*) "   Open file : ", H8
         OPEN (IOF,FILE=H8)
      ELSE
         STOP
      ENDIF
      WRITE (IOF,*) MU, GAMT, ASB
      WRITE (IOF,*) NE2, NP2
      WRITE (IOF,*) EMIN2, EMAX2
      WRITE (IOF,*) RR
      DO IE = 0,NE2,1
         DO IP = 0,NP2,1
            WRITE (IOF,*) DREAL(GRNEP(IE,IP)), DIMAG(GRNEP(IE,IP))
         ENDDO
         WRITE (IOF,*) "--- ", IE, " ---"
      ENDDO
      CLOSE (IOF)
      RETURN
      END
C
C     --------------------------------------------------------------
C     RUNGE_KUTTA: one classical 4th-order Runge-Kutta step of size
C     H in x = log10(r) for GG = (u, du/dx): GG -> GGN.
C     K1FUNC..K4FUNC give the four stages (each returns H*f).
C     --------------------------------------------------------------
      SUBROUTINE RUNGE_KUTTA (X,H,GG,GGN)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H
      DOUBLE COMPLEX GG(N),GGN(N)
      DOUBLE COMPLEX KK1(N),KK2(N),KK3(N),KK4(N)
      CALL K1FUNC (X,H,GG,    KK1)
      CALL K2FUNC (X,H,GG,KK1,KK2)
      CALL K3FUNC (X,H,GG,KK2,KK3)
      CALL K4FUNC (X,H,GG,KK3,KK4)
      DO I=1,N
         GGN(I) = GG(I) + 1D0/6D0 *
     -        (KK1(I) + 2D0*KK2(I) + 2D0*KK3(I) + KK4(I))
      ENDDO
      RETURN
      END
C
C     --------------------------------------------------------------
C     RKFUNC: FF = H * d(GG)/dx. With r = 10**x,
C       d2u/dx2 = ln10 * du/dx + (ln10*r)**2 * d2u/dr2,
C     and d2u/dr2 = -FFUNC = -m_t*(E + i*Gamma_t - V(r))*u.
C     --------------------------------------------------------------
      SUBROUTINE RKFUNC (X,H,GG,FF)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE PRECISION X,H,R
      DOUBLE COMPLEX GG(N), FF(N)
      DOUBLE COMPLEX FFUNC
      EXTERNAL FFUNC
      R = 10D0**X
      FF(1) = H * ( GG(2) )
      FF(2) = H * ( L10*GG(2) - L10**2*R**2 * FFUNC (X,GG(1)) )
      RETURN
      END
C
C     FFUNC = m_t*(E + i*Gamma_t - V(r))*u at r = 10**x. V = VQCD
C     (Vqcd.f), which depends on ICLR, MUB and ASB.
      DOUBLE COMPLEX FUNCTION FFUNC (X,G1)
      IMPLICIT NONE
      DOUBLE COMPLEX CIMG
      PARAMETER ( CIMG=(0D0,1D0))
      DOUBLE PRECISION X,R
      DOUBLE COMPLEX G1
      DOUBLE PRECISION VQCD
      EXTERNAL VQCD
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      R = 10D0**X
      FFUNC = MU * (E + CIMG*GAMT - VQCD (R)) * G1
      RETURN
      END
C
C     K1FUNC..K4FUNC: Runge-Kutta stages at x, x+H/2, x+H/2, x+H.
      SUBROUTINE K1FUNC (X,H,GG,KK1)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X1
      DOUBLE COMPLEX GG(N),GG1(N),KK1(N)
      X1 = X
      DO I=1,N
         GG1(I) = GG(I)
      ENDDO
      CALL RKFUNC (X1,H,GG1,KK1)
      RETURN
      END
C
      SUBROUTINE K2FUNC (X,H,GG,KK1,KK2)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X2
      DOUBLE COMPLEX GG(N),GG2(N),KK1(N),KK2(N)
      X2 = X + H/2D0
      DO I=1,N
         GG2(I) = GG(I) + KK1(I)/2D0
      ENDDO
      CALL RKFUNC (X2,H,GG2,KK2)
      RETURN
      END
C
      SUBROUTINE K3FUNC (X,H,GG,KK2,KK3)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X3
      DOUBLE COMPLEX GG(N),GG3(N),KK2(N),KK3(N)
      X3 = X + H/2D0
      DO I=1,N
         GG3(I) = GG(I) + KK2(I)/2D0
      ENDDO
      CALL RKFUNC (X3,H,GG3,KK3)
      RETURN
      END
C
      SUBROUTINE K4FUNC (X,H,GG,KK3,KK4)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X4
      DOUBLE COMPLEX GG(N),GG4(N),KK3(N),KK4(N)
      X4 = X + H
      DO I=1,N
         GG4(I) = GG(I) + KK3(I)
      ENDDO
      CALL RKFUNC (X4,H,GG4,KK4)
      RETURN
      END
C
C     --------------------------------------------------------------
C     S1(p,r) = m_t*sin(p r)/p, kernel of the sine transform;
C     m_t*r at p = 0 and a Taylor series for p*r < 1e-3.
C     --------------------------------------------------------------
      DOUBLE PRECISION FUNCTION S1 (P,R,MU)
      IMPLICIT NONE
      DOUBLE PRECISION P,R,MU
      IF (P .EQ. 0D0) THEN
         S1 = MU * R
      ELSEIF (P*R .LT. 1D-3) THEN
         S1 = MU * ( R - R**3*P**2/6D0 + R**5*P**4/120D0 )
      ELSE
         S1 = MU/P * DSIN(P*R)
      ENDIF
      RETURN
      END
C
C
C     --------------------------------------------------------------
C     Grid functions of the high-energy table:
C       EE(IE)   E of grid point IE = 0..NE2 (uniform, EMIN2..EMAX2)
C       EINT(E)  index of the E interval containing E
C       PP(IP), PINT(P): p grid around p0 = POS (see PP)
C     EINT and PINT are not called in this program.
C     --------------------------------------------------------------
      DOUBLE PRECISION FUNCTION EE (IE)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IE
      EE = EMIN2 + DBLE(IE)/DBLE(NE2) * (EMAX2 - EMIN2)
      RETURN
      END
C
      INTEGER FUNCTION EINT (E)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION E
      EINT = DINT(NE2*(E-EMIN2)/(EMAX2-EMIN2))
      RETURN
      END
C
C     --------------------------------------------------------------
C     PP(IP): p grid around p0 = POS = sqrt(m_t*E), IP = 0..NP2,
C     M = (NP2-1)/2, R1 = 1/RR:
C       IP <= M: p0*(1-RR**(IP+1))/(1-RR**(M+1)); PP(M) = p0,
C                PP(0) > 0, spacing shrinks by RR towards p0
C       IP >  M: p0 + A*RR**M*(1-R1**(IP-M+1))/(1-R1),
C                spacing grows by 1/RR away from p0
C     --------------------------------------------------------------
      DOUBLE PRECISION FUNCTION PP (IP)
      IMPLICIT NONE
      INTEGER IP,M
      DOUBLE PRECISION A,R1
      INCLUDE 'tblgrnep.inc'
      M = (NP2-1)/2
      A = (1D0-RR)/(1D0-RR**M) * POS
      R1 = 1D0/RR
      IF (IP.LE.M) THEN
         PP = POS * (1D0-RR**(IP+1)) / (1D0-RR**(M+1))
      ELSE
         PP = POS + A*RR**M*(1D0-R1**(IP-M+1))/(1D0-R1)
      ENDIF
      RETURN
      END
C
C     --------------------------------------------------------------
C     PINT(P): index IP of the p interval [PP(IP), PP(IP+1)]
C     containing P (inverse of PP).
C     Note (docs/REVIEW.md P7): for P in (p0, p0 + A*RR**M) it
C     returns M-1 instead of M, inconsistent with PP. The readers
C     (readgrnep.f) use the same formula.
C     --------------------------------------------------------------
      INTEGER FUNCTION PINT (P)
      IMPLICIT NONE
      INTEGER M
      DOUBLE PRECISION P,A,R1
      INCLUDE 'tblgrnep.inc'
      M = (NP2-1)/2
      A = (1D0-RR)/(1D0-RR**M) * POS
      R1 = 1D0/RR
      IF (P.LE.POS) THEN
         PINT = -1 + DINT(DLOG(1D0-P/POS*(1D0-RR**(M+1)))/DLOG(RR))
      ELSE
         PINT = DINT(DLOG(1D0-(P-POS)/(A*RR**M)*(1D0-R1))/DLOG(R1))
     -        + M - 1
      ENDIF
      RETURN
      END
