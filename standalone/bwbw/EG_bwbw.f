C=======================================================================
C     EG_bwbw : stand-alone event generator (BASES/SPRING) for
C       p p -> t t~ -> b W+ b~ W-   (W bosons on shell, not decayed)
C     including the t t~ threshold (bound-state) effects of
C     Y. Sumino and H. Yokoya, JHEP 09 (2010) 034 [arXiv:1007.0075].
C
C     With IGRN>=1 the resonant t t~ amplitude is corrected with the
C     Green-function factor Gamma(E,p) = G/G0 of the colour-singlet
C     or -octet t t~ pair, interpolated from tables
C     (docs/green-function-tables.md). Matrix elements: MadGraph II /
C     HELAS (me_*.f, generated code). Unweighted events: SPRING.
C
C     Inputs (read from the run directory):
C       param_card_sm.dat   SM masses, widths and couplings (SETPARA)
C       grnep1thre.tbl, grnep1high.tbl, grnep8thre.tbl,
C       grnep8high.tbl      Green-function tables (needed if IGRN>=1)
C       topbs.nml           optional namelist /TOPBS/ that overrides
C                           the defaults set below (common/steer.f)
C       LHAPDF set cteq6ll.LHpdf (CTEQ6L1, LO alpha_s)
C     Outputs:
C       stdout, fort.(50+IBSS) settings and SIGMA [pb]; BASES summary
C                           and histograms (BSINFO, BHPLOT)
C       fort.54             SPRING summary and histograms
C       fort.55             BASES grid (written with IBSS=1, read
C                           with IBSS=2)
C       events_bwbw.lhe     unweighted events, Les Houches format
C                           (lheout_bwbw.f): 6 particles per event
C       fort.(31+IGRN)      the two b W invariant masses MT1, MT2
C                           of every event
C     Defaults: LHC 7 TeV, all subprocesses, IGRN=ITR=INR=IKF=0.
C     Run: scripts/run_standalone.sh bwbw <tag> [KEY=VALUE ...],
C     which builds the program, copies the inputs and writes
C     topbs.nml (IBSS=1 unless given). Units: GeV, pb.
C=======================================================================
      PROGRAM MAIN
      IMPLICIT NONE
C     Main COMMON blocks, shared with the integrand in sigh_*.f:
C       /MASS/ m_t, m_b, m_W, Gamma_t [GeV]
C       /INV/  m_tt integration range MIN, MAX [GeV]
C       /CUT/  cut switch ICUT
C       /PRC/  number and list of the active subprocesses
C       /KFAC/ K factors KF(3) and switch IKF
C       /SCL/  renormalisation and factorisation scales [GeV]
C       /THR/  (thr.inc) threshold switches IGRN, ITR, INR
C       run.inc, coupl_sm.inc: MadGraph scale, beam and coupling
C       blocks used by the PDF and matrix-element routines
      INTEGER I
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION RS
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141593D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,GS,ASR,GF, RGAMT
      EXTERNAL ASQCD, RGAMT
      INTEGER IPDF,IPORD,ISUB,IERR
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      DOUBLE PRECISION MUB,ASB
      COMMON/BOHR/MUB,ASB
      INTEGER IBSS,ISPR
      INTEGER IDYN
      COMMON/DYN/IDYN
      INTEGER ICUT
      COMMON/CUT/ICUT
      INTEGER IPART,ICLR,FPROC(4)
      DATA FPROC/1,1,1,1/
      INTEGER NPART,OPROC(4)
      COMMON/PRC/NPART,OPROC
      INTEGER IKF
      COMMON/KFAC/KF,IKF
C     K factors that normalise each subprocess to NLO (IKF=1), from
C     arXiv:1007.0075: KF(1) g g singlet, KF(2) g g octet,
C     KF(3) q q~. kfac.inc selects the LHC 14 TeV values.
      include 'kfac.inc'
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
C     Integrator: IBSS = 0 VEGAS (intvegas.f), 1 BASES (integrate,
C     then save the grid to fort.55), 2 read the grid of an earlier
C     run from fort.55 and skip the integration.
C     ISPR = 1: generate unweighted events with SPRING afterwards.
      IBSS  = 1                 ! 0:VEGAS, 1:BASES, 2:BSREAD
      ISPR  = 1                 ! SPRING Event Generation 0:OFF, 1:ON
      WRITE (50+IBSS,*) "IBSS, ISPR =", IBSS, ISPR
      WRITE (6,*) "IBSS, ISPR =", IBSS, ISPR
C.....
C     PDFs through LHAPDF (CTEQ6L1). IPDF, IPORD, ISUB, IERR are only
C     used by the old SETPDF interface (commented out).
C     Collider: IPRC = 1 p p, 2 p p~; RS = sqrt(s) [GeV].
      IPDF  = 4                 ! 1:GRV98,2:Alekhin02,3:MRST2004NLO,4:CTEQ6
      IPORD = 1                 ! 1:LO, 2:NLO, 3:NNLO
      ISUB  = 0                 !
      IERR  = 0                 !
*     CALL SETPDF (IPDF,IPORD,ISUB,IERR)
*     CALL INITPDFSETBYNAME ('cteq6l.LHpdf') ! NLO as
      CALL INITPDFSETBYNAME ('cteq6ll.LHpdf') ! LO as
      CALL INITPDF (0)
      IPRC  = 1                 ! 1:p-p, 2:p-pbar
*     RS = 1.96D3               ! Tevatron
      RS = 7D3                  ! LHChalf
*     RS = 14D3                 ! LHC
      WRITE (50+IBSS,*) "IPRC, RS =", IPRC, RS
      WRITE (6,*) "IPRC, RS =", IPRC, RS
C...  'run.inc'
C     beam energies [GeV] and beam types (1: proton) in run.inc
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      lpp(1) = 1
      lpp(2) = 1
C.....
C     Physics switches (defaults; topbs.nml may override them):
C       IPART 0 sum, 1 q q~ only, 2 g g only
C       ICLR  0 sum, 1 colour singlet (g g only), 2 colour octet
C       IGRN  0 no bound-state correction, 1 overall correction
C             prescription (OCP), 2 total-angular-momentum one (TAP)
C       ITR   0 naive diagrammatic separation, 1 remove the phase-
C             space suppression in the resonant amplitude (CORR)
C       INR   non-factorizable (non-resonant) diagrams: 0 none,
C             1 fixed-width scheme, 2 overall-factor scheme
C       IDYN  0 fixed scale mu = m_t, 1 or 2 dynamical scale
C             (set in sigh_*.f)
C       IKF   1 multiply each subprocess by its K factor KF
C       ICUT  1 keep |m(bW) - m_t| < 20 GeV for both tops (KINECUT)
      IPART = 0                 ! 0:Sum,1:qq,2:gg
      ICLR  = 0                 ! 0:Sum,1:singlet,2:octet
      IGRN  = 0                 ! 0:No binding, 1:Green function
      ITR   = 0                 ! 0:Naive, 1:Time Retardation
      INR   = 0                 ! 0:NO NON-RES DIGRAMS, 1:INCLUDE
      IDYN  = 0                 ! 0:Fixed scale,1:Dynamical scale
      IKF   = 0                 ! 0:No, 1:Yes
      ICUT  = 0                 ! 0: no cut, 1: knematical cut
C     Optional overrides from ./topbs.nml (also NCALL, ITMX1, ITMX2,
C     NEVENT). The beam energies are set again in case RS changed.
      CALL STEER (IBSS,ISPR,RS,IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT)  ! optional ./topbs.nml
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      WRITE (50+IBSS,*) "IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT =",
     -     IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT
      WRITE (6,*) "IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT =",
     -     IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT
C.....
C     Active subprocesses (index used in OPROC and in ISPRC):
C       1 q q~ (octet), 2 q~ q (octet; q from the second beam),
C       3 g g colour singlet, 4 g g colour octet.
C     The integrand picks one of the NPART active ones at random
C     and multiplies its weight by NPART.
      IF (IPART.EQ.1) THEN
         FPROC(3) = 0
         FPROC(4) = 0
      ELSEIF (IPART.EQ.2) THEN
         FPROC(1) = 0
         FPROC(2) = 0
      ENDIF
      IF (ICLR.EQ.1) THEN
         FPROC(1) = 0
         FPROC(2) = 0
         FPROC(4) = 0
      ELSEIF (ICLR.EQ.2) THEN
         FPROC(3) = 0
      ENDIF
      NPART = 0
      DO I = 1, 4
         IF (FPROC(I).EQ.1) THEN
            NPART = NPART + 1
            OPROC(NPART) = I
         ENDIF
      ENDDO
      IF (NPART.EQ.0) THEN
         WRITE (6,*) "No process"
         STOP
      ENDIF
C.....
C     Masses and widths [GeV] from param_card_sm.dat. The tables
C     grnep*.tbl were made for m_t = 173 GeV, Gamma_t = 1.4911 GeV;
C     the card should agree with them when IGRN>=1.
      call setpara ('param_card_sm.dat',.true.)
      MT = TMASS
      MB = BMASS
      MW = WMASS
C.....
      WRITE (50+IBSS,*) "MT, MB, MW =", MT, MB, MW
      WRITE (6,*) "MT, MB, MW =", MT, MB, MW
      GAMT = TWIDTH
C.....
C     Scales [GeV]: mu_R = mu_F = m_t. MUB = 20 GeV is the scale mu_B
C     of alpha_s in the binding (Coulomb) potential. SCALE: run.inc.
      MUR = MT
      MUF = MUR
      MUB = 20D0
      SCALE = MUR
C.....
C     alpha_s(mu_R) from the PDF set (5 flavours). For a fixed scale
C     it is stored once in run.inc/coupl_sm.inc; for IDYN>=1 the
C     integrand recomputes it at each point.
      INF = 5
      NF = DBLE (INF)
      ASR = ASQCD (MUR,INF)
      IF (IDYN.EQ.0) THEN
         FIXED_REN_SCALE=.TRUE.
         FIXED_FAC_SCALE=.TRUE.
         SCALE = MUR
         Q2FACT(1) = MUF**2
         Q2FACT(2) = Q2FACT(1)
         ALFAS = ASR
         WRITE (50+IBSS,*) "MUR, ASR =", MUR, ASR
         WRITE (6,*) "MUR, ASR =", MUR, ASR
      ELSE
         FIXED_REN_SCALE=.FALSE.
         FIXED_FAC_SCALE=.FALSE.
         WRITE (50+IBSS,*) "Dynamical scale"
         WRITE (6,*) "Dynamical scale"
      ENDIF
C     ASB = alpha_s(mu_B) = 0.1534, the value used to make the tables.
C     COMMON/BOHR is not read by any other routine.
      ASB = 0.1534D0            ! ASQCD (MUB,INF)
C     strong couplings of the HELAS vertices (GG, G) from g_s
      GS = DSQRT (4D0*PI*ASR)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C.....
C     Top width from the card (RGAMT would give the LO width;
C     commented out).
      GF = gfermi
*     GAMT = RGAMT (MT,MW,GF,0D0,MB)
      WRITE (50+IBSS,*) "GAM_T =", GAMT
      WRITE (6,*) "GAM_T =", GAMT, TWIDTH
*     TWIDTH = GAMT
C.....
C     Integration variables (INT12 in sigh_bwbwE.f): X(1), X(2)
C     the two b W invariant masses (Breit-Wigner mapping),
C     X(3)-X(8) decay angles, X(9) rapidity of t t~, X(10) m_tt,
C     X(11) subprocess, X(12) helicity combination (of 144).
      NDIM = 12
C.....
C     m_tt range [GeV]: from the b W+ b~ W- threshold (on-shell
C     W) up to sqrt(s).
      MIN = 2D0 * (MB + MW) + 1D-4 ! God-Hand
      MAX = RS
*     MIN = 330D0
*     MAX = 380D0
      WRITE (50+IBSS,*) "MIN, MAX =", MIN, MAX
      WRITE (6,*) "MIN, MAX =", MIN, MAX
C.....
C     BASES initialisation (BSINIT) and set-up (USERIN). For IBSS=2
C     the grid of an earlier run is read from fort.55.
      IF ( IBSS.EQ.1 ) THEN
         CALL BSINIT
         CALL USERIN
      ELSEIF (IBSS.EQ.2) THEN
         CALL BSINIT
         CALL USERIN
         OPEN (55,FILE="fort.55",FORM="unformatted")
      ENDIF
C.....
C     integrate, then generate events
      CALL SIGBWBW (RS,IBSS,ISPR,NDIM)
C.....
      STOP
      END
C
C
C-----------------------------------------------------------------------
C     SIGBWBW: integrate the cross section, then generate events.
C       RS   : sqrt(s) [GeV]   IBSS, ISPR : see the main program
C       IVAL : number of integration variables (NDIM)
C     The integrand is INT12 (sigh_bwbwE.f), in pb. The result
C     goes to COMMON/BSS/ SIG, ERR [pb], which LHEOUT writes into
C     the <init> block. COMMON/V10/ passes RS, IVAL, IBSS to the
C     integrand; COMMON/RESULT/ receives the BASES/VEGAS result.
C-----------------------------------------------------------------------
      SUBROUTINE SIGBWBW (RS,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER NEVSTR
      EXTERNAL NEVSTR
      INTEGER I,IVAL,JVAL
      DOUBLE PRECISION RS,SIG,ERR
      COMMON/BSS/SIG,ERR
      DOUBLE PRECISION INT12
      EXTERNAL INT12
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS
      INTEGER IBSS,ISPR,JBSS
      COMMON/V10/VRS,JVAL,JBSS
      DOUBLE PRECISION MT1,MT2
      COMMON/BWM/MT1,MT2
      include 'parameter.inc'
      include 'thr.inc'
      INTEGER IH
      COMMON/HST/IH
C     MKHIST fills histograms only when called from INT12
      IH = 12
      SIG = 0D0
      ERR = 0D0
      VRS = RS
      JVAL = IVAL
      JBSS = IBSS
      NVP = 25000
C     VEGAS (intvegas.f) or BASES; S1 = sigma, S2 = its error [pb]
      IF (IBSS.EQ.0) THEN
         CALL VEGAS (INT12,1D-4,IVAL,NVP,10,0,0)
      ELSEIF (IBSS.EQ.1) THEN
         CALL BASES (INT12,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (50+IBSS)
         CALL BHPLOT (50+IBSS)
         SIG = S1
         ERR = S2
         WRITE (50+IBSS,*) "SIGMA =", SIG, " (pb) +-", ERR
         WRITE (6,*) "SIGMA =", SIG, " (pb) +-", ERR
*         CALL XHSAVE2 (10,27)
*         CALL XHSAVE2 (11,28)
*         CALL XHSAVE2 (12,25)
*         CALL XHSAVE2 (13,21)
*         CALL XHSAVE2 (11,32)
*         CALL XHSAVE2 (12,33)
*         CALL XHSAVE2 (13,34)
C     save the BASES grid to fort.55 for later runs with IBSS=2
         CALL BSWRIT (55)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BSREAD (55)
         CLOSE (55)
      ENDIF
C     Event generation: each SPRING call returns one unweighted
C     event (at most 50 trials), LHEOUT(0) writes it to events_bwbw.lhe.
C     LHEOUT(-1) opens the file and writes <init>, LHEOUT(1) closes
C     it. Number of events: 1000000, or NEVENT of topbs.nml.
      IF (ISPR.EQ.1) THEN
         CALL LHEOUT (-1)
         DO I = 1, NEVSTR(1000000)
            CALL SPRING (INT12,50)
C     b W masses MT1, MT2 [GeV] of this event
            WRITE (31+IGRN,*) MT1,MT2
            CALL LHEOUT (0)
         ENDDO
         CALL LHEOUT (1)
         CALL SPINFO (54)
         CALL SHPLOT (54)
      ENDIF
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     USERIN: BASES set-up, called once after BSINIT.
C     All NDIM variables run over [0,1] and are grid-adapted (IG=1).
C     NCALL points per iteration, ITMX1/ITMX2 iterations of the grid
C     and accumulation steps (STEERB applies NCALL, ITMX1, ITMX2 of
C     topbs.nml), ACC1/ACC2 target accuracies [%]. Then the
C     histograms are booked (XHINIT id, xmin, xmax, nbins, title).
C-----------------------------------------------------------------------
      subroutine userin
      implicit none
      include 'parameter.inc'
      integer maxdim
      parameter (maxdim=50)
      double precision xl(maxdim),xu(maxdim)
      integer ig(maxdim)
      integer ncall,itmx1,itmx2
      double precision acc1,acc2
      integer i,idim
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      INTEGER NPART,OPROC(4)
      COMMON/PRC/NPART,OPROC
      nwild = ndim
      do idim = 1, ndim
         xl(idim) = 0d0
         xu(idim) = 1d0
         ig(idim) = 1
      enddo
      call bssetd (ndim,nwild,xl,xu,ig)
      ncall = 1 000 000
      itmx1 = 15
      itmx2 = 15
      call steerb (ncall,itmx1,itmx2)
      acc1  = .1D0
      acc2  = .1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
C     Histograms: 1-10 X(1)-X(10), 11 subprocess, 12 X(12);
C     21-22 MT1, MT2; 23 cos(theta) of t in the parton c.m.
C     frame; 24-26 p_T, |p| and rapidity of t in the lab;
C     27-29 m_tt in three ranges; 31 MT1 vs MT2;
C     32-34 b b~ invariant mass, Delta R, Delta phi.
      do i = 1,10
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      call xhinit (11,0D0,1D0,NPART,'dS/dX')
      call xhinit (12,0D0,1D0,50,'dS/dX')
      CALL XHINIT (21,168D0,178D0,50,'MT1')
      CALL XHINIT (22,168D0,178D0,50,'MT2')
      CALL XHINIT (23,-1D0,1D0,50,'COST')
      CALL XHINIT (24,0D0,500D0,50,'PT_T')
      CALL XHINIT (25,0D0,1000D0,50,'P_T')
      CALL XHINIT (26,-5D0,5D0,50,'ETA_T')
      CALL XHINIT (27,330D0,380D0,50,'MTT')
      CALL XHINIT (28,300D0,800D0,50,'MTT')
      CALL XHINIT (29,MIN,MAX,50,'MTT')
      CALL DHINIT (31,166D0,180D0,50,166D0,180D0,50,'MT1 vs MT2')
      call xhinit (32,0D0,500D0,50,'dS/dMbb')
      call xhinit (33,0D0,5D0,50,'ds/dRbb')
      call xhinit (34,0D0,3.1416D0,50,'dS/dDeltaPhi_bb')
      return
      end
c
C-----------------------------------------------------------------------
C     MKHIST: fill the BASES histograms for one integration point.
C       IS   : id of the calling integrand. Only the call from the
C              top-level one (IS = IH = 12) fills, so each point is
C              counted once
C       IVAL : number of variables X, IBSS : integrator switch
C       SIG  : integrand value = weight [pb], X : variables in [0,1]
C     COMMON/MOM/ P: momenta in the parton c.m. frame,
C       1, 2 incoming, 3 b, 4 W+, 5 b~, 6 W-.
C     PCM2HCM boosts them to the lab frame (PH) with X1, X2 (/XBJ/).
C-----------------------------------------------------------------------
      SUBROUTINE MKHIST (IS,IVAL,IBSS,SIG,X)
      IMPLICIT NONE
      INTEGER I,IS,IVAL,IBSS
      DOUBLE PRECISION SIG,X(10)
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MT1,MT2
      COMMON/BWM/MT1,MT2
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      DOUBLE PRECISION P(0:3,6),PH(0:3,6)
      COMMON/MOM/P
      DOUBLE PRECISION PTOP(0:3),PTB(0:3),PTH(0:3)
      DOUBLE PRECISION COST,PTPT,PABS,ETAT
      DOUBLE PRECISION MBB,RBB,PHIBB
      INTEGER IH
      COMMON/HST/IH
      include 'kin_functions.inc'
      IF (IH.EQ.IS .AND. IBSS.GE.1) THEN
         CALL PCM2HCM (X1,X2,P,PH)
         DO I = 0,3
            PTOP(I) = P(I,3) + P(I,4)
            PTB (I) = P(I,5) + P(I,6)
            PTH (I) = PH(I,3) + PH(I,4)
         ENDDO
C     cos(theta) of t in the parton c.m. frame; p_T, |p| and
C     rapidity of t in the lab frame
         COST = PTOP(3) / DSQRT(PTOP(1)**2+PTOP(2)**2+PTOP(3)**2)
         PTPT = DSQRT(PTH(1)**2+PTH(2)**2)
         PABS = DSQRT(PTH(1)**2+PTH(2)**2+PTH(3)**2)
         ETAT = .5D0 * DLOG((PTH(0)+PTH(3))/(PTH(0)-PTH(3)))
C     b b~ invariant mass [GeV], Delta R and Delta phi
         MBB = DSQRT(SumDot(P(0,3),P(0,5),+1D0))
         RBB = DSQRT(    R2(P(0,3),P(0,5)))
         PHIBB = DELTA_PHI(P(0,3),P(0,5))
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
         CALL XHFILL (21,MT1,SIG)
         CALL XHFILL (22,MT2,SIG)
         CALL XHFILL (23,COST,SIG)
         CALL XHFILL (24,PTPT,SIG)
         CALL XHFILL (25,PABS,SIG)
         CALL XHFILL (26,ETAT,SIG)
         CALL XHFILL (27,MTT,SIG)
         CALL XHFILL (28,MTT,SIG)
         CALL XHFILL (29,MTT,SIG)
         CALL DHFILL (31,MT1,MT2,SIG)
         CALL XHFILL (32,MBB,SIG)
         CALL XHFILL (33,RBB,SIG)
         CALL XHFILL (34,PHIBB,SIG)
      ENDIF
      RETURN
      END
