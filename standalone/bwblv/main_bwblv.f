C=======================================================================
C     main_bwblv : t t~ invariant-mass distribution
C       d(sigma)/d(m_tt) [pb/GeV]  for
C       p p -> t t~ -> b W+ b~ mu- nu_mu~   (LHC 14 TeV;
C       W+ on shell and not decayed, W- -> mu- nu_mu~),
C     computed at fixed m_tt with BASES, including the bound-state
C     effects of arXiv:1007.0075 (IGRN>=1).
C     Integrand INT12 of sigh_bwblv.f (m_tt fixed); matrix elements
C     me_ggbwblv.f, me_qqbwblv.f (helicities summed inside).
C     A calculator without steering: the settings are fixed below
C     (default: g g colour singlet only, IGRN=1); no K factors and
C     no event file. In its present form the scan loop is commented
C     out and only m_tt = 344 GeV is computed (see below).
C
C     Inputs: param_card_sm.dat, grnep{1,8}{thre,high}.tbl (IGRN>=1),
C             LHAPDF set cteq6l.LHpdf (via SETPDF).
C     Outputs: stdout   m_tt, d(sigma)/d(m_tt), error
C              fort.7   m_tt, d(sigma)/d(m_tt)
C              fort.2   settings, BASES summary
C     Build: make PROC=bwblv (standalone/Makefile); run the program
C     in a directory that holds the input files. Units: GeV, pb.
C=======================================================================
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I
C     COMMON blocks shared with sigh_bwblv.f: /MASS/ masses and
C     widths [GeV], /ID/ IPART, ICLR, /SCL/ scales [GeV], /CUT/ cut,
C     /THR/ (thr.inc) IGRN, ITR, INR.
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS,MTT,DSDM,ERR
      DOUBLE PRECISION MIN,MAX,DM
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141593D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,GS,AS0,GF, RGAMT
      EXTERNAL ASQCD, RGAMT
      INTEGER IPDF,IPORD,ISUB,IERR
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      DOUBLE PRECISION MUB,ASB
      COMMON/BOHR/MUB,ASB
      INTEGER IBSS,ISPR,IDYN
      INTEGER ICUT
      DOUBLE PRECISION MWCUT
      COMMON/CUT/MWCUT,ICUT
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
C     Integrator: IBSS = 1 VEGAS, 2 BASES (note: not the numbering
C     of EG_*.f). ISPR = 1 would also generate 1000 SPRING events
C     per integration, used only for the SPRING histograms.
      IBSS  = 2                 ! 1:VEGAS, 2:BASES
      ISPR  = 0                 ! SPRING Event Generation 0:OFF, 1:ON
      WRITE (6,*) "IBSS, ISPR =", IBSS, ISPR
C.....
C     PDFs through SETPDF (pdfport_lhapdf.f): IPDF=4, IPORD=1 loads
C     cteq6l.LHpdf (CTEQ6L, NLO alpha_s). Collider: IPRC = 1 p p,
C     2 p p~; RS = sqrt(s) [GeV].
      IPDF  = 4                 ! 1:GRV98,2:Alekhin02,3:MRST2004NLO,4:CTEQ6
      IPORD = 1                 ! 1:LO, 2:NLO, 3:NNLO
      ISUB  = 0                 !
      IERR  = 0                 !
      CALL SETPDF (IPDF,IPORD,ISUB,IERR)
      IPRC  = 1                 ! 1:p-p, 2:p-pbar
      RS = 14D3                 ! LHC
C...  'run.inc'
C     beam energies [GeV] and beam types (1: proton) in run.inc
      ebeam(1) = rs/2.
      ebeam(2) = rs/2.
      lpp(1) = 1
      lpp(2) = 1
C.....
C     Physics switches (no topbs.nml steering in this program):
C       IPART 0 sum of q q~ and g g, 1 q q~ only, 2 g g only
C       ICLR  0 singlet + octet, 1 colour singlet, 2 colour octet
C             (selects the g g colour state)
C       IGRN  0 no bound-state correction, 1 OCP, 2 TAP (thr.inc)
C       ITR, INR  see thr.inc (used for g g in GGTTRES, ggtt.f)
C       IDYN  0 fixed scale mu = m_t, 1 mu_R = mu_F = m_tt
C       ICUT  1 KINECUT (sigh_bwblv.f): keep Delta R(b~,mu-)**2 <= 0.2
C     No K factors are applied.
      IPART = 2                 ! 0:Sum, 1:QQ, 2:GG
      ICLR  = 1                 ! 0:Sum, 1:Singlet, 2:Octet
      IGRN  = 1                 ! 0:No binding, 1:OCP,2:TAP
      ITR   = 0                 ! 0:
      INR   = 0                 ! 0:
      IDYN  = 0                 ! 0:Fixed scale, 1:Dynamical scale
      ICUT  = 0                 ! 0: no cut, 1: kinematical cut
      WRITE (2,*) "IPART,ICLR,IGRN,IDYN,ICUT =",
     -     IPART,ICLR,IGRN,IDYN,ICUT
      WRITE (6,*) "IPART,ICLR,IGRN,IDYN,ICUT =",
     -     IPART,ICLR,IGRN,IDYN,ICUT
C.....
C     Masses and widths [GeV] from param_card_sm.dat (SETPARA)
      call setpara ('param_card_sm.dat',.true.)
      MT = TMASS
      MB = BMASS
      MW = WMASS
      WRITE (2,*) "MT, MB, MW =", MT, MB, MW
      WRITE (6,*) "MT, MB, MW =", MT, MB, MW
      GAMT = TWIDTH
      GAMW = WWIDTH
C.....
C     Scales [GeV]: mu_R = mu_F = m_t; MUB = 20 GeV is the scale
C     mu_B of the binding (Coulomb) alpha_s.
      MUR = MT
      MUF = MUR
      MUB = 20D0
C.....
C     alpha_s(mu_R) and alpha_s(mu_B) from the PDF set (5 flavours);
C     ASB (COMMON/BOHR) is not read by any other routine. g_s goes
C     into the HELAS couplings GG, G.
      INF = 5
      NF = DBLE (INF)
      AS0 = ASQCD (MUR,INF)
      ASB = ASQCD (MUB,INF)
      GS = DSQRT (4D0*PI*AS0)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C
C     Top width: the LO width RGAMT (gamt.f) for the masses of the
C     card replaces the card value, also in HELAS (TWIDTH).
*      GAMT = 1.4911D0
      GF = gfermi
      GAMT = RGAMT (MT,MW,GF,0D0,MB)
*      GAMT = RGAMT (MT,MW,GF,0D0,0D0)
      WRITE (6,*) "GAM_T, GAM_W =", GAMT, GAMW
      TWIDTH = GAMT
C
C     Integration variables (INT12): X(1) b W+ invariant mass
C     (W+ on shell), X(2) t~ and X(3) W- virtualities,
C     X(4)-X(11) decay angles, X(12) rapidity.
      NDIM = 12
C
C     MWCUT [GeV]: lower limit of the W- virtuality; the t~
C     virtuality starts at m_b + MWCUT (SIGHAT).
      MWCUT = 5D0
      MIN = 2D0*MB + MW + MWCUT
      MAX = RS
C
*      MIN = 2D0*MT - 15D0
*      MAX = 360D0
C     m_tt range [GeV] of the scan (overrides the MIN, MAX above)
      MIN = 330D0
*      MIN = 400D0
      MAX = 3550D0
      DM = .5D0
*      MIN = 2D0*MT
*     MAX = 1000D0
*      DM = 1D0
      MTT = MIN
C
C     BASES initialisation and set-up (USERIN)
      IF ( IBSS.EQ.2 ) THEN
         CALL BSINIT
         CALL USERIN
      ENDIF
C
C     Single point m_tt = 344 GeV: the DO statement of the scan
C     is commented out, so the block below runs once (labels
C     99 and 100 are left over from the loop).
C     With the loop active: at each m_tt integrate over the remaining
C     variables (SIGBWBLV) and print d(sigma)/d(m_tt). The step
C     DM shrinks to 0.5, 0.2, 0.1 GeV from 2m_t-15, 2m_t-10,
C     2m_t-3 GeV and grows to 2, 5, 20, 50 GeV from 2m_t+4, 400,
C     500, 1000 GeV.
      MTT = 344D0
*      DO 99 I = 1, 500
C     IDYN=1: scales and couplings at mu = m_tt
         IF (IDYN.EQ.1) THEN
            MUR = MTT
            MUF = MUR
            AS0 = ASQCD (MUR,INF)
            GS = DSQRT (4D0*PI*AS0)
            GG(1) = DCMPLX (GS,0D0)
            GG(2) = GG(1)
            G = -GS
         ENDIF
         WRITE (2,*) "MTT =", MTT
         CALL SIGBWBLV (RS,MTT,DSDM,ERR,IBSS,ISPR,NDIM)
C     m_tt [GeV], d(sigma)/d(m_tt) [pb/GeV] (and error)
         WRITE (6,*) MTT, DSDM, ERR
         WRITE (7,*) MTT, DSDM
         IF (MTT.GE.2D0*MT-15D0) DM = .5D0
         IF (MTT.GE.2D0*MT-10D0) DM = .2D0
         IF (MTT.GE.2D0*MT-3D0) DM = .1D0
         IF (MTT.GE.2D0*MT+4D0) DM = 2D0
         IF (MTT.GE.400D0) DM = 5D0
         IF (MTT.GE.500D0) DM = 20D0
         IF (MTT.GE.1000D0) DM = 50D0
         IF ( MTT.GE.MAX ) GOTO 100
         MTT = MTT + DM
 99   CONTINUE
 100  CONTINUE
      STOP
      END
C
C
C-----------------------------------------------------------------------
C     SIGBWBLV: d(sigma)/d(m_tt) [pb/GeV] at fixed m_tt.
C       RS, MTT : sqrt(s) and m_tt [GeV] (input)
C       DSDM, ERR : result and its error [pb/GeV] (output)
C       IBSS : 1 VEGAS, 2 BASES; ISPR : 1 also run SPRING
C       IVAL : number of integration variables (NDIM)
C     Integrand INT12 (sigh_bwblv.f) through COMMON/V12/;
C     COMMON/RESULT/ receives the BASES/VEGAS result S1, S2.
C     JJ in COMMON/V12/ is not used by INT12.
C-----------------------------------------------------------------------
      SUBROUTINE SIGBWBLV (RS,MTT,DSDM,ERR,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER I,IVAL,JVAL
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS,MTT, DSDM,ERR
      DOUBLE PRECISION GFN,INT12
      EXTERNAL GFN,INT12
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS,VMTT
      INTEGER IBSS,ISPR,JBSS,JJ
      COMMON/V12/VRS,VMTT,JVAL,JBSS,JJ
      INTEGER IH
      COMMON/HST/IH
      include 'parameter.inc'
C     MKHIST fills histograms only when called from INT12
      IH = 12
      DSDM = 0D0
      ERR = 0D0
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      JJ = 9
      NVP = 10000
C     VEGAS (intvegas.f) or BASES; S1 = result, S2 = error
      IF (IBSS.EQ.1) THEN
         CALL VEGAS (INT12,1D-3,IVAL,NVP,8,0,0)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BASES (INT12,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (2)
         CALL BHPLOT (2)
*         CALL XHSAVE (3,24)
         IF (ISPR.EQ.1) THEN
            DO I = 1,1000
               CALL SPRING (INT12,50)
            ENDDO
            CALL SPINFO (2)
            CALL SHPLOT (2)
         ENDIF
      ENDIF
      DSDM = S1
      ERR = S2
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     USERIN: BASES set-up, called once after BSINIT. All NDIM
C     variables in [0,1], grid-adapted; NCALL points per iteration,
C     ITMX1/ITMX2 iterations, ACC1/ACC2 target accuracies [%].
C     Histograms 1..NDIM: the integration variables X(i).
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
      nwild = ndim
      do idim = 1, ndim
         xl(idim) = 0d0
         xu(idim) = 1d0
         ig(idim) = 1
      enddo
      call bssetd (ndim,nwild,xl,xu,ig)
      ncall = 25000
      itmx1 = 20
      itmx2 = 20
      acc1  = 1D0
      acc2  = 1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
      do i = 1,ndim
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      return
      end
c
C-----------------------------------------------------------------------
C     MKHIST: fill the histograms of the integration variables for
C     one point (only for IBSS=2 and calls from the top-level
C     integrand, IS = IH = 12). SIG is the weight [pb/GeV or pb].
C     COMMON/MOM/ P in the parton c.m. frame: 1, 2 incoming,
C     3 b, 4 W+, 5 b~, 6 mu-, 7 nu_mu~. COST, PTPT, PABS, ETA (top
C     kinematics) are computed but not filled.
C-----------------------------------------------------------------------
      SUBROUTINE MKHIST (IS,IVAL,IBSS,SIG,X)
      IMPLICIT NONE
      INTEGER I,IS,IVAL,IBSS
      DOUBLE PRECISION SIG,X(12)
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MT1,MT2,MW1
      COMMON/BWM/MT1,MT2,MW1
      DOUBLE PRECISION P(0:3,7),PH(0:3,7)
      COMMON/MOM/P
      DOUBLE PRECISION PT(0:3),PTB(0:3),PTH(0:3)
      DOUBLE PRECISION COST,PTPT,PABS,ETA
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      INTEGER IH
      COMMON/HST/IH
      IF (IH.EQ.IS .AND. IBSS.EQ.2) THEN
         CALL PCM2HCM (X1,X2,P,PH)
         DO I = 0,3
            PT (I) = P(I,3) + P(I,4)
            PTB(I) = P(I,5) + P(I,6)
            PTH (I) = PH(I,3) + PH(I,4)
         ENDDO
         COST = PT(3) / DSQRT(PT(1)**2+PT(2)**2+PT(3)**2)
         PTPT = DSQRT(PTH(1)**2+PTH(2)**2)
         PABS = DSQRT(PTH(1)**2+PTH(2)**2+PTH(3)**2)
         ETA = .5D0 * DLOG((PTH(0)+PTH(3))/(PTH(0)-PTH(3)))
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
      ENDIF
      RETURN
      END
