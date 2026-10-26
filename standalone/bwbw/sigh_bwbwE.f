C=======================================================================
C     sigh_bwbwE.f : integrand of EG_bwbw.f for p p -> b W+ bbar W-
C     (off-shell t and tbar, on-shell W bosons) with the bound-state
C     correction of Sumino & Yokoya (arXiv:1007.0075). BASES
C     integrates INT12; SPRING generates events from it.
C
C     Call chain. The integrand is d(sigma) [pb] per unit volume of
C     the BASES hypercube [0,1]^12:
C       INT12 (12 dim) : Monte-Carlo sum over subprocess and helicity
C       INT10 (10 dim) : m_tt = sqrt(s-hat) [GeV] from X(10)
C       INT9  ( 9 dim) : tau = m_tt^2/s, luminosity
C       SIGHAT         : virtualities, angles, |M|^2, K factor, flux
C                        and phase space
C       LUMINO, PTFL   : rapidity and parton luminosity
C       KINECUT        : optional cut (ICUT=1)
C
C     Integration variables X(1..12), all in [0,1]:
C       X(1), X(2)  : top and antitop virtualities (Breit-Wigner)
C       X(3), X(8)  : cos(theta), phi of the top (parton c.m. frame)
C       X(4), X(7)  : b in the top rest frame
C       X(5), X(6)  : bbar in the antitop rest frame
C       X(9)        : rapidity of the t tbar system (LUMINO)
C       X(10)       : m_tt
C       X(11)       : subprocess,  X(12) : helicity combination
C     The decay angles are taken with respect to the axes of the
C     parton c.m. frame (MOM2CX + BOOSTX in kinema4.f).
C=======================================================================
C
C-----------------------------------------------------------------------
C     INT12(X): X(11) chooses one of the NPART active subprocesses
C     OPROC(1..NPART) (COMMON/PRC/, set in EG_bwbw.f from IPART and
C     ICLR):
C       ISPRC = 1 : q qbar (quark from beam 1), colour octet
C               2 : qbar q (antiquark from beam 1), colour octet
C               3 : g g, colour singlet
C               4 : g g, colour octet
C     ISPRC is passed in COMMON/PROC/ to the LHE writer
C     (lheout_bwbw.f). X(12) chooses one helicity combination IHEL
C     (COMMON/HEL/, read by the matrix elements). The factor
C     NPART*NCOMB makes the result a Monte-Carlo estimate of the sum
C     over subprocesses and helicities.
C     MKHIST (EG_bwbw.f) fills the histograms with this weight.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT12 (X)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION X(12),XX(10)
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      INTEGER IHEL,IPROC,NCOMB
      COMMON/HEL/IHEL
      INTEGER ISPRC
      COMMON/PROC/ISPRC
      DOUBLE PRECISION RS
      INTEGER IVAL,IBSS
      COMMON/V10/RS,IVAL,IBSS
      DOUBLE PRECISION INT10
      EXTERNAL INT10
      INTEGER NPART,OPROC
      COMMON/PRC/NPART,OPROC(4)
      INT12 = 0D0
      IPROC = INT(NPART*X(11)) + 1
      ISPRC = OPROC(IPROC)
      IF (ISPRC.EQ.1) THEN
         IPART = 1
         ICLR = 2
      ELSEIF (ISPRC.EQ.2) THEN
         IPART = -1
         ICLR = 2
      ELSEIF (ISPRC.EQ.3) THEN
         IPART = 2
         ICLR = 1
      ELSEIF (ISPRC.EQ.4) THEN
         IPART = 2
         ICLR = 2
      ELSE
         RETURN
      ENDIF
C     144 = 2 x 2 (initial partons) x 2 (b) x 3 (W+) x 2 (bbar)
C           x 3 (W-) helicity combinations
      NCOMB = 144
      IHEL = INT(NCOMB*X(12)) + 1
      DO I = 1, 10
         XX(I) = X(I)
      ENDDO
      INT12 = INT10(XX) * DBLE(NPART) * DBLE(NCOMB)
      CALL MKHIST (12,IVAL,IBSS,INT12,X)
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT10(X): t tbar invariant mass
C       MTT = MIN + (MAX-MIN) z,  z = X(10)**A,  A = 2   [GeV]
C     with MIN, MAX from COMMON/INV/ (EG_bwbw.f: 2(m_b+m_W) + 1D-4
C     and sqrt(s)). Jacobian d(MTT)/dX(10) = (MAX-MIN) A X(10)**(A-1);
C     A = 2 puts more points at small MTT.
C     IDYN=1: mu_R = mu_F = m_t/2 + MTT/4; alpha_s, g_s (coupl_sm.inc)
C     and SCALE, Q2FACT, ALFAS (run.inc) are reset.
C     MTT goes to INT9 through COMMON/V9/ (and COMMON/TT/).
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT10 (X)
      IMPLICIT NONE
      INTEGER I,A
      DOUBLE PRECISION X(10),XX(9),Z,ZJ
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      DOUBLE PRECISION INT9
      EXTERNAL INT9
      DOUBLE PRECISION RS
      INTEGER IVAL,IBSS
      COMMON/V10/RS,IVAL,IBSS
      DOUBLE PRECISION VRS,VMTT
      INTEGER JVAL,JBSS
      COMMON/V9/VRS,VMTT,JVAL,JBSS
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      INTEGER IDYN
      COMMON/DYN/IDYN
      INTEGER INF
      DOUBLE PRECISION PI,ASQCD,ASR,GS
      PARAMETER (PI=3.141592654D0,INF=5)
      EXTERNAL ASQCD
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      include 'coupl_sm.inc'
      include 'run.inc'
      INT10 = 0D0
      A = 2
      Z = X(10)**A
      ZJ = A * X(10)**(A-1)
      MTT = (MAX-MIN)*Z + MIN
      IF (IDYN.EQ.1) THEN
         MUR = MT/2. + MTT/4.
         MUF = MUR
         ASR = ASQCD (MUR,INF)
         GS = DSQRT (4D0*PI*ASR)
         GG(1) = DCMPLX (GS,0D0)
         GG(2) = GG(1)
         G = -GS
         SCALE = MUR
         Q2FACT(1) = MUF**2
         Q2FACT(2) = Q2FACT(1)
         ALFAS = ASR
      ENDIF
      DO I = 1,9
         XX(I) = X(I)
      ENDDO
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      INT10 = INT9(XX) * ZJ * (MAX-MIN)
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT9(X): s = RS^2 (RS = collider energy [GeV]),
C     tau = T = MTT^2/s, s-hat = tau s. JAC = 1/s; together with the
C     factor 2*MTT in SIGHAT it gives d(tau)/d(MTT). X(9) -> rapidity
C     in LUMINO. INT9 = d(sigma)/d(MTT) per unit X(1..9) [pb/GeV].
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT9 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(9),Z
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION JAC,LUM,SIGH
      DOUBLE PRECISION RS,MTT
      INTEGER IVAL,IBSS
      COMMON/V9/RS,MTT,IVAL,IBSS
      DOUBLE PRECISION S,SH,RSH,T
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      INT9 = 0D0
      S = RS**2
      T = (MTT/RS)**2
      SH = T * S
      RSH = DSQRT(SH)
      JAC = 1D0/S
      Z = X(9)
C.....
      CALL SIGHAT (X,RSH,MTT,SIGH,IPART,ICLR,IBSS)
      CALL LUMINO (Z,T,IPART,LUM)
      INT9 = SIGH * LUM * JAC
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     SIGHAT: partonic part. Input X(1..8), RSH = sqrt(s-hat) = MTT
C     [GeV], IPART (1 q qbar, -1 qbar q, 2 g g), ICLR (1 singlet,
C     2 octet). Output SIGH = |M|^2 x flux x phase space x 2*MTT,
C     converted to pb/GeV (multiplied by 1/s and the luminosity in
C     INT9).
C
C     Breit-Wigner mapping of the top virtualities s_i = MT1^2, MT2^2:
C       s = m_t^2 + m_t Gamma_t tan(z), z uniform between the atan of
C       the limits, ds/dz = [(s-m_t^2)^2 + m_t^2 Gamma_t^2]/(m_t
C       Gamma_t). JAC below is the product for t and tbar.
C     The momenta P(0:3,1..6) = g1, g2, b, W+, bbar, W- (parton c.m.
C     frame, KINEMA4) are stored in COMMON/MOM/; MT1, MT2 [GeV] in
C     COMMON/BWM/.
C-----------------------------------------------------------------------
      SUBROUTINE SIGHAT (X,RSH,MTT,SIGH,IPART,ICLR,IBSS)
      IMPLICIT NONE
      INTEGER IPART,ICLR,IBSS
      DOUBLE PRECISION X(9),RSH,MTT,SIGH
      DOUBLE PRECISION COS1,PHI1,COS2,PHI2,COS3,PHI3
      DOUBLE PRECISION LMD1,LMD2,LMD3,LAM,JAC
      DOUBLE PRECISION SH,ME2
      DOUBLE PRECISION MIN1,MAX1,MIN2,MAX2,MT1,MT2
      COMMON/BWM/MT1,MT2
      DOUBLE PRECISION Z1,Z2,Z1MN,Z1MX,Z2MX,Z2MN
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION P(0:3,6)
      COMMON/MOM/P
      INTEGER ICUT,IFLG
      COMMON/CUT/ICUT
      INTEGER IKF
      DOUBLE PRECISION KF(3)
      COMMON/KFAC/KF,IKF
      INTEGER IDYN
      COMMON/DYN/IDYN
      DOUBLE PRECISION PT2
      INTEGER INF
      DOUBLE PRECISION PI,ASQCD,ASR,GS
      PARAMETER (PI=3.141592654D0,INF=5)
      EXTERNAL ASQCD
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      include 'coupl_sm.inc'
      include 'run.inc'
      ME2 = 0D0
      SIGH = 0D0
      SH = RSH**2
C     Top: MT1 in [m_b+m_W, sqrt(s-hat)-m_b-m_W] from X(1)
      MIN1 = MB + MW
      MAX1 = RSH - MIN1
*     IF (MAX1.LE.MIN1) MAX1 = MIN1
      Z1MN = DATAN( (MIN1**2-MT**2)/MT/GAMT )
      Z1MX = DATAN( (MAX1**2-MT**2)/MT/GAMT )
      Z1 = (Z1MX-Z1MN)*X(1) + Z1MN
      MT1 = DSQRT( MT**2 + MT*GAMT*DTAN(Z1) )
      IF (MT1.LT.MB+MW) RETURN
C     Antitop: MT2 in [m_b+m_W, sqrt(s-hat)-MT1] from X(2)
      MIN2 = MB + MW
      MAX2 = RSH - MT1
*     IF (MAX2.LE.MIN2) MAX2 = MIN2
      Z2MN = DATAN( (MIN2**2-MT**2)/MT/GAMT )
      Z2MX = DATAN( (MAX2**2-MT**2)/MT/GAMT )
      Z2 = (Z2MX-Z2MN)*X(2) + Z2MN
      MT2 = DSQRT( MT**2 + MT*GAMT*DTAN(Z2) )
      IF (MT2.LT.MB+MW) RETURN
C     Angles: cos(theta) = 2X-1, phi = 2 pi X. With the normalisation
C     dOmega/(4 pi) of the two-body phase space below, the Jacobian is
C     1. COS1, PHI1: top direction; COS2, PHI2: b in the top rest
C     frame; COS3, PHI3: bbar in the antitop rest frame.
      COS1 = 2D0*X(3) - 1D0     !
      COS2 = 2D0*X(4) - 1D0     !
      COS3 = 2D0*X(5) - 1D0     !
      PHI3 = 2D0*PI*X(6)        !
      PHI2 = 2D0*PI*X(7)        !
      PHI1 = 2D0*PI*X(8)        !
      CALL KINEMA4 (RSH,MT1,MT2,COS1,PHI1,COS2,PHI2,COS3,PHI3,P)
C.....
C     IDYN=2: mu_R = mu_F = sqrt(m_t^2 + pT^2), pT of the top (b + W+);
C     couplings reset as in INT10
      IF (IDYN.EQ.2) THEN
         PT2 = (P(1,3)+P(1,4))**2 + (P(2,3)+P(2,4))**2
         MUR = DSQRT(MT**2 + PT2)
         MUF = MUR
         ASR = ASQCD (MUR,INF)
         GS = DSQRT (4D0*PI*ASR)
         GG(1) = DCMPLX (GS,0D0)
         GG(2) = GG(1)
         G = -GS
         SCALE = MUR
         Q2FACT(1) = MUF**2
         Q2FACT(2) = Q2FACT(1)
         ALFAS = ASR
      ENDIF
C.....
C     Optional cut (ICUT=1, COMMON/CUT/), see KINECUT
      IF (ICUT.EQ.1) THEN
         CALL KINECUT (IFLG)
         IF (IFLG.EQ.1) RETURN
      ENDIF
C.....
C     |M|^2 for helicity combination IHEL, summed over colours and
C     divided by the initial spin and colour average (IDEN in the
C     matrix element). QQBWBW is written for the quark along +z. For
C     IPART = -1 (antiquark from beam 1) the z components of the two
C     incoming momenta are reversed; the reversed momenta stay in
C     COMMON/MOM/. The q qbar matrix element (u ubar) is used for all
C     quark flavours of the luminosity.
      IF ( IPART.EQ.1 ) THEN
         CALL QQBWBW (P,ME2)
      ELSEIF ( IPART.EQ.-1 ) THEN
         P(3,1) = -P(3,1)
         P(3,2) = -P(3,2)
         CALL QQBWBW (P,ME2)
      ELSEIF ( IPART.EQ.2 ) THEN
         CALL GGBWBW (P,ME2,ICLR)
      ENDIF
C.....
C     K factors (IKF=1, COMMON/KFAC/, kfac.inc): KF(1) g g singlet,
C     KF(2) g g octet, KF(3) q qbar. The values are those for LHC
C     14 TeV (docs/REVIEW.md, P4). IPART = -1 (qbar q) is not covered
C     by this IF block, so that orientation keeps K = 1.
      IF (IKF.EQ.1) THEN
         IF (IPART.EQ.1) THEN
            ME2 = ME2 * KF(3)
         ELSEIF (IPART.EQ.2 .AND. ICLR.EQ.1) THEN
            ME2 = ME2 * KF(1)
         ELSEIF (IPART.EQ.2 .AND. ICLR.EQ.2) THEN
            ME2 = ME2 * KF(2)
         ENDIF
      ENDIF
C.....
C     Cross section: flux 1/(2 s-hat) and phase space
C       dPhi = ds_t/(2 pi) ds_tbar/(2 pi) dPhi2(t tbar)
C              dPhi2(t -> b W+) dPhi2(tbar -> bbar W-),
C       dPhi2 = lambda^(1/2)/(8 pi) dOmega/(4 pi),
C     LMD = lambda^(1/2) (LAM in kinema4.f). The 1/(2 pi)^2 is written
C     explicitly. 2*MTT = d(s-hat)/d(MTT); 389429.57D3 converts
C     GeV^-2 to pb (about 0.3894 mb GeV^2).
      JAC  = (Z1MX-Z1MN) * ((MT1**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (MT*GAMT)
     -     * (Z2MX-Z2MN) * ((MT2**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (MT*GAMT)
      LMD1 = LAM (1D0,MT1**2/SH,   MT2**2/SH   )
      LMD2 = LAM (1D0,MB**2/MT1**2,MW**2/MT1**2)
      LMD3 = LAM (1D0,MB**2/MT2**2,MW**2/MT2**2)
C.....
      SIGH = ME2/(2D0*SH) * LMD1/(8D0*PI) * LMD2/(8D0*PI)
     -     * LMD3/(8D0*PI) / (2D0*PI)**2 * JAC
     -     * 2D0*MTT * 389429.57D3 ! pb/GeV
      RETURN
      END
C
C-----------------------------------------------------------------------
C     LUMINO: parton luminosity at tau = T = x1 x2 for IPART. X (the
C     variable X(9)) is mapped linearly to the t tbar rapidity
C       y in [ln(T)/2, -ln(T)/2],  x1 = sqrt(T) e^y,  x2 = sqrt(T) e^-y,
C     so that dx1 dx2 = d(tau) dy; the Jacobian DY is included in LUM.
C     x1, x2 are stored in COMMON/XBJ/ and in xbk (run.inc). The PDFs
C     are taken at mu_F (COMMON/SCL/). LUM = QQB (IPART=1), QBQ (-1)
C     or GG (2), times DY; negative values are set to 0.
C-----------------------------------------------------------------------
      SUBROUTINE LUMINO (X,T,IPART,LUM)
      IMPLICIT NONE
      INTEGER IPART
      DOUBLE PRECISION X,T,LUM
      DOUBLE PRECISION YMIN,YMAX,DY,Y
      DOUBLE PRECISION X1,X2,QQB,QBQ,GG,QG,GQ
      DOUBLE PRECISION JAC
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      COMMON/XBJ/X1,X2
      include 'run.inc'
      YMIN = DLOG(T)/2D0
      YMAX = -YMIN
      DY = YMAX - YMIN
      Y = DY*X + YMIN
      X1 = DSQRT(T) * DEXP(Y)
      X2 = DSQRT(T) * DEXP(-Y)
C...  'run.inc'
      xbk(1) = x1
      xbk(2) = x2
C.....
      CALL PTFL (X1,X2,QQB,QBQ,GG,QG,GQ,MUF**2,IPRC)
      IF (IPART.EQ.1) THEN
         LUM = QQB
      ELSEIF (IPART.EQ.-1) THEN
         LUM = QBQ
      ELSEIF (IPART.EQ.2) THEN
         LUM = GG
      ENDIF
      JAC = DY
      LUM = LUM * JAC
      IF (LUM.LT.0D0) LUM = 0D0
      RETURN
      END
C
C-----------------------------------------------------------------------
C     PTFL (version of the EG_* generators; common/ptfl.f has another
C     interface): products of parton densities at scale^2 MU2 [GeV^2]
C       GG  : g(X1) g(X2)
C       QQB : sum over q = d,u,s,c of q(X1) qbar(X2)
C       QBQ : sum over q = d,u,s,c of qbar(X1) q(X2)
C     b quarks are not included. IPRC = 1 pp, 2 p pbar (beam 2 is an
C     antiproton: its qbar density is the proton density of q).
C     QG and GQ are not set. CALLPDF returns x f(x), hence the
C     division by X1 X2.
C-----------------------------------------------------------------------
      SUBROUTINE PTFL (X1,X2,QQB,QBQ,GG,QG,GQ,MU2,IPRC)
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2,QQB,QBQ,GG,QG,GQ,MU2
      DOUBLE PRECISION MU
      INTEGER IPRC, IEXP
      DOUBLE PRECISION PDF1(-6:6), PDF2(-6:6)
      EXTERNAL CALLPDF
      INTEGER IPDF,IPORD,ISUB,IERR
      COMMON/PDF/IPDF,IPORD,ISUB,IERR
      MU = DSQRT (MU2)
      CALL CALLPDF (MU,X1,PDF1)
      CALL CALLPDF (MU,X2,PDF2)
      GG = PDF1(0) * PDF2(0) / X1 / X2
      IF(IPRC.EQ.1) THEN        ! P-P
         QQB = ( PDF1( 1)*PDF2(-1) + PDF1( 2)*PDF2(-2)
     -        + PDF1( 3)*PDF2(-3) + PDF1( 4)*PDF2(-4) ) / X1 / X2
         QBQ = ( PDF1(-1)*PDF2( 1) + PDF1(-2)*PDF2( 2)
     -        + PDF1(-3)*PDF2( 3) + PDF1(-4)*PDF2( 4) ) / X1 / X2
      ELSEIF (IPRC.EQ.2) THEN   ! P-PBAR
         QQB = ( PDF1( 1)*PDF2( 1) + PDF1( 2)*PDF2( 2)
     -        + PDF1( 3)*PDF2( 3) + PDF1( 4)*PDF2( 4) ) / X1 / X2
         QBQ = ( PDF1(-1)*PDF2(-1) + PDF1(-2)*PDF2(-2)
     -        + PDF1(-3)*PDF2(-3) + PDF1(-4)*PDF2(-4) ) / X1 / X2
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C     KINECUT (ICUT=1): IFLG = 1 (point rejected) if a top virtuality
C     MT1 or MT2 differs from m_t by more than 20 GeV.
C-----------------------------------------------------------------------
      SUBROUTINE KINECUT (IFLG)
      IMPLICIT NONE
      INTEGER IFLG
      DOUBLE PRECISION P(0:3,6),PH(0:3,6)
      COMMON/MOM/P
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION MMVJ
      COMMON/MVJ/MMVJ
      DOUBLE PRECISION MT1,MT2
      COMMON/BWM/MT1,MT2
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION SUMDOT
      EXTERNAL SUMDOT
      IFLG = 0
*      CALL PCM2HCM (X1,X2,P,PH)
      IF (DABS(MT1-MT).GT.20D0) IFLG = 1
      IF (DABS(MT2-MT).GT.20D0) IFLG = 1
      RETURN
      END
