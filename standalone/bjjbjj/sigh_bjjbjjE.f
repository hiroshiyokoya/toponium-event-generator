C=======================================================================
C     sigh_bjjbjjE.f : integrand of EG_bjjbjj.f for
C     p p -> b j j bbar j j (t, tbar, W+, W- off shell; both W decay
C     to two quarks). The matrix elements are those of blvblv with
C     the colour factors of the hadronic W decays (me_*bjjbjjE.f).
C     with the bound-state correction of Sumino & Yokoya
C     (arXiv:1007.0075). BASES integrates INT18; SPRING generates
C     events from it.
C
C     Call chain. The integrand is d(sigma) [pb] per unit volume of
C     the BASES hypercube [0,1]^18:
C       INT18 (18 dim) : Monte-Carlo sum over subprocess and helicity
C       INT16 (16 dim) : m_tt = sqrt(s-hat) [GeV] from X(16)
C       INT15 (15 dim) : tau = m_tt^2/s, luminosity, pb conversion
C       SIGHAT         : virtualities, angles, |M|^2, K factor, flux
C                        and phase space
C       LUMINO, PTFL   : rapidity and parton luminosity
C       KINECUT        : optional cut (ICUT=1)
C
C     Integration variables X(1..18), all in [0,1]:
C       X(1), X(2)  : top and antitop virtualities (Breit-Wigner)
C       X(3), X(4)  : W+ and W- virtualities (Breit-Wigner)
C       X(5), X(14) : cos(theta), phi of the top (parton c.m. frame)
C       X(6), X(13) : b in the top rest frame
C       X(7), X(12) : bbar in the antitop rest frame
C       X(8), X(11) : W+ decay (particle 4) in the W+ rest frame
C       X(9), X(10) : W- decay (particle 7) in the W- rest frame
C       X(15)       : rapidity of the t tbar system (LUMINO)
C       X(16)       : m_tt
C       X(17)       : subprocess,  X(18) : helicity combination
C     The decay angles are taken with respect to the axes of the
C     parton c.m. frame (MOM2CX + BOOSTX in kinema6.f).
C=======================================================================
C
C-----------------------------------------------------------------------
C     INT18(X): X(17) chooses one of the NPART active subprocesses
C     OPROC(1..NPART) (COMMON/PRC/, set in EG_bjjbjj.f from IPART and
C     ICLR):
C       code 1 : q qbar (quark from beam 1), colour octet
C            2 : qbar q (antiquark from beam 1), colour octet
C            3 : g g, colour singlet
C            4 : g g, colour octet
C     X(18) chooses one helicity combination IHEL (COMMON/HEL/, read
C     by the matrix elements). The factor NPART*NCOMB makes the result
C     a Monte-Carlo estimate of the sum over subprocesses and
C     helicities. MKHIST (EG_bjjbjj.f) fills the histograms.
C     Unlike blvblv and blvbjj, COMMON/PROC/ here holds IPROC, the
C     position in OPROC, not the code OPROC(IPROC). lheout.f reads
C     it as the subprocess code; the two agree only when all four
C     subprocesses are active (IPART=0, ICLR=0).
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT18 (X)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION X(18),XX(16)
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      INTEGER IHEL,IPROC,NCOMB
      COMMON/HEL/IHEL
      COMMON/PROC/IPROC
      DOUBLE PRECISION RS
      INTEGER IVAL,IBSS
      COMMON/V16/RS,IVAL,IBSS
      DOUBLE PRECISION INT16
      EXTERNAL INT16
      INTEGER NPART,OPROC
      COMMON/PRC/NPART,OPROC(4)
      INT18 = 0D0
      IPROC = INT(NPART*X(17)) + 1
      IF (OPROC(IPROC).EQ.1) THEN
         IPART = 1
         ICLR = 2
      ELSEIF (OPROC(IPROC).EQ.2) THEN
         IPART = -1
         ICLR = 2
      ELSEIF (OPROC(IPROC).EQ.3) THEN
         IPART = 2
         ICLR = 1
      ELSEIF (OPROC(IPROC).EQ.4) THEN
         IPART = 2
         ICLR = 2
      ELSE
         RETURN
      ENDIF
C     16 helicity combinations: those listed in the NHEL DATA
C     statements of the matrix elements
      NCOMB = 16
      IHEL = INT(NCOMB*X(18)) + 1
      DO I = 1, 16
         XX(I) = X(I)
      ENDDO
      INT18  = INT16(XX) * DBLE(NPART) * DBLE(NCOMB)
      CALL MKHIST (18,IVAL,IBSS,INT18,X)
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT16(X): t tbar invariant mass MTT = exp(L) [GeV], L uniform in
C     [ln MIN, ln MAX] (COMMON/INV/; EG_bjjbjj.f: 2 m_b + 2 MWCUT1 and
C     sqrt(s)). Jacobian d(MTT)/dX(16) = (ln MAX - ln MIN) MTT. The
C     power-law mapping is kept as comments.
C     IDYN=1: mu_R = mu_F = m_t/2 + MTT/4; alpha_s, g_s (coupl_sm.inc)
C     and SCALE, Q2FACT, ALFAS (run.inc) are reset.
C     MTT goes to INT15 through COMMON/V15/ (and COMMON/TT/).
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT16 (X)
      IMPLICIT NONE
      INTEGER I,A
      DOUBLE PRECISION X(16),XX(15),Z,ZJ,JAC
      DOUBLE PRECISION INT15
      EXTERNAL INT15
      DOUBLE PRECISION RS,MTT
      COMMON/TT/MTT
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      INTEGER IVAL,IBSS
      COMMON/V16/RS,IVAL,IBSS
      DOUBLE PRECISION VRS,VMTT
      INTEGER JVAL,JBSS
      COMMON/V15/VRS,VMTT,JVAL,JBSS
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      DOUBLE PRECISION LMX,LMN,L
      INTEGER IDYN
      COMMON/DYN/IDYN
      INTEGER INF
      DOUBLE PRECISION PI,ASQCD,ASR,GS
      PARAMETER (PI=3.1415923D0,INF=5)
      EXTERNAL ASQCD
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      include 'coupl_sm.inc'
      include 'run.inc'
      INT16 = 0D0
C.....
*     A = 3
*     Z = X(16)**A              !
*     ZJ = A * X(16)**(A-1)
*     MTT = (MAX-MIN)*Z + MIN
*     JAC = (MAX-MIN) * ZJ
C.....
      LMX = DLOG(MAX)
      LMN = DLOG(MIN)
      L = (LMX-LMN)*X(16) + LMN
      MTT = DEXP(L)
      JAC = (LMX-LMN) * MTT
C.....
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
      DO I = 1, 15
         XX(I) = X(I)
      ENDDO
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      INT16 = INT15(XX) * JAC
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT15(X): s = RS^2 (RS = collider energy [GeV]),
C     tau = T = MTT^2/s, s-hat = tau s; JAC = 1/s. X(15) -> rapidity in
C     LUMINO. INT15 = d(sigma)/d(MTT) per unit X(1..15) [pb/GeV].
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT15 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(15),Z
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION JAC,LUM,SIGH
      DOUBLE PRECISION RS,MTT
      INTEGER IVAL,IBSS
      COMMON/V15/RS,MTT,IVAL,IBSS
      DOUBLE PRECISION S,SH,RSH,T
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      INT15 = 0D0
      S = RS**2
      T = (MTT/RS)**2
      SH = T * S
      RSH = DSQRT(SH)
      JAC = 1D0/S
      Z = X(15)
C.....
      CALL SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
      CALL LUMINO (Z,T,IPART,LUM)
      INT15 = SIGH * LUM * JAC
C.....
C     2*MTT/s = d(tau)/d(MTT); 389429.57D3 converts GeV^-2 to pb
C     (about 0.3894 mb GeV^2)
      INT15 = INT15 * 2D0*MTT * 389429.57D3 ! pb/GeV
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     SIGHAT: partonic part. Input X(1..14), RSH = sqrt(s-hat) = MTT
C     [GeV], IPART (1 q qbar, -1 qbar q, 2 g g), ICLR (1 singlet,
C     2 octet). Output SIGH = |M|^2 x flux x phase space [GeV^-2].
C
C     Breit-Wigner mapping of a virtuality s (mass m, width w):
C       s = m^2 + m w tan(z), z uniform between the atan of the
C       limits; JAC_i = (ds/dX)/(2 pi) = (z_max-z_min)
C       [(s-m^2)^2 + m^2 w^2]/(2 pi m w), i.e. it includes the
C       1/(2 pi) of ds/(2 pi). The flat mappings are kept as comments.
C     MWCUT1 = m_W + 0.1 GeV, MWCUT2 = 0.1 GeV (COMMON/CUT/, set in
C     EG_bjjbjj.f): the top virtualities stay above m_b + MWCUT1,
C     so that the width in CORR (ITR=1) does not vanish; MWCUT2 is
C     the lower limit of the W virtualities.
C     The momenta P(0:3,1..8) (parton c.m. frame, KINEMA6) are stored
C     in COMMON/MOM/: g1, g2, b, W+ decay products (4, 5), bbar,
C     W- decay products (7, 8); MR1..MR4 [GeV] in COMMON/BWM/.
C-----------------------------------------------------------------------
      SUBROUTINE SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
      IMPLICIT NONE
      INTEGER IPART,ICLR,IBSS
      DOUBLE PRECISION X(15),RSH,SIGH,MTT
      DOUBLE PRECISION COS1,PHI1,COS2,PHI2,COS3,PHI3,COS4,PHI4,COS5,PHI5
      DOUBLE PRECISION LMD1,LMD2,LMD3,LMD4,LMD5,LAM
      DOUBLE PRECISION JAC,JAC1,JAC2,JAC3,JAC4
      DOUBLE PRECISION SH,ME2
      DOUBLE PRECISION MIN1,MAX1,MIN2,MAX2,MIN3,MAX3,MIN4,MAX4
      DOUBLE PRECISION MR1,MR2,MR3,MR4
      COMMON/BWM/MR1,MR2,MR3,MR4
      DOUBLE PRECISION Z1,Z2,Z3,Z4
      DOUBLE PRECISION Z1MN,Z1MX,Z2MX,Z2MN,Z3MX,Z3MN,Z4MX,Z4MN
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      INTEGER ICUT
      DOUBLE PRECISION MWCUT1,MWCUT2
      COMMON/CUT/MWCUT1,MWCUT2,ICUT
      DOUBLE PRECISION P(0:3,8)
      COMMON/MOM/P
      COMMON/TT/MTT
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
      INTEGER IFLAG
      include 'coupl_sm.inc'
      include 'run.inc'
      ME2 = 0D0
      SIGH = 0D0
      SH = RSH**2
C...  Resonance mass 1
C     Top: MR1 in [m_b+MWCUT1, sqrt(s-hat)-m_b-MWCUT1] from X(1)
      MIN1 = MB + MWCUT1
      MAX1 = RSH - MB - MWCUT1
      Z1MN = DATAN2 (MIN1**2-MT**2, MT*GAMT)
      Z1MX = DATAN2 (MAX1**2-MT**2, MT*GAMT)
      Z1 = (Z1MX-Z1MN)*X(1) + Z1MN !
      MR1 = DSQRT( MT**2 + MT*GAMT*DTAN(Z1) )
      JAC1 = (Z1MX-Z1MN) * ((MR1**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (2D0*PI * MT*GAMT)
*      MR1 = (MAX1-MIN1)*X(1) + MIN1
*      JAC1 = MAX1-MIN1
C...  Resonance mass 2
C     Antitop: MR2 in [m_b+MWCUT1, sqrt(s-hat)-MR1] from X(2)
      MIN2 = MB + MWCUT1
      MAX2 = RSH - MR1
      Z2MN = DATAN2 (MIN2**2-MT**2, MT*GAMT)
      Z2MX = DATAN2 (MAX2**2-MT**2, MT*GAMT)
      Z2 = (Z2MX-Z2MN)*X(2) + Z2MN !
      MR2 = DSQRT( MT**2 + MT*GAMT*DTAN(Z2) )
      JAC2 = (Z2MX-Z2MN) * ((MR2**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (2D0*PI * MT*GAMT)
*      MR2 = (MAX2-MIN2)*X(2) + MIN2
*      JAC2 = MAX2-MIN2
C...  Resonance mass 3
C     W+: MR3 in [MWCUT2, MR1-m_b] from X(3)
      MIN3 = MWCUT2
      MAX3 = MR1 - MB
      Z3MN = DATAN2 (MIN3**2-MW**2, MW*GAMW)
      Z3MX = DATAN2 (MAX3**2-MW**2, MW*GAMW)
      Z3 = (Z3MX-Z3MN)*X(3) + Z3MN !
      MR3 = DSQRT( MW**2 + MW*GAMW*DTAN(Z3) )
      JAC3 = (Z3MX-Z3MN) * ((MR3**2-MW**2)**2 + MW**2*GAMW**2)
     -     / (2D0*PI * MW*GAMW)
*      MR3 = (MAX3-MIN3)*X(3) + MIN3
*      JAC3 = MAX3-MIN3
C...  Resonance mass 4
C     W-: MR4 in [MWCUT2, MR2-m_b] from X(4)
      MIN4 = MWCUT2
      MAX4 = MR2 - MB
      Z4MN = DATAN2 (MIN4**2-MW**2, MW*GAMW)
      Z4MX = DATAN2 (MAX4**2-MW**2, MW*GAMW)
      Z4 = (Z4MX-Z4MN)*X(4) + Z4MN !
      MR4 = DSQRT( MW**2 + MW*GAMW*DTAN(Z4) )
      JAC4 = (Z4MX-Z4MN) * ((MR4**2-MW**2)**2 + MW**2*GAMW**2)
     -     / (2D0*PI * MW*GAMW)
*      MR4 = (MAX4-MIN4)*X(4) + MIN4
*      JAC4 = MAX4-MIN4
C.....
C     Angles: cos(theta) = 2X-1, phi = 2 pi X; with dOmega/(4 pi) in
C     the two-body phase space the Jacobian is 1. 1: top, 2: b in t,
C     3: bbar in tbar, 4: W+ decay, 5: W- decay (rest frames).
      COS1 = 2D0*X(5) - 1D0     !
      COS2 = 2D0*X(6) - 1D0     !
      COS3 = 2D0*X(7) - 1D0     !
      COS4 = 2D0*X(8) - 1D0     !
      COS5 = 2D0*X(9) - 1D0     !
      PHI5 = 2D0*PI*X(10)       !
      PHI4 = 2D0*PI*X(11)       !
      PHI3 = 2D0*PI*X(12)       !
      PHI2 = 2D0*PI*X(13)       !
      PHI1 = 2D0*PI*X(14)       !
      CALL KINEMA6 (RSH, MR1,MR2,MR3,MR4,
     -     COS1,PHI1,COS2,PHI2,COS3,PHI3,COS4,PHI4,COS5,PHI5, P)
C.....
C     IDYN=2: mu_R = mu_F = sqrt(m_t^2 + pT^2), pT of the top
C     (particles 3+4+5); couplings reset as in INT16
      IF (IDYN.EQ.2) THEN
         PT2 = (P(1,3)+P(1,4)+P(1,5))**2 + (P(2,3)+P(2,4)+P(2,5))**2
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
C...  Kinematical Cut
C     Optional cut (ICUT=1), see KINECUT
      IF (ICUT.EQ.1) THEN
         CALL KINECUT (P,IFLAG)
         IF (IFLAG.EQ.1) RETURN
      ENDIF
C.....
C     |M|^2 for helicity combination IHEL, summed over colours and
C     divided by the initial spin and colour average (IDEN in the
C     matrix element). The q qbar matrix element is written for the
C     quark along +z. For IPART = -1 (antiquark from beam 1) the z
C     components of the two incoming momenta are reversed; the
C     reversed momenta stay in COMMON/MOM/. The q qbar matrix element
C     (u ubar) is used for all quark flavours of the luminosity.
      IF ( IPART.EQ.1 ) THEN
         CALL QQBJJBJJ (P,ME2)
      ELSEIF ( IPART.EQ.-1 ) THEN
         P(3,1) = -P(3,1)
         P(3,2) = -P(3,2)
         CALL QQBJJBJJ (P,ME2)
      ELSEIF ( IPART.EQ.2 ) THEN
         CALL GGBJJBJJ (P,ME2,ICLR)
      ENDIF
C.....
C     K factors (IKF=1, COMMON/KFAC/, DATA in EG_bjjbjj.f): KF(1) g g
C     singlet, KF(2) g g octet, KF(3) q qbar; values for LHC 14 TeV
C     (docs/REVIEW.md, P4). IPART = -1 (qbar q) is not covered by
C     this IF block, so that orientation keeps K = 1.
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
C       dPhi = prod_i ds_i/(2 pi) x dPhi2(t tbar) dPhi2(t -> b W+*)
C              dPhi2(tbar -> bbar W-*) dPhi2(W+*) dPhi2(W-*),
C       dPhi2 = lambda^(1/2)/(8 pi) dOmega/(4 pi),
C     LMD = lambda^(1/2) (LAM in kinema6.f); LMD4 = LMD5 = 1 for the
C     massless W decay products.
      JAC  = JAC1 * JAC2 * JAC3 * JAC4
      LMD1 = LAM (1D0,MR1**2/SH,   MR2**2/SH    )
      LMD2 = LAM (1D0,MB**2/MR1**2,MR3**2/MR1**2)
      LMD3 = LAM (1D0,MB**2/MR2**2,MR4**2/MR2**2)
      LMD4 = 1D0
      LMD5 = 1D0
      MTT = RSH
      SIGH = ME2 / (2D0*SH) * LMD1/(8D0*PI) * LMD2/(8D0*PI)
     -     * LMD3/(8D0*PI) * LMD4/(8D0*PI) * LMD5/(8D0*PI) * JAC
      RETURN
      END
C
C-----------------------------------------------------------------------
C     LUMINO: parton luminosity at tau = T = x1 x2 for IPART. X (the
C     variable X(15)) is mapped linearly to the t tbar rapidity
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
C     KINECUT (ICUT=1): all cuts are commented out, IFLAG stays 0.
C     The particle indices are those of sigh_bwblv.f; with the
C     ordering of KINEMA6, P(0,5) is a W+ decay product and P(0,6)
C     the bbar. R2 (kin_functions.f) is Delta R squared.
C-----------------------------------------------------------------------
      SUBROUTINE KINECUT (P,IFLAG)
      IMPLICIT NONE
      INTEGER IFLAG
      DOUBLE PRECISION P(0:3,8), PH(0:3,8)
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MR1,MR2,MR3,MR4
      COMMON/BWM/MR1,MR2,MR3,MR4
      DOUBLE PRECISION RBL1,RBL2,ETB1,ETB2,ETL
      DOUBLE PRECISION R2,ETA
      EXTERNAL R2,ETA
      IFLAG = 0
*     CALL PCM2HCM (X1,X2,P,PH)
      RBL1 = R2(P(0,3),P(0,6))
      RBL2 = R2(P(0,5),P(0,6))
      ETL  = ETA(P(0,6))
      ETB1 = ETA(P(0,3))
      ETB2 = ETA(P(0,5))
*     IF (RBL1.LT.0.2D0) IFLAG = 1
*     IF (RBL2.LT.0.2D0 .OR. RBL2.GT.0.4D0) IFLAG = 1
*     IF (RBL2.GT.0.2D0) IFLAG = 1
      RETURN
      END
