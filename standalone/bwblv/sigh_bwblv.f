C=======================================================================
C     sigh_bwblv.f : integrand for p p -> b W+ bbar mu- vm~ (on-shell
C     W+, off-shell t, tbar and W-) with the bound-state correction.
C     Older version, used by main_bwblv.f (d(sigma)/d(m_tt) at fixed
C     m_tt, integrand INT12) and tot_bwblv.f (sigma integrated over
C     m_tt, integrand INT13). Helicities and colours are summed inside
C     the matrix elements. For q qbar the luminosity QQ of
C     common/ptfl.f (quark from either beam) multiplies the matrix
C     element with the quark along +z.
C
C     Integration variables, all in [0,1]:
C       X(1), X(2)  : top and antitop virtualities (Breit-Wigner)
C       X(3)        : W- virtuality (Breit-Wigner)
C       X(4), X(11) : cos(theta), phi of the top (parton c.m. frame)
C       X(5), X(10) : b in the top rest frame
C       X(6), X(9)  : bbar in the antitop rest frame
C       X(7), X(8)  : mu- in the W- rest frame
C       X(12)       : rapidity of the t tbar system (LUMINO)
C       X(13)       : m_tt (INT13 only)
C=======================================================================
C
C-----------------------------------------------------------------------
C     INT13(X): sigma [pb] integrand over X(1..13) (tot_bwblv.f).
C     MTT = exp(L), L uniform in [ln MIN, ln MAX] (COMMON/INV/);
C     Jacobian d(MTT)/dX(13) = (ln MAX - ln MIN) MTT. The power-law
C     mapping is kept as comments.
C     IDYN=1: mu_R = mu_F = MTT. GG and G set here belong to a local
C     COMMON/C1/, not to COUPL_QCD (coupl_sm.inc) used by the matrix
C     elements, so only the PDF scale changes.
C     MTT is passed to INT12 through COMMON/V12/.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT13 (X)
      IMPLICIT NONE
      INTEGER I,A
      DOUBLE PRECISION X(13),XX(12),Z,ZJ,JAC
      DOUBLE PRECISION INT12
      EXTERNAL INT12
      DOUBLE PRECISION RS,MTT
      COMMON/TT/MTT
      INTEGER IVAL,IBSS
      COMMON/V13/RS,IVAL,IBSS
      DOUBLE PRECISION VRS,VMTT
      INTEGER JVAL,JBSS,JJ
      COMMON/V12/VRS,VMTT,JVAL,JBSS,JJ
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
      DOUBLE COMPLEX GWF(2),GG(2)
      DOUBLE PRECISION G
      COMMON/C1/GWF,GG,G
      JJ = 0
      INT13 = 0D0
C.....
*      A = 3
*      Z = X(13)**A              !
*      ZJ = A * X(13)**(A-1)
*      MTT = (MAX-MIN)*Z + MIN
*      JAC = (MAX-MIN) * ZJ
C.....
      LMX = DLOG(MAX)
      LMN = DLOG(MIN)
      L = (LMX-LMN)*X(13) + LMN
      MTT = DEXP(L)
      JAC = (LMX-LMN) * MTT
C.....
      IF (IDYN.EQ.1) THEN
         MUR = MTT
         MUF = MUR
         ASR = ASQCD (MUR,INF)
         GS = DSQRT (4D0*PI*ASR)
         GG(1) = DCMPLX (GS,0D0)
         GG(2) = GG(1)
         G = -GS
      ENDIF
      DO I = 1,12
         XX(I) = X(I)
      ENDDO
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      INT13 = INT12(XX) * JAC
      CALL MKHIST (13,IVAL,IBSS,INT13,X)
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT12(X): d(sigma)/d(m_tt) [pb/GeV] at MTT = m_tt from
C     COMMON/V12/ (set by main_bwblv.f or INT13), integrand over
C     X(1..12). s = RS^2, tau = T = MTT^2/s; JAC = 1/s, with the
C     factor 2*MTT in SIGHAT = d(tau)/d(MTT). IPART = 0 sums q qbar
C     (JPART=1) and g g (JPART=2) explicitly; IPART = 1 or 2 selects
C     one channel. MKHIST (main program) fills the histograms.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT12 (X)
      IMPLICIT NONE
      DOUBLE PRECISION X(12),Z
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION JAC,LUM,SIGH
      DOUBLE PRECISION RS,MTT
      INTEGER IVAL,IBSS,II
      COMMON/V12/RS,MTT,IVAL,IBSS,II
      DOUBLE PRECISION S,SH,RSH,T
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      INTEGER JPART
      INT12 = 0D0
      S = RS**2
      T = (MTT/RS)**2
      SH = T * S
      RSH = DSQRT(SH)
      JAC = 1D0/S
      Z = X(12)                 !
      IF (IPART.EQ.0) THEN
         DO JPART = 1,2,1
            CALL SIGHAT (X,RSH,SIGH,JPART,ICLR,IBSS)
            CALL LUMINO (Z,T,JPART,LUM)
            INT12 = INT12 + SIGH * LUM * JAC
         ENDDO
      ELSE
         CALL SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
         CALL LUMINO (Z,T,IPART,LUM)
         INT12 = SIGH * LUM * JAC
      ENDIF
      CALL MKHIST (12,IVAL,IBSS,INT12,X)
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     SIGHAT: partonic part. Input X(1..11), RSH = sqrt(s-hat) [GeV],
C     IPART (1 q qbar, 2 g g), ICLR (g g: 0 singlet + octet,
C     1 singlet, 2 octet). Output SIGH = |M|^2 x flux x phase space
C     x 2*MTT, in pb/GeV.
C     Momenta P(0:3,1..7) = g1, g2, b, W+, bbar, mu-, vm~ (parton c.m.
C     frame, KINEMA5 in kinema5.f), stored in COMMON/MOM/.
C     Breit-Wigner mapping of a virtuality s (mass m, width w):
C       s = m^2 + m w tan(z), z uniform between the atan of the
C       limits; JAC_i = (ds/dX)/(2 pi), i.e. it includes the 1/(2 pi)
C       of ds/(2 pi).
C     MWCUT (COMMON/CUT/, 5 GeV in the main programs) is the lower
C     limit of the W- virtuality.
C-----------------------------------------------------------------------
      SUBROUTINE SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
      IMPLICIT NONE
      INTEGER IPART,ICLR,IBSS
      DOUBLE PRECISION X(12),RSH,SIGH,MTT
      DOUBLE PRECISION COS1,PHI1,COS2,PHI2,COS3,PHI3,COS4,PHI4
      DOUBLE PRECISION LMD1,LMD2,LMD3,LMD4,LAM
      DOUBLE PRECISION JAC,JAC1,JAC2,JAC3
      DOUBLE PRECISION SH,ME2,ME2S,ME2O
      DOUBLE PRECISION MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      DOUBLE PRECISION MR1,MR2,MR3
      COMMON/BWM/MR1,MR2,MR3
      DOUBLE PRECISION Z1,Z2,Z3,Z1MN,Z1MX,Z2MX,Z2MN,Z3MX,Z3MN
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      INTEGER ICUT
      DOUBLE PRECISION MWCUT
      COMMON/CUT/MWCUT,ICUT
      DOUBLE PRECISION P(0:3,7)
      COMMON/MOM/P
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141593D0)
      INTEGER IFLAG
      ME2 = 0D0
      SIGH = 0D0
      SH = RSH**2
C...  Resonance mass 1
C     Top: MR1 in [m_b+m_W, sqrt(s-hat)-m_b-MWCUT] from X(1)
      MIN1 = MB + MW
      MAX1 = RSH - MB - MWCUT
      Z1MN = DATAN2 (MIN1**2-MT**2, MT*GAMT)
      Z1MX = DATAN2 (MAX1**2-MT**2, MT*GAMT)
      Z1 = (Z1MX-Z1MN)*X(1) + Z1MN !
      MR1 = DSQRT( MT**2 + MT*GAMT*DTAN(Z1) )
      JAC1 = (Z1MX-Z1MN) * ((MR1**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (2D0*PI * MT*GAMT)
C...  Resonance mass 2
C     Antitop: MR2 in [m_b+MWCUT, sqrt(s-hat)-MR1] from X(2); a flat
C     mapping is used if the upper limit MAX2 is below 10 GeV
      MIN2 = MB + MWCUT
      MAX2 = RSH - MR1
      IF (MAX2.GT.10D0) THEN
         Z2MN = DATAN2 (MIN2**2-MT**2, MT*GAMT)
         Z2MX = DATAN2 (MAX2**2-MT**2, MT*GAMT)
         Z2 = (Z2MX-Z2MN)*X(2) + Z2MN !
         MR2 = DSQRT( MT**2 + MT*GAMT*DTAN(Z2) )
         JAC2 = (Z2MX-Z2MN) * ((MR2**2-MT**2)**2 + MT**2*GAMT**2)
     -        / (2D0*PI * MT*GAMT)
      ELSE
         MR2 = (MAX2-MIN2)*X(2) + MIN2
         JAC2 = MAX2-MIN2
      ENDIF
C...  Resonance mass 3
C     W-: MR3 in [MWCUT, MR2-m_b] from X(3); flat if MAX3 < 10 GeV
      MIN3 = MWCUT
      MAX3 = MR2 - MB
      IF (MAX3.GT.10D0) THEN
         Z3MN = DATAN2 (MIN3**2-MW**2, MW*GAMW)
         Z3MX = DATAN2 (MAX3**2-MW**2, MW*GAMW)
         Z3 = (Z3MX-Z3MN)*X(3) + Z3MN !
         MR3 = DSQRT( MW**2 + MW*GAMW*DTAN(Z3) )
         JAC3 = (Z3MX-Z3MN) * ((MR3**2-MW**2)**2 + MW**2*GAMW**2)
     -        / (2D0*PI * MW*GAMW)
      ELSE
         MR3 = (MAX3-MIN3)*X(3) + MIN3
         JAC3 = MAX3-MIN3
      ENDIF
C...
C     Angles: cos(theta) = 2X-1, phi = 2 pi X; with dOmega/(4 pi) in
C     the two-body phase space the Jacobian is 1.
      COS1 = 2D0*X(4) - 1D0     !
      COS2 = 2D0*X(5) - 1D0     !
      COS3 = 2D0*X(6) - 1D0     !
      COS4 = 2D0*X(7) - 1D0     !
      PHI4 = 2D0*PI*X(8)        !
      PHI3 = 2D0*PI*X(9)        !
      PHI2 = 2D0*PI*X(10)       !
      PHI1 = 2D0*PI*X(11)       !
      CALL KINEMA5 (RSH,MR1,MR2,MR3,
     -     COS1,PHI1,COS2,PHI2,COS3,PHI3,COS4,PHI4,P)
C...  Kinematical Cut
      IF (ICUT.EQ.1) THEN
         CALL KINECUT (P,IFLAG)
         IF (IFLAG.EQ.1) RETURN
      ENDIF
C.....
C     |M|^2, summed over helicities and colours and averaged over the
C     initial state; q qbar with the quark along +z
      IF ( IPART.EQ.1 ) THEN
         CALL QQBWBLV (P,ME2)
      ELSEIF ( IPART.EQ.2 ) THEN
         IF ( ICLR.EQ.0 ) THEN
            CALL GGBWBLV (P,ME2S,1)
            CALL GGBWBLV (P,ME2O,2)
            ME2 = ME2S + ME2O
         ELSE
            CALL GGBWBLV (P,ME2,ICLR)
         ENDIF
      ENDIF
C     Flux 1/(2 s-hat); two-body phase spaces lambda^(1/2)/(8 pi) for
C     t tbar, t -> b W+, tbar -> bbar W-* and W-* -> mu- vm~ (LMD4 = 1
C     for massless leptons). 2*MTT = d(s-hat)/d(MTT); 389429.57D3
C     converts GeV^-2 to pb.
      JAC  = JAC1 * JAC2 * JAC3
      LMD1 = LAM (1D0,MR1**2/SH,   MR2**2/SH    )
      LMD2 = LAM (1D0,MB**2/MR1**2,MW**2 /MR1**2)
      LMD3 = LAM (1D0,MB**2/MR2**2,MR3**2/MR2**2)
      LMD4 = 1D0
      MTT = RSH
      SIGH = ME2/(2D0*SH) * LMD1/(8D0*PI) * LMD2/(8D0*PI)
     -     * LMD3/(8D0*PI) * LMD4/(8D0*PI) * JAC
     -     * 2D0*MTT * 389429.57D3 ! pb/GeV
      RETURN
      END
C
C-----------------------------------------------------------------------
C     LUMINO: X (= X(12)) is mapped linearly to the t tbar rapidity
C     y in [ln(T)/2, -ln(T)/2]; x1 = sqrt(T) e^y, x2 = sqrt(T) e^-y,
C     dx1 dx2 = d(tau) dy. LUM = luminosity x DY, with QQ (q qbar +
C     qbar q, flavours d..b) for IPART=1 and GG for IPART=2 from
C     common/ptfl.f, at mu_F. x1, x2 go to COMMON/XBJ/ and xbk.
C-----------------------------------------------------------------------
      SUBROUTINE LUMINO (X,T,IPART,LUM)
      IMPLICIT NONE
      INTEGER IPART
      DOUBLE PRECISION X,T,LUM
      DOUBLE PRECISION YMIN,YMAX,DY,Y
      DOUBLE PRECISION X1,X2,QQ,GG,QG,GQ
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
      CALL PTFL (X1,X2,QQ,GG,QG,GQ,MUF**2,IPRC)
      IF (IPART.EQ.1) THEN
         LUM = QQ
      ELSEIF (IPART.EQ.2) THEN
         LUM = GG
      ENDIF
      JAC = DY
      LUM = LUM * JAC
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     KINECUT (ICUT=1). R2 (kin_functions.f) is Delta R squared in
C     (eta, phi). RBL1: b and mu-, RBL2: bbar and mu-; computed in the
C     parton c.m. frame (the boost PCM2HCM is commented out).
C     Active cut: IFLAG = 1 (rejected) if RBL2 > 0.2, i.e. only events
C     with the bbar close to the mu- are kept. ETL, ETB1, ETB2 unused.
C-----------------------------------------------------------------------
      SUBROUTINE KINECUT (P,IFLAG)
      IMPLICIT NONE
      INTEGER IFLAG
      DOUBLE PRECISION P(0:3,7), PH(0:3,7)
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MR1,MR2,MR3
      COMMON/BWM/MR1,MR2,MR3
      DOUBLE PRECISION RBL1,RBL2,ETB1,ETB2,ETL
      DOUBLE PRECISION R2,ETA
      EXTERNAL R2,ETA
      IFLAG = 0
*      CALL PCM2HCM (X1,X2,P,PH)
      RBL1 = R2(P(0,3),P(0,6))
      RBL2 = R2(P(0,5),P(0,6))
      ETL  = ETA(P(0,6))
      ETB1 = ETA(P(0,3))
      ETB2 = ETA(P(0,5))
*      IF (RBL1.LT.0.2D0) IFLAG = 1
*      IF (RBL2.LT.0.2D0 .OR. RBL2.GT.0.4D0) IFLAG = 1
      IF (RBL2.GT.0.2D0) IFLAG = 1
      RETURN
      END
