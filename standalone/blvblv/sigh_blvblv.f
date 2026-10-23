C=======================================================================
C     sigh_blvblv.f : integrand for p p -> b mu+ vm bbar mu- vm~ (off-
C     shell t, tbar, W+, W-) with the bound-state correction. Older
C     version, used by main_blvblv.f (d(sigma)/d(m_tt) at fixed m_tt,
C     integrand INT15) and tot_blvblv.f (sigma, integrand INT16);
C     EG_blvblv.f uses sigh_blvblvE.f instead. Helicities and colours
C     are summed inside the matrix elements. For q qbar the luminosity
C     QQ of common/ptfl.f (both beam orientations) multiplies the
C     matrix element with the quark along +z.
C
C     Integration variables, all in [0,1]:
C       X(1), X(2)  : top and antitop virtualities (Breit-Wigner)
C       X(3), X(4)  : W+ and W- virtualities (Breit-Wigner)
C       X(5), X(14) : cos(theta), phi of the top (parton c.m. frame)
C       X(6), X(13) : b in the top rest frame
C       X(7), X(12) : bbar in the antitop rest frame
C       X(8), X(11) : mu+ in the W+ rest frame
C       X(9), X(10) : mu- in the W- rest frame
C       X(15)       : rapidity of the t tbar system (LUMINO)
C       X(16)       : m_tt (INT16 only)
C=======================================================================
C
C-----------------------------------------------------------------------
C     INT16(X): sigma [pb] integrand over X(1..16) (tot_blvblv.f).
C     MTT = exp(L), L uniform in [ln MIN, ln MAX] (COMMON/INV/);
C     Jacobian d(MTT)/dX(16) = (ln MAX - ln MIN) MTT. The power-law
C     mapping is kept as comments. IDYN=1: mu_R = mu_F = MTT, g_s in
C     coupl_sm.inc reset. MTT goes to INT15 through COMMON/V15/.
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INT16 (X)
      IMPLICIT NONE
      INTEGER I,A
      DOUBLE PRECISION X(16),XX(15),Z,ZJ,JAC
      DOUBLE PRECISION INT15
      EXTERNAL INT15
      DOUBLE PRECISION RS,MTT
      COMMON/TT/MTT
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
         MUR = MTT
         MUF = MUR
         ASR = ASQCD (MUR,INF)
         GS = DSQRT (4D0*PI*ASR)
         GG(1) = DCMPLX (GS,0D0)
         GG(2) = GG(1)
         G = -GS
      ENDIF
      DO I = 1, 15
         XX(I) = X(I)
      ENDDO
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      INT16 = INT15(XX) * JAC
      CALL MKHIST (16,IVAL,IBSS,INT16,X)
      RETURN
      END
C
C-----------------------------------------------------------------------
C     INT15(X): d(sigma)/d(m_tt) [pb/GeV] at MTT = m_tt from
C     COMMON/V15/ (set by main_blvblv.f or INT16), integrand over
C     X(1..15). s = RS^2, tau = T = MTT^2/s, JAC = 1/s. IPART = 0
C     sums q qbar (JPART=1) and g g (JPART=2) explicitly.
C     MKHIST (main program) fills the histograms.
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
      INTEGER JPART
      INT15 = 0D0
      S = RS**2
      T = (MTT/RS)**2
      SH = T * S
      RSH = DSQRT(SH)
      JAC = 1D0/S
      Z = X(15)                 !
      IF (IPART.EQ.0) THEN
         DO JPART = 1,2,1
            CALL SIGHAT (X,RSH,SIGH,JPART,ICLR,IBSS)
            CALL LUMINO (Z,T,JPART,LUM)
            INT15 = INT15 + SIGH * LUM * JAC
         ENDDO
      ELSE
         CALL SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
         CALL LUMINO (Z,T,IPART,LUM)
         INT15 = SIGH * LUM * JAC
      ENDIF
C     2*MTT/s = d(tau)/d(MTT); 389429.57D3 converts GeV^-2 to pb
      INT15 = INT15 * 2D0*MTT * 389429.57D3 ! pb/GeV
      CALL MKHIST (15,IVAL,IBSS,INT15,X)
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     SIGHAT: partonic part. Input X(1..14), RSH = sqrt(s-hat) [GeV],
C     IPART (1 q qbar, 2 g g), ICLR (g g: 0 singlet + octet,
C     1 singlet, 2 octet). Output SIGH = |M|^2 x flux x phase space
C     [GeV^-2].
C     Momenta P(0:3,1..8) = g1, g2, b, mu+, vm, bbar, mu-, vm~ (parton
C     c.m. frame, KINEMA6 in kinema6.f), stored in COMMON/MOM/.
C     Breit-Wigner mapping of a virtuality s (mass m, width w):
C       s = m^2 + m w tan(z), z uniform between the atan of the
C       limits; JAC_i = (ds/dX)/(2 pi), i.e. it includes the 1/(2 pi)
C       of ds/(2 pi). The flat mappings are kept as comments.
C     MWCUT (COMMON/CUT/, 5 GeV in the main programs) is the lower
C     limit of the W virtualities.
C-----------------------------------------------------------------------
      SUBROUTINE SIGHAT (X,RSH,SIGH,IPART,ICLR,IBSS)
      IMPLICIT NONE
      INTEGER IPART,ICLR,IBSS
      DOUBLE PRECISION X(15),RSH,SIGH,MTT
      DOUBLE PRECISION COS1,PHI1,COS2,PHI2,COS3,PHI3,COS4,PHI4,COS5,PHI5
      DOUBLE PRECISION LMD1,LMD2,LMD3,LMD4,LMD5,LAM
      DOUBLE PRECISION JAC,JAC1,JAC2,JAC3,JAC4
      DOUBLE PRECISION SH,ME2,ME2S,ME2O
      DOUBLE PRECISION MIN1,MAX1,MIN2,MAX2,MIN3,MAX3,MIN4,MAX4
      DOUBLE PRECISION MR1,MR2,MR3,MR4
      COMMON/BWM/MR1,MR2,MR3,MR4
      DOUBLE PRECISION Z1,Z2,Z3,Z4
      DOUBLE PRECISION Z1MN,Z1MX,Z2MX,Z2MN,Z3MX,Z3MN,Z4MX,Z4MN
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      INTEGER ICUT
      DOUBLE PRECISION MWCUT
      COMMON/CUT/MWCUT,ICUT
      DOUBLE PRECISION P(0:3,8)
      COMMON/MOM/P
      COMMON/TT/MTT
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141593D0)
      INTEGER IFLAG
      ME2 = 0D0
      SIGH = 0D0
      SH = RSH**2
C...  Resonance mass 1
C     Top: MR1 in [m_b+MWCUT, sqrt(s-hat)-m_b-MWCUT] from X(1)
      MIN1 = MB + MWCUT
      MAX1 = RSH - MB - MWCUT
      Z1MN = DATAN2 (MIN1**2-MT**2, MT*GAMT)
      Z1MX = DATAN2 (MAX1**2-MT**2, MT*GAMT)
      Z1 = (Z1MX-Z1MN)*X(1) + Z1MN !
      MR1 = DSQRT( MT**2 + MT*GAMT*DTAN(Z1) )
      JAC1 = (Z1MX-Z1MN) * ((MR1**2-MT**2)**2 + MT**2*GAMT**2)
     -     / (2D0*PI * MT*GAMT)
*      MR1 = (MAX1-MIN1)*X(1) + MIN1
*      JAC1 = MAX1-MIN1
C...  Resonance mass 2
C     Antitop: MR2 in [m_b+MWCUT, sqrt(s-hat)-MR1] from X(2)
      MIN2 = MB + MWCUT
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
C     W+: MR3 in [MWCUT, MR1-m_b] from X(3)
      MIN3 = MWCUT
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
C     W-: MR4 in [MWCUT, MR2-m_b] from X(4)
      MIN4 = MWCUT
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
C     3: bbar in tbar, 4: mu+ in W+, 5: mu- in W- rest frame.
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
C...  Kinematical Cut
      IF (ICUT.EQ.1) THEN
         CALL KINECUT (P,IFLAG)
         IF (IFLAG.EQ.1) RETURN
      ENDIF
C.....
C     |M|^2, summed over helicities and colours and averaged over the
C     initial state; q qbar with the quark along +z
      IF ( IPART.EQ.1 ) THEN
         CALL QQBLVBLV (P,ME2)
      ELSEIF ( IPART.EQ.2 ) THEN
         IF ( ICLR.EQ.0 ) THEN
            CALL GGBLVBLV (P,ME2S,1)
            CALL GGBLVBLV (P,ME2O,2)
            ME2 = ME2S + ME2O
         ELSE
            CALL GGBLVBLV (P,ME2,ICLR)
         ENDIF
      ENDIF
C     Flux 1/(2 s-hat); two-body phase spaces lambda^(1/2)/(8 pi) for
C     t tbar, t -> b W+*, tbar -> bbar W-*, and W -> l nu (LMD4 = LMD5
C     = 1 for massless leptons). LAM is in kinema6.f.
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
C     LUMINO: X (= X(15)) is mapped linearly to the t tbar rapidity
C     y in [ln(T)/2, -ln(T)/2]; x1 = sqrt(T) e^y, x2 = sqrt(T) e^-y,
C     dx1 dx2 = d(tau) dy. LUM = luminosity x DY, with QQ (q qbar +
C     qbar q, flavours d..b) for IPART=1 and GG for IPART=2 from
C     common/ptfl.f, at mu_F; negative values are set to 0.
C     x1, x2 go to COMMON/XBJ/ and xbk (run.inc).
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
      IF (LUM.LT.0D0) LUM = 0D0
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C     KINECUT (ICUT=1): all cuts are commented out, IFLAG stays 0.
C     The particle indices are those of sigh_bwblv.f; with the
C     ordering of KINEMA6, P(0,5) is the vm and P(0,6) the bbar.
C     R2 (kin_functions.f) is Delta R squared.
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
