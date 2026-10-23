C=======================================================================
C     GRNMOD: bound-state correction factor Gamma = G/G0 at one
C     phase-space point (Sumino & Yokoya, arXiv:1007.0075, Eq. (3.1);
C     docs/green-function-tables.md, Sec. 3). Called from GGTTRES
C     (ggtt.f, g g channel) and from the q qbar matrix elements.
C
C     Input
C       ST, SX : squared invariant masses of the top and the antitop,
C                s_t = (p_b+p_W+)^2, s_tbar = (p_bbar+p_W-)^2 [GeV^2]
C       SH     : partonic s-hat [GeV^2]
C       MT     : top mass [GeV]
C       ICLR   : colour state of the t tbar pair, 1 singlet, 2 octet
C       ITHR   : switch; the callers pass IGRN. ITHR >= 1 applies the
C                correction, otherwise GAMMA = 1.
C     Output
C       GAMMA  : complex correction factor (dimensionless)
C
C     Kinematics, as coded below:
C       E = s-hat/(4 m_t) - m_t                              [GeV]
C       p = sqrt(s-hat)/2 * lambda^(1/2)(1, s_t/s-hat, s_tbar/s-hat)
C     p [GeV] is the top momentum in the t tbar rest frame; near
C     threshold E ~ sqrt(s-hat) - 2 m_t. Gamma(E,p) is read from the
C     tables by bilinear interpolation (READGRN -> HOKAN2D,
C     readgrnep.f).
C=======================================================================
      SUBROUTINE GRNMOD (GAMMA,ST,SX,SH,MT,ICLR,ITHR)
      IMPLICIT NONE
      DOUBLE COMPLEX GAMMA
      DOUBLE PRECISION MT
      INTEGER ICLR,ITHR
      DOUBLE PRECISION E,P
      DOUBLE PRECISION ST,SX,SH,SBT
      DOUBLE PRECISION LAM,A,B,C
C     Statement function: square root of the Kallen function
C     lambda(A,B,C). The MAX(0,...) acts on the result of DSQRT, so it
C     does not guard against a negative argument.
      LAM(A,B,C) = MAX(0D0,DSQRT(A**2+B**2+C**2-2D0*(A*B+B*C+C*A)))
      IF (ITHR.GE.1) THEN
         SBT = 1D0 - 4D0*MT**2/SH
C.....
C     E [GeV]. Active definition: s-hat/(4 m_t) - m_t (= p_on^2/m_t,
C     p_on = on-shell top momentum). Older choices kept as comments:
C     m_t beta^2 and sqrt(s-hat) - 2 m_t.
*         E = MT * SBT
*         E = DSQRT(SH) - 2D0*MT
         E = SH/(4D0*MT) - MT
C.....
C     p [GeV]: top momentum in the t tbar rest frame
         P = DSQRT(SH)/2D0 * LAM(1D0,ST/SH,SX/SH)
         CALL READGRN (GAMMA,E,P,ICLR)
      ELSE
         GAMMA = (1D0,0D0)
      ENDIF
      RETURN
      END
