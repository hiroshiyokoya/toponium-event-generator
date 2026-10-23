C=======================================================================
C     CORR: factor COR for ITR=1 (used in GGTTRES, ggtt.f, and in the
C     q qbar matrix elements me_qq*E.f); docs/green-function-tables.md:
C
C       COR = prod_{i=t,tbar} [ m_t G_t / (sqrt(s_i) G_t(s_i)) ]^(1/2)
C
C     G_t(s) is the LO width of a top quark of virtuality s (function
C     RUN below), G_t = GT the width used in the propagators.
C
C     Input: S1, S2 = s_t, s_tbar [GeV^2]; MT, GT = top mass and width
C            [GeV]; GF = Fermi constant [GeV^-2]; MW, MB [GeV].
C     There is no guard for sqrt(s_i) <= m_W + m_b, where RUN vanishes
C     (or is undefined). The EG_*.f generators keep the top
C     virtualities above m_b + m_W (bwbw) or m_b + MWCUT1 =
C     m_b + m_W + 0.1 GeV (blvblv, blvbjj, bjjbjj).
C=======================================================================
      DOUBLE PRECISION FUNCTION CORR (S1,S2,MT,GT,GF,MW,MB)
      IMPLICIT NONE
      DOUBLE PRECISION S1,S2,MT,GT,GF,MW,MB
      DOUBLE PRECISION RUN,GT1,GT2
      EXTERNAL RUN
      GT1 = RUN(S1,MW,MB,GF)
      GT2 = RUN(S2,MW,MB,GF)
      CORR = DSQRT( MT*GT/(DSQRT(S1)*GT1) * MT*GT/(DSQRT(S2)*GT2) )
      RETURN
      END
C
C     RUN: LO width t -> b W [GeV] of a top quark of virtuality
C     P2 [GeV^2] (|V_tb| = 1, on-shell W):
C       G = G_F P^3/(8 sqrt(2) pi) * A0 * lambda^(1/2)(1,x,y),
C       A0 = (1-y)^2 + x(1+y) - 2x^2,  x = m_W^2/P2,  y = m_b^2/P2,
C     with P^3 = P2**1.5. Same formula as RGAMT (gamt.f) for MB > 0.
      DOUBLE PRECISION FUNCTION RUN (P2,MW,MB,GF)
      IMPLICIT NONE
      DOUBLE PRECISION P2,MW,MB,GF
      DOUBLE PRECISION X,Y,GT0,A0
      DOUBLE PRECISION PI,SQ2
      PARAMETER ( PI=3.141592654D0,SQ2=1.41421356D0)
      X = MW**2/P2
      Y = MB**2/P2
      GT0 = GF/(8D0*SQ2*PI) * P2**1.5
      A0 = (1D0-Y)**2 + X*(1D0+Y) - 2D0*X**2
      RUN = GT0 * A0 * DSQRT( 1D0 + X**2 + Y**2 - 2D0*(X+Y+X*Y) )
      RETURN
      END
