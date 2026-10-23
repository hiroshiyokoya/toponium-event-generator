C=======================================================================
C     RGAMT: top-quark width t -> b W [GeV] for a top of mass P [GeV]
C     (|V_tb| = 1, on-shell W).
C     Input: P, MW, MB [GeV]; GF = Fermi constant [GeV^-2];
C            AS = alpha_s (used only when MB = 0).
C       MB > 0          : LO width with the b mass (same formula as
C                         RUN in corr.f); AS is ignored
C       MB = 0, AS > 0  : massless b, LO term A0 plus the first-order
C                         QCD term AS/pi * C_F * A1 (C_F = 4/3)
C       MB = 0, AS <= 0 : LO width with a massless b
C     main_bwblv.f and main_blvblv.f set GAMT = RGAMT(MT,MW,GF,0D0,MB).
C     In the EG_*.f generators this call is commented out and the
C     width is TWIDTH from param_card_sm.dat.
C     The commented-out PROGRAM TEST below is a small test driver.
C=======================================================================
C
*      PROGRAM TEST
*      IMPLICIT NONE
*      DOUBLE PRECISION P,MT,MW,MB,GF,AS
*      DOUBLE PRECISION GAMT,RGAMT
*      EXTERNAL RGAMT
*      PARAMETER (MT=171D0,MW=80.4D0,MB=5D0,GF=1.16639D-5,AS=0.d0)
*      P = MT
*      GAMT = RGAMT(P,MW,GF,AS,MB)
*      WRITE (6,*) GAMT
*      STOP
*      END
C
      DOUBLE PRECISION FUNCTION RGAMT (P,MW,GF,AS,MB)
      IMPLICIT NONE
      DOUBLE PRECISION P,MW,MB,GF,AS
      DOUBLE PRECISION MW2,MB2,P2
      DOUBLE PRECISION X,Y,GT0,A0,A1
      DOUBLE PRECISION PI,SQ2,CF
      PARAMETER ( PI=3.141592654D0,SQ2=1.41421356D0,CF=1.3333333D0 )
C     x = m_W^2/P^2, y = m_b^2/P^2, GT0 = G_F P^3/(8 sqrt(2) pi)
      P2 = P**2
      MW2 = MW**2
      MB2 = MB**2
      X = MW2/P2
      Y = MB2/P2
      GT0 = GF/(8D0*SQ2*PI) * P**3
      IF (MB.GT.0D0) THEN
         A0 = (1D0-Y)**2 + X*(1D0+Y) - 2D0*X**2
         RGAMT = GT0 * A0 * DSQRT( 1D0 + X**2 + Y**2 - 2D0*(X+Y+X*Y) )
      ELSEIF (AS.GT.0D0) THEN
         A0 = 1D0 - 3D0*X**2 + 2D0*X**3
         A1 = 5D0/4D0 - PI**2/3D0 + 3D0/2D0*X +
     -        X**2 * (-6D0 + PI**2 + 3D0/2D0*DLOG(X)) +
     -        X**3 * (46D0/9D0 - 2D0/3D0*PI**2 - 2D0/3D0*DLOG(X))
         RGAMT = GT0 * (A0 + AS/PI*CF * A1)
      ELSE
         A0 = 1D0 - 3D0*X**2 + 2D0*X**3
         RGAMT = GT0 * A0
      ENDIF
      RETURN
      END
