C     |M|^2 of the 2010 MadEvent singlet subprocess (P1, g g -> mu+ vm
C     mu- vm~ b b~) summed over all helicities, for the points in
C     points.dat (MadEvent order). alpha_s = 0.118. Switches come from the
C     compile-time thr.inc of the process directory.
      PROGRAM T2010
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'coupl.inc'
      DOUBLE PRECISION P(0:3,NEXTERNAL),MATRIX,SUM,PI,AS
      EXTERNAL MATRIX
      INTEGER NHEL(NEXTERNAL),IC(NEXTERNAL),I,J,K,IHEL,IOS,NPT
      PARAMETER (PI=3.14159265358979D0, AS=0.118D0)
      CALL SETPARA ('param_card.dat',.TRUE.)
      G = DSQRT(4D0*PI*AS)
      GG(1) = -G
      GG(2) = -G
      DO I = 1, NEXTERNAL
         IC(I) = 1
      ENDDO
      OPEN (11,FILE='points.dat',STATUS='OLD')
      NPT = 0
  10  CONTINUE
      DO I = 1, NEXTERNAL
         READ (11,*,IOSTAT=IOS) (P(J,I),J=0,3)
         IF (IOS.NE.0) GOTO 20
      ENDDO
      NPT = NPT + 1
      SUM = 0D0
      DO IHEL = 0, 2**NEXTERNAL-1
         DO K = 1, NEXTERNAL
            NHEL(K) = 2*MOD(IHEL/2**(K-1),2) - 1
         ENDDO
         SUM = SUM + MATRIX(P,NHEL,IC)
      ENDDO
      WRITE (6,'(I5,1PE20.12)') NPT, SUM
      GOTO 10
  20  CONTINUE
      END
