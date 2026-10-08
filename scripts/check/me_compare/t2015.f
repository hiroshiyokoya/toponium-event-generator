C     |M|^2 of the 2015 stand-alone singlet matrix element (MEGGBLVBLV,
C     ICLR=1) summed over all helicities, for the points in points.dat
C     (MadEvent order, reordered here to g g b mu+ vm b~ mu- vm~).
C     alpha_s = 0.118.  usage: t2015 IGRN ITR INR
      PROGRAM T2015
      IMPLICIT NONE
      INCLUDE 'thr.inc'
      INCLUDE 'coupl_sm.inc'
      DOUBLE PRECISION P(0:3,8),Q(0:3,8),MEGGBLVBLV,SUM,PI,AS,GS
      EXTERNAL MEGGBLVBLV
      INTEGER NHEL(8),IC(8),I,J,K,IHEL,IOS,NPT,MAP(8)
      PARAMETER (PI=3.14159265358979D0, AS=0.118D0)
C     2015 slot <- 2010 slot:  b<-7, mu+<-3, vm<-4, b~<-8, mu-<-5, vm~<-6
      DATA MAP /1,2,7,3,4,8,5,6/
      CHARACTER*8 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ (ARG,*) IGRN
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ (ARG,*) ITR
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      READ (ARG,*) INR
      CALL SETPARA ('param_card_sm.dat',.TRUE.)
      GS = DSQRT(4D0*PI*AS)
      GG(1) = DCMPLX(GS,0D0)
      GG(2) = GG(1)
      G = -GS
      CALL READHEAD (1)
      DO I = 1, 8
         IC(I) = 1
      ENDDO
      OPEN (11,FILE='points.dat',STATUS='OLD')
      NPT = 0
  10  CONTINUE
      DO I = 1, 8
         READ (11,*,IOSTAT=IOS) (Q(J,I),J=0,3)
         IF (IOS.NE.0) GOTO 20
      ENDDO
      DO I = 1, 8
         DO J = 0, 3
            P(J,I) = Q(J,MAP(I))
         ENDDO
      ENDDO
      NPT = NPT + 1
      SUM = 0D0
      DO IHEL = 0, 255
         DO K = 1, 8
            NHEL(K) = 2*MOD(IHEL/2**(K-1),2) - 1
         ENDDO
         SUM = SUM + MEGGBLVBLV(P,NHEL,IC,1)
      ENDDO
      WRITE (6,'(I5,1PE20.12)') NPT, SUM
      GOTO 10
  20  CONTINUE
      END
