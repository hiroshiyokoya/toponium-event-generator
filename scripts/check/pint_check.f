C     Checks PINT (grid index for p) of the high-energy table against PP
C     (grid points) in standalone/common/readgrnep.f (docs/REVIEW.md).
C     For p in [PP(i), PP(i+1)) PINT should return i.
      PROGRAM PINTCK
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION MT,GAMT,ASB
      COMMON/GrnPAR/MT,GAMT,ASB
      DOUBLE PRECISION PP,P,EGEV(3)
      INTEGER PINT,M,IP,K,IE,NBAD
      EXTERNAL PP,PINT
      DATA EGEV/40D0,200D0,800D0/
      MT = 173D0
      M = (NP2-1)/2
      DO IE = 1, 3
         POS = DSQRT(MT*EGEV(IE))
         WRITE (6,'(A,F7.1,A,F9.4)') ' E =',EGEV(IE),'  p0 =',POS
         DO IP = M-2, M+2
            WRITE (6,'(A,I4,A,F12.6)') '   PP(',IP,') =',PP(IP,2)
         ENDDO
         NBAD = 0
         DO IP = 0, NP2-1
            DO K = 1, 9
               P = PP(IP,2) + (PP(IP+1,2)-PP(IP,2))*DBLE(K)/10D0
               IF (PINT(P,2).NE.IP) THEN
                  NBAD = NBAD + 1
                  IF (NBAD.LE.5) WRITE (6,'(A,F12.6,A,I4,A,I4)')
     .                 '   p =',P,'  PINT =',PINT(P,2),'  expected',IP
               ENDIF
            ENDDO
         ENDDO
         WRITE (6,'(A,I5,A)') '   ',NBAD,' of 909 test points misplaced'
      ENDDO
      END
