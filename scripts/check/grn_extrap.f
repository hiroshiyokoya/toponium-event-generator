C     Values of the Green-function correction returned by the 2015 table
C     reader (standalone/common/readgrnep.f) inside and outside the table
C     range (issue #9). Run in a directory that contains grnep*.tbl.
      PROGRAM GRNEXT
      IMPLICIT NONE
      DOUBLE COMPLEX GAMMA
      DOUBLE PRECISION E(10),P(4)
      INTEGER I,J,ICLR
      DATA E/-100D0,-60D0,-40D0,0D0,10D0,100D0,800D0,830D0,1000D0,
     .       2000D0/
      DATA P/0D0,20D0,100D0,300D0/
      DO ICLR = 1, 2
         CALL READHEAD (ICLR)
         WRITE (6,'(A,I2)') ' colour ',ICLR
         DO I = 1, 10
            DO J = 1, 4
               CALL READGRN (GAMMA,E(I),P(J),ICLR)
               WRITE (6,'(2F9.1,2F12.5)') E(I),P(J),GAMMA
            ENDDO
         ENDDO
      ENDDO
      END
