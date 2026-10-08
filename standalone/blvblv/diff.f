      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K,N
      DOUBLE PRECISION MH(1 000 000)
      INTEGER H(1 000 000), DH(1 000 000)
      OPEN(10,FILE="fort.11")
      DO I=1,1 000 000
         READ(10,*,END=100) MH(I), H(I)
      ENDDO
 100  CONTINUE
      CLOSE(10)
      N=I-1
      WRITE (6,*) N
      DO I=2,N-2
         WRITE (12,*) MH(I), (H(I+1)-H(I-1))/(MH(I+1)-MH(I-1))
     -        - (MH(I+1)-MH(I))**2/6D0*H(I)
      ENDDO
      STOP
      END
