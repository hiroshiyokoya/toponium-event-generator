C     Code to enumerate all the helicty combination,
C     and write the list in Madgraph format.
C     In 'helcomb.in', list the number of helicity followed by all the
C     values of helicity in each line for each particle
C     in the same ordering as the Madgraph numbering.
C
      PROGRAM MAIN
      INTEGER I,N
      INTEGER NCOMB,NEXTERNAL
      COMMON/COMB/NCOMB,NEXTERNAL
      NCOMB = 1
      OPEN(10,FILE="helcomb.in")
      DO I = 1, 10
         READ(10,'(I1)',END=1) N
         IF (N.EQ.0) GOTO 1
         NCOMB = NCOMB * N
         NEXTERNAL = I
      ENDDO
 1    CONTINUE
      CLOSE(10)
      WRITE (6,*) 'NEXTERNAL =',NEXTERNAL,', NCOMB= ',NCOMB
      CALL HELCOMB
      STOP
      END
C
      SUBROUTINE HELCOMB
      IMPLICIT NONE
      INTEGER I,J,K,IO
      INTEGER NCOMB,NEXTERNAL
      COMMON/COMB/NCOMB,NEXTERNAL
      INTEGER NH(NEXTERNAL), HL(NEXTERNAL,0:5)
      INTEGER NHEL(NEXTERNAL,NCOMB)
      INTEGER FLOAT,FMIN,DEG(NEXTERNAL)
      CHARACTER FMT1*48,FMT2*17,FMT3*4,SL*1,FMT*80
      DATA SL/'/'/
      FMT1 = '(''      DATA (NHEL(IHEL,'',I3,''), IHEL=1,'',I2,'')'
      FMT2 = ''',00(I2,'',''),I2,'''
      FMT3 = ''')'
      WRITE (FMT2(3:4),'(I2)') NEXTERNAL-1
      FMT = FMT1//SL//FMT2//SL//FMT3
      OPEN(10,FILE="helcomb.in")
      DO I = 1, NEXTERNAL
         READ (10,*) NH(I)
         BACKSPACE (10)
         READ (10,*) K, (HL(I,J),J=0,NH(I)-1)
         DEG(I) = 0
      ENDDO
      CLOSE(10)
C
      FLOAT = NEXTERNAL
      DO K = 1, NCOMB
         DO J = 1, NEXTERNAL
            NHEL(J,K) = HL(J,DEG(J))
         ENDDO
         WRITE (6,FMT) K, NEXTERNAL, (NHEL(J,K),J=1,NEXTERNAL)
C.....
         DO J = FLOAT, 1, -1
            IF (DEG(FLOAT).EQ.NH(FLOAT)-1) THEN
               DO I = FLOAT, NEXTERNAL
                  DEG(I) = 0
               ENDDO
               FLOAT = FLOAT - 1
            ENDIF
         ENDDO
         DEG(FLOAT) = DEG(FLOAT) + 1
         FLOAT = NEXTERNAL
      ENDDO
C...
      RETURN
      END
