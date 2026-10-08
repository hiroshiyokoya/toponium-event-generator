      SUBROUTINE DILEP (ID)
      IMPLICIT NONE
      INTEGER ID,I,J,K
      DOUBLE PRECISION MT,MW,MB
      PARAMETER (MT=173D0, MW=80.4D0, MB=5D0)
      DOUBLE PRECISION MR1,MR2,MR3,MR4
      COMMON/BWM/MR1,MR2,MR3,MR4
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION P(0:3,8),PH(0:3,8)
      COMMON/MOM/P
      DOUBLE PRECISION MTTORG
      COMMON/TT/MTTORG
      DOUBLE PRECISION PBT(3),PBX(3),PLP(3),PLM(3)
      DOUBLE PRECISION PXMISS,PYMISS
      DOUBLE PRECISION SOL(16,6)
      INTEGER N,NSOL
      COMMON/SOL/SOL,NSOL
      DOUBLE PRECISION PB1(0:3),PB2(0:3),PL1(0:3),PL2(0:3)
      DOUBLE PRECISION PN1(0:3),PN2(0:3),PW1(0:3),PW2(0:3)
      DOUBLE PRECISION PT1(0:3),PT2(0:3),PTT(0:3)
      DOUBLE PRECISION MW1,MW2,MT1,MT2,MTT
      DOUBLE PRECISION ENERGY,MASCAL
      EXTERNAL ENERGY,MASCAL
      INTEGER ICNT
      DOUBLE PRECISION MTTSOL(8),MEAN,MIN,MAX
      CALL PCM2HCM (X1,X2,P,PH)
      DO I = 1, 3
         PBT(I) = PH(I,3)
         PLP(I) = PH(I,4)
         PBX(I) = PH(I,6)
         PLM(I) = PH(I,7)
      ENDDO
      PXMISS = PH(1,5) + PH(1,8)
      PYMISS = PH(2,5) + PH(2,8)
      WRITE (6,*) "ANS"
      WRITE (6,*) (PH(I,5),I=1,3)
      WRITE (6,*) (PH(I,8),I=1,3)
      WRITE (6,*) MR3,MR4,MR1,MR2,MTT
      ICNT = 0
      NSOL = 0
      CALL DILEPSOL (PBT,PLP,PBX,PLM,PXMISS,PYMISS,MT,MW,MB)
      WRITE (6,*) NSOL," SOLUTION"
      DO I = 1, NSOL
         WRITE (6,*) I
         WRITE (6,*) (SOL(I,J),J=1,3)
         WRITE (6,*) (SOL(I,J),J=4,6)
         DO J = 1, 3
            PB1(J) = PBT(J)
            PB2(J) = PBX(J)
            PL1(J) = PLP(J)
            PL2(J) = PLM(J)
            PN1(J) = SOL(I,J)
            PN2(J) = SOL(I,J+3)
         ENDDO
         PB1(0) = ENERGY(PBT, MB)
         PB2(0) = ENERGY(PBX, MB)
         PL1(0) = ENERGY(PLP,0D0)
         PL2(0) = ENERGY(PLM,0D0)
         PN1(0) = ENERGY(PN1(1),0D0)
         PN2(0) = ENERGY(PN2(1),0D0)
         DO J = 0, 3
            PW1(J) = PL1(J) + PN1(J)
            PW2(J) = PL2(J) + PN2(J)
            PT1(J) = PB1(J) + PW1(J)
            PT2(J) = PB2(J) + PW2(J)
            PTT(J) = PT1(J) + PT2(J)
         ENDDO
         MW1 = MASCAL (PW1)
         MW2 = MASCAL (PW2)
         MT1 = MASCAL (PT1)
         MT2 = MASCAL (PT2)
         MTT = MASCAL (PTT)
         WRITE (6,*) MW1, MW2, MT1, MT2, MTT
         MTTSOL(ICNT+1) = MTT
         ICNT = ICNT + 1
      ENDDO
      N = NSOL
      NSOL = 0
      CALL DILEPSOL (PBT,PLM,PBX,PLP,PXMISS,PYMISS,MT,MW,MB)
      WRITE (6,*) NSOL," COMB"
      DO I = 1, NSOL
         WRITE (6,*) I
         WRITE (6,*) (SOL(I,J),J=1,3)
         WRITE (6,*) (SOL(I,J),J=4,6)
         DO J = 1, 3
            PB1(J) = PBT(J)
            PB2(J) = PBX(J)
            PL1(J) = PLM(J)
            PL2(J) = PLP(J)
            PN1(J) = SOL(I,J)
            PN2(J) = SOL(I,J+3)
         ENDDO
         PB1(0) = ENERGY(PBT, MB)
         PB2(0) = ENERGY(PBX, MB)
         PL1(0) = ENERGY(PLM,0D0)
         PL2(0) = ENERGY(PLP,0D0)
         PN1(0) = ENERGY(PN1(1),0D0)
         PN2(0) = ENERGY(PN2(1),0D0)
         DO J = 0, 3
            PW1(J) = PL1(J) + PN1(J)
            PW2(J) = PL2(J) + PN2(J)
            PT1(J) = PB1(J) + PW1(J)
            PT2(J) = PB2(J) + PW2(J)
            PTT(J) = PT1(J) + PT2(J)
         ENDDO
         MW1 = MASCAL (PW1)
         MW2 = MASCAL (PW2)
         MT1 = MASCAL (PT1)
         MT2 = MASCAL (PT2)
         MTT = MASCAL (PTT)
         WRITE (6,*) MW1, MW2, MT1, MT2, MTT
         MTTSOL(ICNT+1) = MTT
         ICNT = ICNT + 1
      ENDDO
      WRITE (1,*) N,NSOL
      IF (ICNT.GT.0) THEN
         MEAN = 0D0
         MIN = 1D6
         MAX = 0D0
         DO I = 1,ICNT
            MEAN = MEAN + MTTSOL(I)
            IF (MTTSOL(I).LT.MIN) MIN = MTTSOL(I)
            IF (MTTSOL(I).GT.MAX) MAX = MTTSOL(I)
         ENDDO
         MEAN = MEAN / DBLE(ICNT)
         WRITE (2,*) MTTORG, MEAN, MIN, MAX
      ELSE
         WRITE (3,*) MTTORG
      ENDIF
      WRITE (4,*) MTTORG
      RETURN
      END
