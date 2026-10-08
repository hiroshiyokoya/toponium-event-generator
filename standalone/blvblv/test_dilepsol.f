      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592654D0)
      DOUBLE PRECISION MT,MW,MB
      DOUBLE PRECISION PBT(3),PBX(3),PLP(3),PLM(3)
      DOUBLE PRECISION PXMISS,PYMISS
      DOUBLE PRECISION SOL(16,6)
      INTEGER NSOL
      COMMON/SOL/SOL,NSOL
      DOUBLE PRECISION PB1(0:3),PB2(0:3),PL1(0:3),PL2(0:3)
      DOUBLE PRECISION PN1(0:3),PN2(0:3),PW1(0:3),PW2(0:3)
      DOUBLE PRECISION PT1(0:3),PT2(0:3),PTT(0:3)
      DOUBLE PRECISION MW1,MW2,MT1,MT2,MTT
      DOUBLE PRECISION ENERGY,MASCAL
      EXTERNAL ENERGY,MASCAL
      PARAMETER ( MT=173D0, MW=80.4D0, MB=5D0 )
      DATA PBT/ -24.17D0,  39.08D0,   9.24D0/
      DATA PBX/ -99.28D0,-108.68D0,-177.50D0/
      DATA PLP/ 103.16D0, 134.10D0,-117.90D0/
      DATA PLM/  -2.63D0, -19.65D0, -64.19D0/
      DATA PXMISS,PYMISS/ 22.92D0,-44.84D0/
      NSOL = 0
      CALL DILEPSOL (PBT,PLP,PBX,PLM,PXMISS,PYMISS,MT,MW,MB)
      WRITE (6,*) "SOLVED", NSOL
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
      ENDDO
      NSOL = 0
      CALL DILEPSOL (PBT,PLM,PBX,PLP,PXMISS,PYMISS,MT,MW,MB)
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
      ENDDO
      STOP
      END
