      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K,L
      INTEGER IEVT
C.....
      INTEGER IO
      DATA IO/60/
      DOUBLE PRECISION MT,MW,MB
      PARAMETER (MT=173D0, MW=80.4D0, MB=5D0)
C.....
      include 'nexternal.inc'
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(NEXTERNAL),ISTUP(NEXTERNAL),MOTHUP(2,NEXTERNAL),
     &     ICOLUP(2,NEXTERNAL),PUP(5,NEXTERNAL),VTIMUP(NEXTERNAL),
     &     SPINUP(NEXTERNAL)
C.....
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
      DOUBLE PRECISION MTTORG
C.....
      CALL LHEIN (-1)
      IEVT = 0
      DO I = 1, 100 000 000
         CALL LHEIN (IEVT)
         IF (IEVT.EQ.1) GOTO 100
         MTTORG = 4D0*SCALUP - 2D0*MT
         PXMISS = 0D0
         PYMISS = 0D0
         DO J = 1, NUP
            IF (IDUP(J).EQ.5) THEN
               DO K = 1, 3
                  PBT(K) = PUP(K,J)
               ENDDO
            ELSEIF (IDUP(J).EQ.-5) THEN
               DO K = 1, 3
                  PBX(K) = PUP(K,J)
               ENDDO
            ELSEIF (IDUP(J).EQ.-13) THEN
               DO K = 1, 3
                  PLP(K) = PUP(K,J)
               ENDDO
            ELSEIF (IDUP(J).EQ.13) THEN
               DO K = 1, 3
                  PLM(K) = PUP(K,J)
               ENDDO
            ELSEIF (ABS(IDUP(J)).EQ.14) THEN
               PXMISS = PXMISS + PUP(1,J)
               PYMISS = PYMISS + PUP(2,J)
            ENDIF
         ENDDO
         ICNT = 0
         NSOL = 0
         CALL DILEPSOL (PBT,PLP,PBX,PLM,PXMISS,PYMISS,MT,MW,MB)
         WRITE (IO,*) NSOL," SOLUTION"
         DO K = 1, NSOL
            WRITE (IO,*) K
            WRITE (IO,*) (SOL(K,L),L=1,3)
            WRITE (IO,*) (SOL(K,L),L=4,6)
            DO L = 1, 3
               PB1(L) = PBT(L)
               PB2(L) = PBX(L)
               PL1(L) = PLP(L)
               PL2(L) = PLM(L)
               PN1(L) = SOL(K,L)
               PN2(L) = SOL(K,L+3)
            ENDDO
            PB1(0) = ENERGY(PBT, MB)
            PB2(0) = ENERGY(PBX, MB)
            PL1(0) = ENERGY(PLP,0D0)
            PL2(0) = ENERGY(PLM,0D0)
            PN1(0) = ENERGY(PN1(1),0D0)
            PN2(0) = ENERGY(PN2(1),0D0)
            DO L = 0, 3
               PW1(L) = PL1(L) + PN1(L)
               PW2(L) = PL2(L) + PN2(L)
               PT1(L) = PB1(L) + PW1(L)
               PT2(L) = PB2(L) + PW2(L)
               PTT(L) = PT1(L) + PT2(L)
            ENDDO
            MW1 = MASCAL (PW1)
            MW2 = MASCAL (PW2)
            MT1 = MASCAL (PT1)
            MT2 = MASCAL (PT2)
            MTT = MASCAL (PTT)
            WRITE (IO,*) MW1, MW2, MT1, MT2, MTT
            ICNT = ICNT + 1
            MTTSOL(ICNT) = MTT
         ENDDO
         N = NSOL
         NSOL = 0
         CALL DILEPSOL (PBT,PLM,PBX,PLP,PXMISS,PYMISS,MT,MW,MB)
         WRITE (IO,*) NSOL," COMB"
         DO K = 1, NSOL
            WRITE (IO,*) K
            WRITE (IO,*) (SOL(K,L),L=1,3)
            WRITE (IO,*) (SOL(K,L),L=4,6)
            DO L = 1, 3
               PB1(L) = PBT(L)
               PB2(L) = PBX(L)
               PL1(L) = PLM(L)
               PL2(L) = PLP(L)
               PN1(L) = SOL(K,L)
               PN2(L) = SOL(K,L+3)
            ENDDO
            PB1(0) = ENERGY(PBT, MB)
            PB2(0) = ENERGY(PBX, MB)
            PL1(0) = ENERGY(PLM,0D0)
            PL2(0) = ENERGY(PLP,0D0)
            PN1(0) = ENERGY(PN1(1),0D0)
            PN2(0) = ENERGY(PN2(1),0D0)
            DO L = 0, 3
               PW1(L) = PL1(L) + PN1(L)
               PW2(L) = PL2(L) + PN2(L)
               PT1(L) = PB1(L) + PW1(L)
               PT2(L) = PB2(L) + PW2(L)
               PTT(L) = PT1(L) + PT2(L)
            ENDDO
            MW1 = MASCAL (PW1)
            MW2 = MASCAL (PW2)
            MT1 = MASCAL (PT1)
            MT2 = MASCAL (PT2)
            MTT = MASCAL (PTT)
            WRITE (IO,*) MW1, MW2, MT1, MT2, MTT
            ICNT = ICNT + 1
            MTTSOL(ICNT) = MTT
         ENDDO
         WRITE (1,*) N,NSOL
         IF (ICNT.GT.0) THEN
            MEAN = 0D0
            MIN = 1D6
            MAX = 0D0
            DO K = 1,ICNT
               MEAN = MEAN + MTTSOL(K)
               IF (MTTSOL(K).LT.MIN) MIN = MTTSOL(K)
               IF (MTTSOL(K).GT.MAX) MAX = MTTSOL(K)
               WRITE (10,*) MTTSOL(K)
            ENDDO
            MEAN = MEAN / DBLE(ICNT)
            WRITE (2,*) MTTORG, MEAN, MIN, MAX
         ELSE
            WRITE (3,*) MTTORG
         ENDIF
         WRITE (4,*) MTTORG
      ENDDO
 100  CONTINUE
      CALL LHEIN (1)
      STOP
      END
