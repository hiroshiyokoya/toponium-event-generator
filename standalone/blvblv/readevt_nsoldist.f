      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,J,K,IEVT
C.....
      DOUBLE PRECISION MT,MW,MB
      PARAMETER (MT=173D0, MW=80.4D0, MB=5D0)
      INTEGER BIN
      PARAMETER (BIN=100)
      DOUBLE PRECISION HMIN,HMAX,HSTP
      PARAMETER (HMIN=160D0,HMAX=186D0)
      INTEGER H(BIN)
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
      INTEGER NSOL
      COMMON/SOL/SOL,NSOL
      DOUBLE PRECISION TRY
C.....
      HSTP = (HMAX-HMIN)/BIN
      CALL LHEIN (-1)
      IEVT = 0
      DO J = 1,BIN
         H(J) = 0
      ENDDO
      DO I = 1, 100 00!0 000
         CALL LHEIN (IEVT)
         IF (IEVT.EQ.1) GOTO 100
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
         DO J = 1, BIN
            TRY = HMIN + HSTP*J
            NSOL = 0
            CALL DILEPSOL(PBT,PLP,PBX,PLM,PXMISS,PYMISS,TRY,MW,MB)
            CALL DILEPSOL(PBT,PLM,PBX,PLP,PXMISS,PYMISS,TRY,MW,MB)
            H(J) = H(J) + NSOL

         ENDDO
      ENDDO
 100  CONTINUE
      CALL LHEIN (1)
      DO J = 1, BIN
         WRITE (11,*) HMIN + HSTP*J, H(J)
      ENDDO
      STOP
      END
