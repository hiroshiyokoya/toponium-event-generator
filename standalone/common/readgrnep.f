      SUBROUTINE READHEAD (ICLR)
      IMPLICIT NONE
      INTEGER IOF,ICLR
      PARAMETER (IOF=13)
      DOUBLE PRECISION MT,GAMT,ASB
      COMMON/GrnPAR/MT,GAMT,ASB
      INCLUDE 'filedata.inc'
      IF (ICLR.EQ.1) THEN
         OPEN (IOF,FILE=T1)
      ELSEIF (ICLR.EQ.2) THEN
         OPEN (IOF,FILE=T8)
      ELSE
         WRITE (6,*) "ERROR in READHEAD", ICLR
         STOP
      ENDIF
      READ (IOF,*) MT, GAMT, ASB
      CLOSE(IOF)
      RETURN
      END
C
      SUBROUTINE READGRN (GAMMA,E,P,ICLR)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER I,J,K,ICLR
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      DOUBLE PRECISION E,P
      INTEGER IE,IP,EINT,PINT
      DOUBLE PRECISION EE,PP
      EXTERNAL EE,PP,EINT,PINT
      DOUBLE COMPLEX GRN,GAMMA
      INTEGER ITBL(2)
      COMMON/GrnCHK/ITBL
      IF (ITBL(ICLR).NE.33) THEN
         CALL READTABLE (ICLR)
         ITBL (ICLR) = 33
      ENDIF
      CALL HOKAN2D (E,P,GRN,ICLR)
      GAMMA = GRN
      RETURN
      END
C
      SUBROUTINE READTABLE (ICLR)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IOF,ICLR
      PARAMETER (IOF=14)
      INTEGER ME,MP,IE,IP
      DOUBLE PRECISION EMN,EMX
      DOUBLE PRECISION PMN,PMX
      DOUBLE PRECISION R
      DOUBLE PRECISION RE,IM
      DOUBLE COMPLEX GRNEPT(0:NE1,0:NP1,2),GRNEPH(0:NE2,0:NP2,2)
      COMMON/GrnTBL/GRNEPT,GRNEPH
      DOUBLE PRECISION MT,GAMT,ASB
      COMMON/GrnPAR/MT,GAMT,ASB
      INCLUDE 'filedata.inc'
      IF (ICLR.EQ.1) THEN
         WRITE (6,*) "   Open file : ", T1
         OPEN (IOF,FILE=T1)
      ELSEIF (ICLR.EQ.2) THEN
         WRITE (6,*) "   Open file : ", T8
         OPEN (IOF,FILE=T8)
      ELSE
         STOP
      ENDIF
      READ (IOF,*) MT,GAMT,ASB
      WRITE (6,*) "MT, GAMT, ASB =",MT,GAMT,ASB
      READ (IOF,*) ME, MP
      READ (IOF,*) EMN, EMX
      READ (IOF,*) PMN, PMX
      IF ( ME.EQ.NE1 .AND. MP.EQ.NP1 ) THEN
         DO IE = 0,NE1,1
            DO IP = 0,NP1,1
               READ (IOF,*) RE,IM
               GRNEPT(IE,IP,ICLR) = DCMPLX(RE,IM)
            ENDDO
            READ(IOF,*)
         ENDDO
      ELSE
         WRITE (6,*) "READGRN (thre): Incompatible Table Size"
         WRITE (6,*) "Table : ", ME, MP
         WRITE (6,*) "Inc :   ", NE1, NP1
         STOP
      ENDIF
      CLOSE (IOF)
C
      IF (ICLR.EQ.1) THEN
         WRITE (6,*) "   Open file : ", H1
         OPEN (IOF,FILE=H1)
      ELSEIF (ICLR.EQ.2) THEN
         WRITE (6,*) "   Open file : ", H8
         OPEN (IOF,FILE=H8)
      ELSE
         STOP
      ENDIF
      READ (IOF,*) MT,GAMT,ASB
      WRITE (6,*) "MT, GAMT, ASB =",MT,GAMT,ASB
      READ (IOF,*) ME, MP
      READ (IOF,*) EMN, EMX
      READ (IOF,*) R
      IF ( ME.EQ.NE2 .AND. MP.EQ.NP2 .AND. R.EQ.RR ) THEN
         DO IE = 0,NE2,1
            DO IP = 0,NP2,1
               READ (IOF,*) RE,IM
               GRNEPH(IE,IP,ICLR) = DCMPLX(RE,IM)
            ENDDO
            READ(IOF,*)
         ENDDO
      ELSE
         WRITE (6,*) "READGRN (high): Incompatible Table Size"
         WRITE (6,*) "Table :", ME, MP, R
         WRITE (6,*) "Inc :  ", NE2, NP2, RR
         STOP
      ENDIF
      CLOSE (IOF)
      WRITE (6,*) ICLR, ": READTABLE"
      RETURN
      END
C
      SUBROUTINE READGRNTBL (GREEN,IE,IP,ICLR,IT)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IE,IP,ICLR
      DOUBLE COMPLEX GREEN
      DOUBLE COMPLEX GRNEPT(0:NE1,0:NP1,2), GRNEPH(0:NE2,0:NP2,2)
      COMMON/GrnTBL/GRNEPT, GRNEPH
      IF (IT.EQ.1) THEN
         GREEN = GRNEPT (IE,IP,ICLR)
      ELSEIF (IT.EQ.2) THEN
         GREEN = GRNEPH (IE,IP,ICLR)
      ENDIF
      RETURN
      END
C
      SUBROUTINE HOKAN2D (EI,PI,GRN,ICLR)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION E,P,EI,PI
      INTEGER ICLR
      INTEGER IE,IP,IE0,IE1,IP0,IP1
      DOUBLE COMPLEX GRN
      INTEGER EINT,PINT
      DOUBLE PRECISION EE,PP
      EXTERNAL EINT,PINT,EE,PP
      DOUBLE PRECISION E0,E1,P0,P1
      DOUBLE COMPLEX G0,G1,G00,G01,G10,G11
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      DOUBLE PRECISION MT,GAMT,ASB
      COMMON/GrnPAR/MT,GAMT,ASB
      GRN = DCMPLX(1D0,0D0)
      IF (EI.LT.EMAX1) THEN
         IT = 1
         IF (EI.LT.EMIN1) THEN
            IE = 0
            E = EI
         ELSEIF (EI.GE.EMIN1) THEN
            E = EI
            IE = EINT(E,1)
         ELSE
            WRITE (6,*) "ERROR; IT = 1"
            STOP
         ENDIF
         PMIN = PMIN1
         PMAX = PMAX1
      ELSEIF (EI.GE.EMIN2) THEN
         IT = 2
         IF (EI.GT.EMAX2) THEN
C...  extrapolation
            IE = NE2-1
            E = EI
C...  stop at the maximum of grid
*            E = EMAX2
*            IE = NE2
         ELSEIF (EI.LE.EMAX2) THEN
            E = EI
            IE = EINT(E,2)
         ELSE
            WRITE (6,*) "ERROR; IT = 2"
            STOP
         ENDIF
         POS = DSQRT(MT*E)
         PMIN = PP (0,IT)
         PMAX = PP (NP2,IT)
      ELSE
         WRITE (6,*) "IT"
         WRITE(6,*) EI,EMIN1,EMAX1,EMIN2,EMAX2
         STOP
      ENDIF
C
      IF (PI.GT.PMAX) THEN
C...  extrapolation to outside grid
*         P = PI
*         IP = PINT(PMAX,IT) - 1
C...  stop at the maxmum
         P = PMAX
         IP = PINT(PMAX,IT)
      ELSEIF (PI.LT.PMIN) THEN
         P = PI
         IP = 0
      ELSE
         P = PI
         IP = PINT(P,IT)
      ENDIF
C
      IF (ABS(EE(IE,IT)-E).LE.1D-4 .AND. ABS(PP(IP,IT)-P).LE.1D-4) THEN
         CALL READGRNTBL (GRN,IE,IP,ICLR,IT)
      ELSEIF ( ABS(PP(IP,IT)-P) .LE. 1D-4 ) THEN
         IE0 = IE
         IE1 = IE0 + 1
         E0 = EE(IE0,IT)
         E1 = EE(IE1,IT)
         CALL READGRNTBL (G0,IE0,IP,ICLR,IT)
         CALL READGRNTBL (G1,IE1,IP,ICLR,IT)
         GRN = (G1-G0)/(E1-E0)*(E-E0) + G0
      ELSEIF ( ABS(EE(IE,IT)-E) .LE. 1D-4 ) THEN
         IP0 = IP
         IP1 = IP0 + 1
         P0 = PP(IP0,IT)
         P1 = PP(IP1,IT)
         CALL READGRNTBL (G0,IE,IP0,ICLR,IT)
         CALL READGRNTBL (G1,IE,IP1,ICLR,IT)
         GRN = (G1-G0)/(P1-P0)*(P-P0) + G0
      ELSE
         IE0 = IE
         IE1 = IE0 + 1
         E0 = EE(IE0,IT)
         E1 = EE(IE1,IT)
         IP0 = IP
         IP1 = IP0 + 1
         P0 = PP(IP0,IT)
         P1 = PP(IP1,IT)
         CALL READGRNTBL (G00,IE0,IP0,ICLR,IT)
         CALL READGRNTBL (G01,IE0,IP1,ICLR,IT)
         G0 = (G01-G00)/(P1-P0)*(P-P0) + G00
         CALL READGRNTBL (G10,IE1,IP0,ICLR,IT)
         CALL READGRNTBL (G11,IE1,IP1,ICLR,IT)
         G1 = (G11-G10)/(P1-P0)*(P-P0) + G10
         GRN = (G1-G0)/(E1-E0)*(E-E0) + G0
      ENDIF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION EE (IE,IT)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IE
      IF (IT.EQ.1) THEN
         NE = NE1
         EMIN = EMIN1
         EMAX = EMAX1
      ELSEIF (IT.EQ.2) THEN
         NE = NE2
         EMIN = EMIN2
         EMAX = EMAX2
      ENDIF
      EE = EMIN + DBLE(IE)/DBLE(NE) * (EMAX - EMIN)
      RETURN
      END
C
      INTEGER FUNCTION EINT (E,IT)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION E
      IF (IT.EQ.1) THEN
         NE = NE1
         EMIN = EMIN1
         EMAX = EMAX1
      ELSEIF (IT.EQ.2) THEN
         NE = NE2
         EMIN = EMIN2
         EMAX = EMAX2
      ENDIF
      EINT = INT(NE*(E-EMIN)/(EMAX-EMIN))
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PP (IP,IT)
      IMPLICIT NONE
      INTEGER IP,M
      DOUBLE PRECISION A,R1
      INCLUDE 'tblgrnep.inc'
      IF (IT.EQ.1) THEN
         PMIN = PMIN1
         PMAX = PMAX1
         PP = PMIN + DBLE(IP)/DBLE(NP1) * (PMAX - PMIN)
      ELSEIF (IT.EQ.2) THEN
         M = (NP2-1)/2
         A = (1D0-RR)/(1D0-RR**M) * POS
         R1 = 1D0/RR
         IF (IP.LE.M) THEN
            PP = POS * (1D0-RR**(IP+1)) / (1D0-RR**(M+1))
         ELSE
            PP = POS + A*RR**M*(1D0-R1**(IP-M+1))/(1D0-R1)
         ENDIF
      ENDIF
      RETURN
      END
C
      INTEGER FUNCTION PINT (P,IT)
      IMPLICIT NONE
      DOUBLE PRECISION P,A,R1
      INTEGER M
      INCLUDE 'tblgrnep.inc'
      IF (IT.EQ.1) THEN
         PMIN = PMIN1
         PMAX = PMAX1
         PINT = INT(NP1*(P-PMIN)/(PMAX-PMIN))
      ELSEIF (IT.EQ.2) THEN
         M = (NP2-1)/2
         A = (1D0-RR)/(1D0-RR**M) * POS
         R1 = 1D0/RR
         IF (P.LE.POS) THEN
            PINT = INT(DLOG(1D0-P/POS*(1D0-RR**(M+1)))/DLOG(RR)) - 1
         ELSE
            PINT = INT(DLOG(1D0-(P-POS)/(A*RR**M)*(1D0-R1))/DLOG(R1))
     -           + M - 1
         ENDIF
      ENDIF
      RETURN
      END
