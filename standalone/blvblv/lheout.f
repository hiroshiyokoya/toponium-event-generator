      SUBROUTINE LHEOUT (ILHE)
      IMPLICIT NONE
      INTEGER I,J,K
      INTEGER ILHE,IO
      COMMON/FILE/IO
      DATA IO/10/
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      INTEGER ISPRC
      COMMON/PROC/ISPRC
      DOUBLE PRECISION RAN
      DOUBLE PRECISION SIG,ERR
      COMMON/BSS/SIG,ERR
      include 'run.inc'
      include 'coupl_sm.inc'
      include 'nexternal.inc'
      INTEGER IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP
      DOUBLE PRECISION EBMUP(2)
      DOUBLE PRECISION XSECUP(10),XERRUP(10),XMAXUP(10)
      INTEGER LPRUP(10)
      INTEGER NUP,IDPRUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP
      INTEGER IDUP(NEXTERNAL),ISTUP(NEXTERNAL),
     -     MOTHUP(2,NEXTERNAL),ICOLUP(2,NEXTERNAL)
      DOUBLE PRECISION MPAR(NEXTERNAL)
      DATA MPAR/0D0,0D0,5D0,0D0,0D0,5D0,0D0,0D0/
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION P(0:3,NEXTERNAL),PH(0:3,NEXTERNAL)
      COMMON/MOM/P
      INTEGER HHEL(NEXTERNAL)
      COMMON/HEL2/HHEL
      DOUBLE PRECISION RTC(2)
      COMMON/COL/RTC
      IF (ILHE.EQ.-1) THEN
         OPEN(IO,FILE="events.lhe")
         WRITE (IO,*) "<LesHouchesEvents version=""1.0"">"
         WRITE (IO,*) "<init>"
         IF (IPRC.EQ.1) THEN
            IDBMUP(1) = 2212
            IDBMUP(2) = 2212
            PDFGUP(1) = 0
            PDFGUP(2) = 0
            PDFSUP(1) = 10042
            PDFSUP(2) = 10042
         ELSEIF (IPRC.EQ.2) THEN
            IDBMUP(1) = 2212
            IDBMUP(2) = -2212
            PDFGUP(1) = 0
            PDFGUP(2) = 0
            PDFSUP(1) = 10042
            PDFSUP(2) = 10042
         ENDIF
         EBMUP(1) = EBEAM(1)
         EBMUP(2) = EBEAM(2)
         IDWTUP = 3
         NPRUP  = 1
         WRITE (IO,1) IDBMUP,EBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP
 1       FORMAT (2I6,2X,2E16.8,1X,2I2,2I7,2I2)
         DO I = 1,NPRUP
            XSECUP(I) = SIG
            XERRUP(I) = ERR
            XMAXUP(I) = 1D0
            LPRUP (I) = 1
            WRITE (IO,2) XSECUP(I),XERRUP(I),XMAXUP(I),LPRUP(I)
 2          FORMAT (3(E16.8,1X),I2)
         ENDDO
         WRITE (IO,*) "</init>"
      ELSEIF (ILHE.EQ.0) THEN
         WRITE (IO,*) "<event>"
         NUP = NEXTERNAL
         IDPRUP = 1
         XWGTUP = 1D0
         SCALUP = SCALE
         AQEDUP = ALPHA
         AQCDUP = ALFAS
         WRITE (IO,3) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
 3       FORMAT (2(I2,1X),4(E16.8,2X))
         CALL PCM2HCM (X1,X2,P,PH)
C.....
         IDUP(3) =  5
         IDUP(4) =-13
         IDUP(5) = 14
         IDUP(6) = -5
         IDUP(7) = 13
         IDUP(8) =-14
C.....
         ISTUP(1) = -1
         ISTUP(2) = -1
         ISTUP(3) =  1
         ISTUP(4) =  1
         ISTUP(5) =  1
         ISTUP(6) =  1
         ISTUP(7) =  1
         ISTUP(8) =  1
C.....
         MOTHUP(1,1) = 0
         MOTHUP(2,1) = 0
         MOTHUP(1,2) = 0
         MOTHUP(2,2) = 0
         MOTHUP(1,3) = 1
         MOTHUP(2,3) = 2
         MOTHUP(1,4) = 1
         MOTHUP(2,4) = 2
         MOTHUP(1,5) = 1
         MOTHUP(2,5) = 2
         MOTHUP(1,6) = 1
         MOTHUP(2,6) = 2
         MOTHUP(1,7) = 1
         MOTHUP(2,7) = 2
         MOTHUP(1,8) = 1
         MOTHUP(2,8) = 2
C.....
         IF (ISPRC.EQ.1) THEN   ! QQB
            IDUP(1) =  1
            IDUP(2) = -1
            ICOLUP(1,1) = 501
            ICOLUP(2,1) = 0
            ICOLUP(1,2) = 0
            ICOLUP(2,2) = 503
            ICOLUP(1,3) = 501
            ICOLUP(2,3) = 0
            ICOLUP(1,4) = 0
            ICOLUP(2,4) = 0
            ICOLUP(1,5) = 0
            ICOLUP(2,5) = 0
            ICOLUP(1,6) = 0
            ICOLUP(2,6) = 503
            ICOLUP(1,7) = 0
            ICOLUP(2,7) = 0
            ICOLUP(1,8) = 0
            ICOLUP(2,8) = 0
         ELSEIF (ISPRC.EQ.2) THEN ! QBQ
            IDUP(1) = -1
            IDUP(2) =  1
            ICOLUP(1,1) = 0
            ICOLUP(2,1) = 503
            ICOLUP(1,2) = 501
            ICOLUP(2,2) = 0
            ICOLUP(1,3) = 501
            ICOLUP(2,3) = 0
            ICOLUP(1,4) = 0
            ICOLUP(2,4) = 0
            ICOLUP(1,5) = 0
            ICOLUP(2,5) = 0
            ICOLUP(1,6) = 0
            ICOLUP(2,6) = 503
            ICOLUP(1,7) = 0
            ICOLUP(2,7) = 0
            ICOLUP(1,8) = 0
            ICOLUP(2,8) = 0
         ELSEIF (ISPRC.EQ.3) THEN
            IDUP(1) = 21
            IDUP(2) = 21
            ICOLUP(1,1) = 501
            ICOLUP(2,1) = 502
            ICOLUP(1,2) = 502
            ICOLUP(2,2) = 501
            ICOLUP(1,3) = 503
            ICOLUP(2,3) = 0
            ICOLUP(1,4) = 0
            ICOLUP(2,4) = 0
            ICOLUP(1,5) = 0
            ICOLUP(2,5) = 0
            ICOLUP(1,6) = 0
            ICOLUP(2,6) = 503
            ICOLUP(1,7) = 0
            ICOLUP(2,7) = 0
            ICOLUP(1,8) = 0
            ICOLUP(2,8) = 0
         ELSEIF (ISPRC.EQ.4) THEN
            IDUP(1) = 21
            IDUP(2) = 21
            RAN = RAND(0)
            IF(RAN.LE.RTC(1)) THEN
               ICOLUP(1,1) = 501
               ICOLUP(2,1) = 502
               ICOLUP(1,2) = 502
               ICOLUP(2,2) = 503
            ELSEIF (RAN.GT.RTC(1)) THEN
               ICOLUP(1,1) = 502
               ICOLUP(2,1) = 503
               ICOLUP(1,2) = 501
               ICOLUP(2,2) = 502
            ENDIF
            ICOLUP(1,3) = 501
            ICOLUP(2,3) = 0
            ICOLUP(1,4) = 0
            ICOLUP(2,4) = 0
            ICOLUP(1,5) = 0
            ICOLUP(2,5) = 0
            ICOLUP(1,6) = 0
            ICOLUP(2,6) = 503
            ICOLUP(1,7) = 0
            ICOLUP(2,7) = 0
            ICOLUP(1,8) = 0
            ICOLUP(2,8) = 0
         ENDIF
         DO J = 1,NEXTERNAL
            WRITE (IO,4) IDUP(J),ISTUP(J),MOTHUP(1,J),MOTHUP(2,J),
     -           ICOLUP(1,J),ICOLUP(2,J),(PH(K,J),K=1,3),PH(0,J),
     -           MPAR(J),0.,HHEL(J)
 4          FORMAT (2(I8,1X),4I5,5E16.8,F4.0,1X,I2)
         ENDDO
         WRITE (IO,*) "</event>"
      ELSEIF (ILHE.EQ.1) THEN
         WRITE (IO,*) "</LesHouchesEvents>"
         CLOSE(IO)
      ENDIF
      RETURN
      END
