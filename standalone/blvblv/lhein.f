      SUBROUTINE LHEIN (ILHE)
      IMPLICIT NONE
      INTEGER I,J
      INTEGER ILHE,IO,IBEG
      COMMON/FILE/IO
      DATA IO/99/
      INTEGER ICNT
      COMMON/EVT/ICNT
      DATA ICNT/0/
      include 'nexternal.inc'
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(NEXTERNAL),ISTUP(NEXTERNAL),MOTHUP(2,NEXTERNAL),
     &     ICOLUP(2,NEXTERNAL),PUP(5,NEXTERNAL),VTIMUP(NEXTERNAL),
     &     SPINUP(NEXTERNAL)
      SAVE /HEPEUP/
C...  Lines to read in assumed never longer than 200 characters.
      INTEGER MAXLEN
      PARAMETER (MAXLEN=100)
      CHARACTER*(MAXLEN) STRING
C...  Format for reading lines.
      CHARACTER*6 STRFMT
      COMMON/STR/STRFMT
      IF (ILHE.EQ.-1) THEN
         STRFMT='(A000)'
         WRITE (STRFMT(3:5),'(I3)') MAXLEN
         OPEN (IO,FILE="events.lhe")
         WRITE (6,*) "Open",IO
      ELSEIF (ILHE.EQ.0) THEN
C...  Loop until finds line beginning with "<event>" or "<event ".
 100     READ (IO,STRFMT,END=130,ERR=130) STRING
         IBEG=0
 110     IBEG=IBEG+1
C...  Allow indentation.
         IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 110
         IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &        STRING(IBEG:IBEG+6).NE.'<event ') GOTO 100
C...  Read first line of event info.
         READ (IO,*,END=130,ERR=130) NUP,IDPRUP,XWGTUP,SCALUP,
     &        AQEDUP,AQCDUP
C...  Read NUP subsequent lines with information on each particle.
         DO 120 I=1,NUP
            READ (IO,*,END=130,ERR=130) IDUP(I),ISTUP(I),
     &           MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
     &           (PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I)
 120     CONTINUE
         ICNT=ICNT+1
         RETURN
C...  Error exit, typically when no more events.
 130     WRITE(6,*) ' Failed to read LHEF event information.'
         WRITE(6,*) ' Will assume end of file has been reached.'
         NUP=0
         ILHE=1
      ELSEIF (ILHE.EQ.1) THEN
         CLOSE(IO)
         WRITE (6,*) "Close",IO
         WRITE (6,*) "Num. of Event",ICNT
      ENDIF
      RETURN
      END
