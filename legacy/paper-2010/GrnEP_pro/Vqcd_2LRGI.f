C...  2-Loop RG-improved (1-loop) QCD potential
C...  Described in Y.Sumino et.al., PRD47,56 (1993).
C...  It call POTQCD.f and EI.f, which are provided by JLCCVS.
C...  Only color-singlet potential is implemented.
      DOUBLE PRECISION FUNCTION VQCD (R)
      IMPLICIT NONE
      DOUBLE PRECISION R,V
      EXTERNAL V
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      IF (ICLR.NE.1) THEN
         WRITE (6,*) "No 2LRGI for Color-Octet", ICLR
         STOP
      ENDIF
      VQCD = V(R)
      RETURN
      END
