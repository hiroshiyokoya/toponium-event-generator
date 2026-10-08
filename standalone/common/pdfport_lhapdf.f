*     -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6
*     TB, BB, CB, SB, UB, DB, G, D, U, S, C, B, T
      SUBROUTINE SETPDF (IPDF,IPORD,ISUB,IERR)
      IMPLICIT NONE
      INTEGER IPDF,IPORD,ISUB,IERR
      INTEGER JPDF,JPORD,JSUB,JERR
      COMMON/PDF/JPDF,JPORD,JSUB,JERR
      EXTERNAL INITPDFSET,INITPDF
      JPDF = IPDF
      JPORD = IPORD
      JSUB = ISUB
      JERR = IERR
      WRITE (6,*) " SETPDF: LHAPDF"
      IF (IPDF.EQ.4) THEN
         IF (IPORD.EQ.1) THEN
*            CALL INITPDFSETBYNAME ('cteq4l.LHgrid')
            CALL INITPDFSETBYNAME ('cteq6l.LHpdf') ! NLO as
*            CALL INITPDFSETBYNAME ('cteq6ll.LHpdf') ! LO as
            CALL INITPDF (0)
         ELSEIF (IPORD.EQ.2) THEN
            CALL INITPDFSETBYNAME ('cteq6m.LHpdf')
            CALL INITPDF (0)
         ELSE
            WRITE (6,*) "no IPORD",IPORD
            STOP
         ENDIF
      ELSEIF (IPDF.EQ.3) THEN
         IF (IPORD.EQ.2) THEN
            CALL INITPDFSETBYNAME ('MRST2004nlo.LHgrid')
            CALL INITPDF (0)
         ELSEIF (IPORD.EQ.3) THEN
            CALL INITPDFSETBYNAME ('MRST2004nnlo.LHgrid')
            CALL INITPDF (0)
         ELSE
            WRITE (6,*) "no IPORD",IPORD
            STOP
         ENDIF
      ELSEIF (IPDF.EQ.2) THEN
         IF (IPORD.EQ.1) THEN
            CALL INITPDFSETBYNAME ('a02m_lo.LHgrid')
            CALL INITPDF (0)
         ELSEIF (IPORD.EQ.2) THEN
            CALL INITPDFSETBYNAME ('a02m_nlo.LHgrid')
            CALL INITPDF (0)
         ELSE
            WRITE (6,*) "no IPORD",IPORD
            STOP
         ENDIF
      ELSEIF (IPDF.EQ.1) THEN
         IF (IPORD.EQ.1) THEN
            CALL INITPDFSETBYNAME ('GRV98lo.LHgrid')
            CALL INITPDF (0)
         ELSEIF (IPORD.EQ.2) THEN
            CALL INITPDFSETBYNAME ('GRV98nlo.LHgrid')
            CALL INITPDF (0)
         ELSE
            WRITE (6,*) "no IPORD",IPORD
            STOP
         ENDIF
      ELSEIF (IPDF.EQ.0) THEN   ! AS setting for RENO
         IF (IPORD.EQ.1) THEN
            CALL INITPDFSETBYNAME ('GRV98lo.LHgrid')
            CALL INITPDF (0)
         ELSEIF (IPORD.EQ.2) THEN
            CALL INITPDFSETBYNAME ('GRV98nlo.LHgrid')
            CALL INITPDF (0)
         ELSE
            WRITE (6,*) "no IPORD",IPORD
            STOP
         ENDIF
      ELSE
         WRITE (6,*) "no IPDF",IPDF
         STOP
      ENDIF
      RETURN
      END
C
      SUBROUTINE CALLPDF (Q,X,PDF) ! Return x*f(x)
      IMPLICIT NONE
      DOUBLE PRECISION Q,X,PDF(-6:6)
      INTEGER I
      EXTERNAL EVOLVEPDF
      CALL EVOLVEPDF (X,Q,PDF)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION ASQCD (Q,INF) ! return alpha_s(Q)
      IMPLICIT NONE
      DOUBLE PRECISION Q
      INTEGER INF
      DOUBLE PRECISION ALPHASPDF
      EXTERNAL ALPHASPDF
      ASQCD = ALPHASPDF (Q)
      RETURN
      END
C
