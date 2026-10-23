C=======================================================================
C     Interface to LHAPDF 6 through its Fortran (LHAGlue) routines
C     INITPDFSETBYNAME, INITPDF, EVOLVEPDF and ALPHASPDF.
C       SETPDF  : choose a PDF set by (IPDF, IPORD) and store the
C                 arguments in COMMON/PDF/. Called by main_bwblv.f,
C                 tot_bwblv.f and main_blvblv.f. The EG_*.f generators
C                 call INITPDFSETBYNAME('cteq6ll.LHpdf') directly.
C       CALLPDF : PDF(-6:6) = x f(x,Q) for all flavours
C       ASQCD   : alpha_s(Q) of the PDF set
C     Flavour index of PDF(-6:6) (the two lines below):
C=======================================================================
*     -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6
*     TB, BB, CB, SB, UB, DB, G, D, U, S, C, B, T
C     SETPDF: IPDF = 1 GRV98, 2 Alekhin02, 3 MRST2004, 4 CTEQ6, 0 GRV98
C     (as 1); IPORD = 1 LO, 2 NLO, 3 NNLO (not every pair exists).
C     Member 0 of the set. ISUB and IERR are only stored. Stops for
C     an unknown choice. The set names are LHAPDF 5 style names.
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
C     CALLPDF: Q [GeV], X; returns PDF(i) = x f_i(x,Q), i = -6..6
      SUBROUTINE CALLPDF (Q,X,PDF) ! Return x*f(x)
      IMPLICIT NONE
      DOUBLE PRECISION Q,X,PDF(-6:6)
      EXTERNAL EVOLVEPDF
      CALL EVOLVEPDF (X,Q,PDF)
      RETURN
      END
C
C     ASQCD: alpha_s(Q), Q [GeV], from the LHAPDF set; INF is unused
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
