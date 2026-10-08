C     PTFL : RETURN THE PRODUCT OF TWO PDFS; F(X1) * F(X2)
C     IPRC : 1:P-P, 2:P-PBAR
C     PDF MUST BE DETERMINED PREVIOUSELY (pdfport)
C
      SUBROUTINE PTFL (X1,X2,QQ,GG,QG,GQ,MU2,IPRC)
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2,QQ,GG,QG,GQ,MU2
      DOUBLE PRECISION MU
      INTEGER IPRC,IORD,IEXP
      DOUBLE PRECISION PDF1(-6:6), PDF2(-6:6)
      EXTERNAL CALLPDF
      INTEGER IPDF,IPORD,ISUB,IERR
      COMMON/PDF/IPDF,IPORD,ISUB,IERR
      MU = DSQRT (MU2)
      CALL CALLPDF (MU,X1,PDF1)
      CALL CALLPDF (MU,X2,PDF2)
      GG = PDF1(0) * PDF2(0) / X1 / X2
      QG = (PDF1(1) + PDF1(2) + PDF1(3) + PDF1(4) + PDF1(5)
     -     + PDF1(-1) + PDF1(-2) + PDF1(-3) + PDF1(-4) + PDF1(-5))
     -     * PDF2(0) / X1 / X2
      GQ = PDF1(0) * (PDF2(1) + PDF2(2) + PDF2(3) + PDF2(4) + PDF2(5)
     -     + PDF2(-1) + PDF2(-2) + PDF2(-3) + PDF2(-4) + PDF2(-5))
     -     / X1 / X2
      IF(IPRC.EQ.1) THEN        ! P-P
         QQ = ( PDF1( 1)*PDF2(-1) + PDF1( 2)*PDF2(-2) + PDF1(3)*PDF2(-3)
     -        + PDF1( 4)*PDF2(-4) + PDF1( 5)*PDF2(-5)
     -        + PDF1(-1)*PDF2( 1) + PDF1(-2)*PDF2( 2) + PDF1(-3)*PDF2(3)
     -        + PDF1(-4)*PDF2( 4) + PDF1(-5)*PDF2( 5) ) / X1 / X2
      ELSEIF (IPRC.EQ.2) THEN   ! P-PBAR
         QQ = (PDF1(1)*PDF2(1) + PDF1(2)*PDF2(2) + PDF1(3)*PDF2(3)
     -        + PDF1(4)*PDF2(4) + PDF1(5)*PDF2(5)
     -        + PDF1(-1)*PDF2(-1) + PDF1(-2)*PDF2(-2)
     -        + PDF1(-3)*PDF2(-3) + PDF1(-4)*PDF2(-4)
     -        + PDF1(-5)*PDF2(-5)) / X1 / X2
      ENDIF
      RETURN
      END
