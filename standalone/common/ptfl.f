C     PTFL : RETURN THE PRODUCT OF TWO PDFS; F(X1) * F(X2)
C     IPRC : 1:P-P, 2:P-PBAR
C     PDF MUST BE DETERMINED PREVIOUSELY (pdfport)
C
C     This version is linked with the older programs main_* and tot_*
C     in bwblv/ and blvblv/ (sigh_bwblv.f, sigh_blvblv.f). The EG_*
C     generators use the PTFL inside their sigh_*E.f, which keeps the
C     two beam orientations of q qbar apart and omits b quarks.
C
C     Input : X1, X2 momentum fractions of beams 1 and 2; MU2 =
C             factorisation scale squared [GeV^2]; IPRC 1 pp, 2 p pbar
C     Output: products of parton densities f_i(X1) f_j(X2)
C       GG : g g
C       QG : (sum of quarks and antiquarks d..b)(X1) * g(X2)
C       GQ : g(X1) * (sum of quarks and antiquarks d..b)(X2)
C       QQ : sum over q = d,u,s,c,b of q qbar + qbar q (both beam
C            orientations); for p pbar the antiproton density of qbar
C            is the proton density of q
C     CALLPDF returns x f(x), hence the division by X1 X2.
      SUBROUTINE PTFL (X1,X2,QQ,GG,QG,GQ,MU2,IPRC)
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2,QQ,GG,QG,GQ,MU2
      DOUBLE PRECISION MU
      INTEGER IPRC
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
