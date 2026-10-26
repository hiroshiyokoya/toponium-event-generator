C=======================================================================
C     GGTTRES: resonant g g -> t tbar amplitudes (t-, u- and s-channel
C     diagrams) with the bound-state correction (arXiv:1007.0075,
C     Sec. 2.1; docs/green-function-tables.md, Sec. 3). Called from
C     the g g matrix elements me_gg*.f.
C
C     Input (HELAS wave functions, 6 complex components each)
C       G1, G2 : incoming gluons
C       TP, TB : off-shell top and antitop, built from the decay
C                products (t -> b W+, tbar -> bbar W-)
C       ICLR   : colour state, 1 singlet, 2 octet (selects the table)
C     Output
C       AMPT, AMPU, AMPS : t-, u- and s-channel amplitudes, modified
C                          according to IGRN, ITR and INR (thr.inc)
C     COMMON/BW/ MT1S, MT2S, S : s_t, s_tbar and s-hat [GeV^2], set by
C                the matrix-element routine before the call
C
C     Switches (thr.inc; docs/versions-2010-2015.md):
C       IGRN=0 : amplitudes unchanged (free top propagators)
C       IGRN=1 : overall correction prescription (OCP). All three
C                amplitudes A are replaced by
C                A * ( COR*Gamma + (1-COR)*FNR(INR) )
C       IGRN=2 : total-angular-momentum prescription (TAP). In the
C                t- and u-channel amplitudes only the part computed
C                with the modified top propagator FVIXXZ (fvixxz.F) is
C                multiplied by Gamma; the s-channel amplitude is not
C                changed. The paper uses IGRN=1.
C       ITR=0  : COR = 1
C       ITR=1  : COR = CORR(s_t, s_tbar, ...) (corr.f)
C       INR    : FNR = 0, 1, 1 for INR = 0, 1, 2; the remainder
C                (1-COR)*A is kept only for INR >= 1. For ITR=0 it
C                vanishes anyway.
C=======================================================================
      SUBROUTINE GGTTRES (G1,G2,TP,TB,AMPT,AMPU,AMPS,ICLR)
      IMPLICIT NONE
      INTEGER ICLR
      DOUBLE COMPLEX G1(6),G2(6),TP(6),TB(6)
      DOUBLE COMPLEX AMPT,AMPU,AMPS,AMPT0,AMPU0
      DOUBLE COMPLEX W(6,3)
      include "coupl_sm.inc"
      INCLUDE "thr.inc"
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      DOUBLE COMPLEX GAMMA
      DOUBLE PRECISION MT1S,MT2S,S
      COMMON/BW/MT1S,MT2S,S
      DOUBLE PRECISION COR,CORR
      EXTERNAL CORR
C     FNR(INR): weight of the uncorrected remainder (1-COR)*A
      DOUBLE PRECISION FNR(0:2)
      DATA FNR/0.,1.,1./
C.....
C     COR: 1 for ITR=0; for ITR=1 the factor of CORR, which takes the
C     squared top virtualities MT1S, MT2S [GeV^2]
      IF (ITR.EQ.0) THEN
         COR = 1D0
      ELSEIF (ITR.EQ.1) THEN
         COR = CORR (MT1S,MT2S,TMASS,TWIDTH,gfermi,WMASS,BMASS)
      ENDIF
C.....
C     Amplitudes with the ordinary (Breit-Wigner) top propagator
      CALL FVIXXX (TB,G2,GG,TMASS,TWIDTH,W(1,1))
      CALL IOVXXX (W(1,1),TP,G1,GG,AMPT) ! t-channel
      CALL FVIXXX (TB,G1,GG,TMASS,TWIDTH,W(1,2))
      CALL IOVXXX (W(1,2),TP,G2,GG,AMPU) ! u-channel
      CALL JVVXXX (G1,G2,G ,ZERO ,ZERO  ,W(1,3))
      CALL IOVXXX (TB,TP,W(1,3),GG,AMPS) ! s-channel
C     IGRN=1 (OCP): the same factor on all three amplitudes
      IF ( IGRN.EQ.1 ) THEN     ! OCP
         CALL GRNMOD (GAMMA,MT1S,MT2S,S,TMASS,ICLR,IGRN)
         AMPT = AMPT * ( COR*GAMMA + (1D0-COR) * FNR(INR) )
         AMPU = AMPU * ( COR*GAMMA + (1D0-COR) * FNR(INR) )
         AMPS = AMPS * ( COR*GAMMA + (1D0-COR) * FNR(INR) )
C     IGRN=2 (TAP): AMPT0, AMPU0 are the t- and u-channel amplitudes
C     computed with FVIXXZ; A -> A0*Gamma + (A - A0) for A = AMPT,
C     AMPU, times COR, plus the remainder. AMPS is left unchanged.
      ELSEIF ( IGRN.EQ.2 ) THEN ! TAP
         CALL FVIXXZ (TB,G2,GG,TMASS,TWIDTH,W(1,1))
         CALL IOVXXX (W(1,1),TP,G1,GG,AMPT0) ! t-channel
         CALL FVIXXZ (TB,G1,GG,TMASS,TWIDTH,W(1,2))
         CALL IOVXXX (W(1,2),TP,G2,GG,AMPU0) ! u-channel
         CALL GRNMOD (GAMMA,MT1S,MT2S,S,TMASS,ICLR,IGRN)
         AMPT = ( AMPT0*GAMMA + AMPT - AMPT0 ) * COR +
     -        AMPT*(1D0-COR) * FNR(INR)
         AMPU = ( AMPU0*GAMMA + AMPU - AMPU0 ) * COR +
     -        AMPU*(1D0-COR) * FNR(INR)
      ENDIF
      RETURN
      END
