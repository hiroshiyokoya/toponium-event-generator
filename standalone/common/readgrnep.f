C=======================================================================
C     Reading and interpolation of the Green-function tables
C     (docs/green-function-tables.md). The correction factor
C     Gamma(E,p) = G(E,p)/G0(E,p) is tabulated, for each colour state,
C     on two grids whose constants are in tblgrnep.inc:
C       threshold   grnep{1,8}thre.tbl : E = EMIN1..EMAX1 = -50..30
C                   GeV, NE1+1 = 801 points, uniform; p = PMIN1..PMAX1
C                   = 0..300 GeV, NP1+1 = 301 points, uniform
C       high energy grnep{1,8}high.tbl : E = EMIN2..EMAX2 = 30..830
C                   GeV, NE2+1 = 401 points, uniform; NP2+1 = 102
C                   p points, geometric around p0 = sqrt(m_t E) (PP)
C     (1 = colour singlet, 8 = octet). File names: filedata.inc; the
C     tables are read from the run directory.
C
C     Table format (list-directed):
C       line 1 : m_t, Gamma_t, alpha_s(mu_B)
C       line 2 : NE, NP (numbers of intervals)
C       line 3 : E_min, E_max [GeV]
C       line 4 : p_min, p_max [GeV] (thre) or the ratio RR (high)
C       then NE+1 blocks (E index outer) of NP+1 lines "Re Im"
C       (p index inner), each block followed by a separator line.
C
C     Entry point: READGRN(GAMMA,E,P,ICLR), called by GRNMOD.
C=======================================================================
C
C     READHEAD: reads only line 1 of the threshold table of colour
C     ICLR into COMMON/GrnPAR/. Not called by the stand-alone
C     generators (READTABLE reads the same values).
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
C     READGRN: Gamma(E,P) for colour ICLR (1 singlet, 2 octet).
C       E, P [GeV] as defined in GRNMOD; GAMMA complex, dimensionless.
C     On the first call for a colour both tables of that colour are
C     read (READTABLE); ITBL(ICLR) = 33 in COMMON/GrnCHK/ marks them
C     as loaded. m_t, Gamma_t and alpha_s(mu_B) of the table go to
C     COMMON/GrnPAR/ and are not compared with the generator's values
C     (docs/REVIEW.md, P3).
      SUBROUTINE READGRN (GAMMA,E,P,ICLR)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER   ICLR
      DOUBLE PRECISION E,P
      INTEGER  EINT,PINT
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
C     READTABLE: reads grnep{1,8}thre.tbl into GRNEPT(0:NE1,0:NP1,ICLR)
C     and grnep{1,8}high.tbl into GRNEPH(0:NE2,0:NP2,ICLR), both in
C     COMMON/GrnTBL/. Stops if the numbers of intervals (and, for the
C     high-energy table, RR) differ from tblgrnep.inc. The E and p
C     ranges of the header (EMN, EMX, PMN, PMX) are read but not
C     checked; the ranges in tblgrnep.inc are used (REVIEW.md, P14).
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
C     skip the separator line after each E block
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
C     High-energy table; its 4th header line is the ratio RR
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
C     READGRNTBL: table entry GREEN at grid indices (IE,IP) for colour
C     ICLR; IT = 1 threshold table, IT = 2 high-energy table (IT is
C     declared in tblgrnep.inc).
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
C     HOKAN2D: Gamma at (EI,PI) [GeV] by bilinear interpolation.
C     (PI is the momentum p here, not pi.) Result in GRN.
C
C     1. Table and E interval IE:
C        EI <  EMAX1 (30 GeV) : threshold table, IT = 1
C           EI < EMIN1 (-50 GeV): IE = 0, i.e. linear extrapolation
C           from the first interval (the 2010 code returns Gamma = 1
C           there; docs/REVIEW.md, P10)
C        EI >= EMIN2 (30 GeV) : high-energy table, IT = 2
C           EI > EMAX2 (830 GeV, sqrt(s-hat) > ~833 GeV): IE = NE2-1,
C           i.e. linear extrapolation from the last interval. Known
C           issue P1 (docs/REVIEW.md): the 2010 code returns Gamma = 1
C           there, and the extrapolated values are not physical (for
C           the singlet at E = 2000 GeV, Gamma = 0.62 - 0.25i). The
C           commented-out lines would freeze E at EMAX2 instead.
C        For IT = 2 the p grid moves with E: POS = p0 = sqrt(m_t E)
C        (COMMON/OS/, m_t from the table header) is set here and used
C        by PP and PINT.
C     2. p interval IP:
C        PI > PMAX : P is clamped to PMAX
C        PI < PMIN : IP = 0, linear extrapolation (only IT = 2, where
C                    PMIN = PP(0) > 0)
C     3. Interpolation. If E and/or P are on grid points (within
C        1D-4 GeV) the table value or a one-dimensional interpolation
C        is used; otherwise linear in p at IE and at IE+1, then
C        linear in E. For IT = 2 both E rows are taken at the same
C        index IP, i.e. the interpolation in E is done at fixed p/p0
C        (PP/p0 depends only on IP).
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
      DOUBLE PRECISION MT,GAMT,ASB
      COMMON/GrnPAR/MT,GAMT,ASB
C     --- 1. table and E interval
      GRN = DCMPLX(1D0,0D0)
      IF (EI.LT.EMAX1) THEN
         IT = 1
C     below the threshold table: extrapolate from the first interval
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
C     above the high-energy table: extrapolate from the last interval
C     (REVIEW.md, P1)
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
C     on-shell momentum p0 [GeV], centre of the IT = 2 p grid
         POS = DSQRT(MT*E)
         PMIN = PP (0,IT)
         PMAX = PP (NP2,IT)
C     reached only if EI is not a number (EMAX1 = EMIN2)
      ELSE
         WRITE (6,*) "IT"
         WRITE(6,*) EI,EMIN1,EMAX1,EMIN2,EMAX2
         STOP
      ENDIF
C
C     --- 2. p interval
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
C     --- 3. interpolation: on a grid point / linear in E / linear in
C     p / bilinear (p first, then E)
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
C     EE: E [GeV] of grid point IE (uniform grid) of table IT
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
C     EINT: index of the E interval [EE(IE),EE(IE+1)] containing E
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
C     PP: p [GeV] of grid point IP of table IT.
C       IT = 1 : uniform, PMIN1..PMAX1
C       IT = 2 : around p0 = POS (set by HOKAN2D), M = (NP2-1)/2 = 50,
C         IP <= M : PP = p0 (1-RR^(IP+1))/(1-RR^(M+1)), so PP(M) = p0
C                   and the spacing shrinks towards p0 (RR = 0.9)
C         IP >  M : PP = p0 + A RR^M (1-R1^k)/(1-R1), k = IP-M+1,
C                   R1 = 1/RR, A = (1-RR)/(1-RR^M) p0; the spacing
C                   grows away from p0
C       PP/p0 depends only on IP.
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
C     PINT: index IP of the p interval [PP(IP),PP(IP+1)] containing P,
C     obtained by inverting PP. Known issue P7 (docs/REVIEW.md): for
C     IT = 2 and p0 < P < p0 + A RR^M (k in (0,1) in the formula of
C     PP) the branch P > POS returns M-1 instead of M, so that HOKAN2D
C     extrapolates from the interval below p0 over a short distance.
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
