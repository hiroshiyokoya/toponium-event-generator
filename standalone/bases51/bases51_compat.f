C=======================================================================
C     Compatibility layer: lets the 2015 stand-alone code run with the
C     BASES/SPRING V5.1 of the CPC Program Library (S. Kawabata,
C     Comput. Phys. Commun. 88 (1995) 309, catalogue AAFW_v2_0).
C
C     The 2015 code was linked against a KEK V5.0 distribution that had
C       BSSETD, BSSETP  (J. Kanzaki, July '94) and
C       XHSAVE2         (H. Yokoya, tabulated histogram output),
C     and V5.1 calls a few CERN Program Library routines
C       DATIME, UCOPY, TIMEX, TIMEST.
C     They are re-implemented here. BASES/SPRING itself is not part of
C     this repository (CPC licence); see scripts/extract_bases51.sh.
C=======================================================================
C
C     BSSETD: integration variables (V5.0 API) -> /BPARM1/ of V5.1
C     (same as BSDIMS + BSGRID; only positive NDIM0/NWILD0 are taken)
      SUBROUTINE BSSETD (NDIM0,NWILD0,XL0,XU0,IG0)
      IMPLICIT NONE
      INTEGER MXDIM
      PARAMETER (MXDIM = 50)
      INTEGER NDIM0,NWILD0,IG0(NDIM0)
      DOUBLE PRECISION XL0(NDIM0),XU0(NDIM0)
      DOUBLE PRECISION XL,XU
      INTEGER NDIM,NWILD,IG,NCALL
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     .               IG(MXDIM),NCALL
      INTEGER I
      IF (NDIM0.GT.0)  NDIM  = NDIM0
      IF (NWILD0.GT.0) NWILD = NWILD0
      DO 100 I = 1, NDIM0
         XL(I) = XL0(I)
         XU(I) = XU0(I)
         IG(I) = IG0(I)
  100 CONTINUE
      RETURN
      END
C
C     BSSETP: sampling points, iterations and accuracies (V5.0 API)
C     -> /BPARM1/, /BPARM2/ of V5.1 (only positive values are taken)
      SUBROUTINE BSSETP (NCALL0,ITMX10,ITMX20,ACC10,ACC20)
      IMPLICIT NONE
      INTEGER MXDIM
      PARAMETER (MXDIM = 50)
      INTEGER NCALL0,ITMX10,ITMX20
      DOUBLE PRECISION ACC10,ACC20
      DOUBLE PRECISION XL,XU
      INTEGER NDIM,NWILD,IG,NCALL
      COMMON /BPARM1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     .               IG(MXDIM),NCALL
      DOUBLE PRECISION ACC1,ACC2
      INTEGER ITMX1,ITMX2
      COMMON /BPARM2/ ACC1,ACC2,ITMX1,ITMX2
      IF (NCALL0.GT.0)  NCALL = NCALL0
      IF (ITMX10.GT.0)  ITMX1 = ITMX10
      IF (ITMX20.GT.0)  ITMX2 = ITMX20
      IF (ACC10.GT.0D0) ACC1  = ACC10
      IF (ACC20.GT.0D0) ACC2  = ACC20
      RETURN
      END
C
C     XHSAVE2: write histogram ID as a table (x, dS/dx, error, and both
C     normalised to the total cross section) on unit LUNIT.
C     Ported from the author's V5.0 add-on; the COMMON blocks of V5.1
C     have the same layout.
      SUBROUTINE XHSAVE2 (LUNIT,ID)
      IMPLICIT NONE
      INTEGER LUNIT,ID
      DOUBLE PRECISION SCALLS,WGT,TI,TSI,TACC
      INTEGER IT
      COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT
      DOUBLE PRECISION AVGI,SD,CHI2A
      REAL STIME
      INTEGER ITG,ITF
      COMMON /BSRSLT/ AVGI,SD,CHI2A,STIME,ITG,ITF
      INTEGER NHS,NSC
      PARAMETER (NHS = 50, NSC = 50)
      INTEGER XHASH,DHASH,IFBASE,NHIST,MAPL,NSCAT,MAPD,NW
      COMMON /PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),
     .               NHIST, MAPL(4,NHS),
     .               NSCAT, MAPD(4,NSC),
     .               NW
      INTEGER IBUF
      REAL BUFF
      COMMON /PLOTB/ IBUF(281*NHS + 2527*NSC)
      DIMENSION BUFF(281*NHS + 2527*NSC)
      EQUIVALENCE (IBUF(1),BUFF(1))
      INTEGER IHIST,NTOTAL,IP1,IP2,IP3,IPF,IPF2,NXBIN,NX,I,J,K
      REAL XMIN,XMAX,DEV,FACT,TX,VLS,DEV2,VER,XX,TOT
      REAL TBL(0:2,51)

      TOT = 0.
      IF (NHIST.GT.0) THEN
         I = IABS(MOD(ID,13)) + 1
         IF (XHASH(1,I).EQ.1) THEN
            IF (ID.EQ.MAPL(1,XHASH(2,I))) THEN
               IHIST = XHASH(2,I)
               GO TO 200
            ENDIF
         ELSEIF (XHASH(1,I).GT.1) THEN
            DO 100 K = 2, XHASH(1,I)+1
               IF (ID.EQ.MAPL(1,XHASH(K,I))) THEN
                  IHIST = XHASH(K,I)
                  GO TO 200
               ENDIF
  100       CONTINUE
         ENDIF
      ENDIF
      WRITE (6,9000) ID
 9000 FORMAT(1X,'************ Warning from XHSAVE2 ************',
     .      /1X,' Histogram ID(',I5,' ) is illegal.',
     .      /1X,' This call is neglected.',
     .      /1X,'**********************************************')
      RETURN

  200 NTOTAL = SCALLS
      IP1   = MAPL(2,IHIST)
      XMIN  = BUFF(IP1)
      XMAX  = BUFF(IP1+1)
      NXBIN = IBUF(IP1+2)
      DEV   = BUFF(IP1+3)
      IP2   = MAPL(3,IHIST)
      IP3   = MAPL(4,IHIST)
      WRITE (6,9200) ID,LUNIT,(BUFF(I),I=IP3+1,IP3+15),NTOTAL,NXBIN,DEV
 9200 FORMAT(/1X,'** Histogram ID(',I5,' ) was saved in Unit(',I2,') **',
     .       /1X,'Title : ',15A4,
     .       /1X,'Entries     =',I10,
     .       /1X,'No. of bins =',I10,'  Width =',G13.4)
      IPF  = IP2 + 156
      IPF2 = IPF + 52
      FACT = 1./(NTOTAL*DEV)
      DO 400 I = 1, NXBIN
         TX  = BUFF(I+IPF)
         NX  = IBUF(I+IP2)
         VLS = TX*FACT
         IF (NX.GT.1) THEN
            DEV2 = NX*BUFF(I+IPF2)-TX*TX
            IF (DEV2.LE.0.0) THEN
               VER = 0.0
            ELSE
               VER = FACT*SQRT(DEV2/(NX-1))
            ENDIF
         ELSEIF (NX.EQ.1) THEN
            VER = VLS
         ELSE
            VER = 0.0
         ENDIF
         XX = XMIN + DEV*(FLOAT(I) - 1.)
         TBL(0,I) = XX
         TBL(1,I) = VLS
         TBL(2,I) = VER
         TOT = TOT + VLS*DEV
  400 CONTINUE
      DO 500 I = 1, NXBIN
         WRITE (LUNIT,9800) (TBL(J,I),J=0,2),(TBL(J,I)/AVGI,J=1,2)
  500 CONTINUE
      WRITE (LUNIT,9800) TBL(0,NXBIN)+DEV,(TBL(J,NXBIN),J=1,2),
     .                   (TBL(J,NXBIN)/AVGI,J=1,2)
 9800 FORMAT(1X,E11.4,4(2X,E14.7))
      RETURN
      END
C
C     ---- CERN Program Library routines used by V5.1 ----
C
C     DATIME: ID = yymmdd, IT = hhmm; also /SLATE/ IS(1:6) =
C     year, month, day, hour, minute, second (as in KERNLIB Z007)
      SUBROUTINE DATIME (ID,IT)
      IMPLICIT NONE
      INTEGER ID,IT
      INTEGER IS
      COMMON /SLATE/ IS(40)
      INTEGER V(8)
      CALL DATE_AND_TIME (VALUES=V)
      IS(1) = V(1)
      IS(2) = V(2)
      IS(3) = V(3)
      IS(4) = V(5)
      IS(5) = V(6)
      IS(6) = V(7)
      ID = MOD(V(1),100)*10000 + V(2)*100 + V(3)
      IT = V(5)*100 + V(6)
      RETURN
      END
C
C     UCOPY: copy N words from A to B (KERNLIB V301)
      SUBROUTINE UCOPY (A,B,N)
      IMPLICIT NONE
      INTEGER N,A(*),B(*),I
      DO 10 I = 1, N
         B(I) = A(I)
   10 CONTINUE
      RETURN
      END
C
C     TIMEST: start the CPU-time clock (the time limit is ignored)
C     TIMEX : CPU time in seconds since TIMEST (KERNLIB Z007)
      SUBROUTINE TIMEST (TLIM)
      IMPLICIT NONE
      REAL TLIM,T0
      COMMON /TBSTIM/ T0
      CALL CPU_TIME (T0)
      RETURN
      END
C
      SUBROUTINE TIMEX (T)
      IMPLICIT NONE
      REAL T,T0,TNOW
      COMMON /TBSTIM/ T0
      CALL CPU_TIME (TNOW)
      T = TNOW - T0
      RETURN
      END
