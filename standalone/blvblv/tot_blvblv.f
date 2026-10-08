      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS,SIGT,ERR
      DOUBLE PRECISION PI,SQ2
      PARAMETER (PI=3.141593D0,SQ2=1.41421356D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,AS0,GS,GF,RGWF,RGAMT
      EXTERNAL ASQCD, RGAMT
      INTEGER IPDF,IPORD,ISUB,IERR
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      INTEGER IPART,ICLR
      COMMON/ID/IPART,ICLR
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      DOUBLE PRECISION MUB,ASB
      COMMON/BOHR/MUB,ASB
      INTEGER IBSS,ISPR
      INTEGER IDYN
      COMMON/DYN/IDYN
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      INTEGER ICUT
      INTEGER IHEL
      COMMON/HEL/IHEL
      DOUBLE PRECISION MWCUT
      COMMON/CUT/MWCUT,ICUT
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
      IBSS  = 2                 ! 1:VEGAS, 2:BASES
      ISPR  = 0                 ! SPRING Event Generation 0:OFF, 1:ON
      WRITE (2,*) "IBSS, ISPR =", IBSS, ISPR
      WRITE (6,*) "IBSS, ISPR =", IBSS, ISPR
C.....
      IPDF  = 4                 ! 1:GRV98,2:Alekhin02,3:MRST2004NLO,4:CTEQ6
      IPORD = 1                 ! 1:LO, 2:NLO, 3:NNLO
      ISUB  = 0                 !
      IERR  = 0                 !
*     CALL SETPDF (IPDF,IPORD,ISUB,IERR)
*     CALL INITPDFSETBYNAME ('cteq6l.LHpdf') ! NLO as
      CALL INITPDFSETBYNAME ('cteq6ll.LHpdf') ! LO as
      CALL INITPDF (0)
      IPRC  = 1                 ! 1:p-p, 2:p-pbar
      RS = 14D3                 ! LHC
C...  'run.inc'
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      lpp(1) = 1
      lpp(2) = 1
C....
      IHEL = 0
C.....
      IPART = 0                 ! 0:Sum, 1:QQ, 2:GG
      ICLR  = 0                 ! 0:Sum, 1:Singlet, 2:Octet
      IGRN  = 0                 ! 0:No binding, 1:OCP,2:TAP
      ITR   = 0                 ! 0 WWIDTH is not included in running TWIDTH
      INR   = 0                 ! 0
      IDYN  = 0                 ! 0:Fixed scale, 1:Dynamical scale
      ICUT  = 0                 ! 0: no cut, 1: kinematical cut
      WRITE (2,*) "IPART,ICLR,IGRN,IDYN,ICUT =",
     -     IPART,ICLR,IGRN,IDYN,ICUT
      WRITE (6,*) "IPART,ICLR,IGRN,IDYN,ICUT =",
     -     IPART,ICLR,IGRN,IDYN,ICUT
C.....
      call setpara ('param_card_sm.dat',.true.)
      MT = TMASS
      MB = BMASS
      MW = WMASS
      WRITE (2,*) "MT, MB, MW =", MT, MB, MW
      WRITE (6,*) "MT, MB, MW =", MT, MB, MW
      GAMT = TWIDTH
      GAMW = WWIDTH
C.....
      MUR = MT
      MUF = MUR
      MUB = 20D0
C.....
      INF = 5
      NF = DBLE (INF)
      AS0 = ASQCD (MUR,INF)
      IF (IDYN.EQ.0) THEN
         WRITE (2,*) "MUR, ASR =", MUR, AS0
         WRITE (6,*) "MUR, ASR =", MUR, AS0
      ENDIF
      ASB = ASQCD (MUB,INF)
      GS = DSQRT (4D0*PI*AS0)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C.....
      GF = gfermi
*      GAMT = RGAMT (MT,MW,GF,0D0,MB)
      WRITE (2,*) "GAM_T, GAM_W =", GAMT, GAMW
      WRITE (6,*) "GAM_T, GAM_W =", GAMT, GAMW
*      TWIDTH = GAMT
C
      NDIM = 16
C
      MWCUT = 5D0
      MIN = 2D0*MB + 2D0*MWCUT
      MAX = RS
*      MIN = 340D0
*      MAX = 346D0
C
      WRITE (2,*) "MIN, MAX =", MIN, MAX
      WRITE (6,*) "MIN, MAX =", MIN, MAX
C
      IF ( IBSS.EQ.2 ) THEN
         CALL BSINIT
         CALL USERIN
      ENDIF
      CALL SIGBWBLV (RS,SIGT,ERR,IBSS,ISPR,NDIM)
      WRITE (2,*) "SIGMA =", SIGT, " (pb) +-", ERR
      WRITE (6,*) "SIGMA =", SIGT, " (pb) +-", ERR
      STOP
      END
C
C
      SUBROUTINE SIGBWBLV (RS,SIG,ERR,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER I,IVAL,JVAL
      DOUBLE PRECISION RS,SIG,ERR
      DOUBLE PRECISION INT16
      EXTERNAL INT16
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS
      INTEGER IBSS,ISPR,JBSS
      COMMON/V16/VRS,JVAL,JBSS
      include 'parameter.inc'
      INTEGER IH
      COMMON/HST/IH
      IH = 16
      SIG = 0D0
      ERR = 0D0
      VRS = RS
      JVAL = IVAL
      JBSS = IBSS
      NVP = 25000
      IF (IBSS.EQ.1) THEN
         CALL VEGAS (INT16,1D-4,IVAL,NVP,10,0,0)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BASES (INT16,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (2)
         CALL BHPLOT (2)
*         CALL XHSAVE2 (11,33)
*         CALL XHSAVE2 (12,34)
*         CALL XHSAVE2 (13,36)
         IF (ISPR.EQ.1) THEN
            DO I = 1,1000
               CALL SPRING (INT16,50)
            ENDDO
            CALL SPINFO (2)
            CALL SHPLOT (2)
         ENDIF
      ENDIF
      SIG = S1
      ERR = S2
      RETURN
      END
C
C
      subroutine userin
      implicit none
      include 'parameter.inc'
      integer maxdim
      parameter (maxdim=50)
      double precision xl(maxdim),xu(maxdim)
      integer ig(maxdim)
      integer ncall,itmx1,itmx2
      double precision acc1,acc2
      integer i,idim
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      nwild = ndim
      do idim = 1, ndim
         xl(idim) = 0d0
         xu(idim) = 1d0
         ig(idim) = 1
      enddo
      call bssetd (ndim,nwild,xl,xu,ig)
      ncall = 50000
      itmx1 = 25
      itmx2 = 25
      acc1  = 1D0
      acc2  = 1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
      do i = 1,16
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      call xhinit (20,300D0,800D0,50,'dS/dmtt')
      call xhinit (21,330D0,360D0,50,'dS/dmtt')
      return
      end
C
      SUBROUTINE MKHIST (IS,IVAL,IBSS,SIG,X)
      IMPLICIT NONE
      INTEGER I,IS,IVAL,IBSS
      DOUBLE PRECISION SIG,X(16)
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MR1,MR2,MR3,MR4
      COMMON/BWM/MR1,MR2,MR3,MR4
      DOUBLE PRECISION P(0:3,8),PH(0:3,8)
      COMMON/MOM/P
      DOUBLE PRECISION PT(0:3),PTB(0:3)
      DOUBLE PRECISION PTL,PTM
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      DOUBLE PRECISION RBL1,RBL2,RLV
      DOUBLE PRECISION R2
      EXTERNAL R2
      INTEGER IH
      COMMON/HST/IH
      IF (IBSS.EQ.2 .AND. IS.EQ.IH) THEN
*         CALL PCM2HCM (X1,X2,P,PH)
*         DO I = 0,3
*            PT (I) = P(I,3) + P(I,4) + P(I,5)
*            PTB(I) = P(I,6) + P(I,7) + P(I,8)
*         ENDDO
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
         CALL XHFILL (20,MTT,SIG)
         CALL XHFILL (21,MTT,SIG)
      ENDIF
      RETURN
      END
