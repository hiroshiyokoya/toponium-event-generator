      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,IVAL
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS,MTT,DSDM,ERR
      DOUBLE PRECISION MIN,MAX,DM
      DOUBLE PRECISION PI,SQ2
      PARAMETER (PI=3.141593D0,SQ2=1.41421356D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,GS,AS0,GF,RGWF,RGAMT
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
      INTEGER IBSS,ISPR,IDYN
      INTEGER ICUT
      DOUBLE PRECISION MWCUT
      COMMON/CUT/MWCUT,ICUT
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
      IBSS  = 2                 ! 1:VEGAS, 2:BASES
      ISPR  = 0                 ! SPRING Event Generation 0:OFF, 1:ON
      WRITE (6,*) "IBSS, ISPR =", IBSS, ISPR
C.....
      IPDF  = 4                 ! 1:GRV98,2:Alekhin02,3:MRST2004NLO,4:CTEQ6
      IPORD = 1                 ! 1:LO, 2:NLO, 3:NNLO
      ISUB  = 0                 !
      IERR  = 0                 !
      CALL SETPDF (IPDF,IPORD,ISUB,IERR)
      IPRC  = 1                 ! 1:p-p, 2:p-pbar
      RS = 14D3                 ! LHC
C...  'run.inc'
      ebeam(1) = rs/2.
      ebeam(2) = rs/2.
      lpp(1) = 1
      lpp(2) = 1
C.....
      IPART = 2                 ! 0:Sum, 1:QQ, 2:GG
      ICLR  = 1                 ! 0:Sum, 1:Singlet, 2:Octet
      IGRN  = 1                 ! 0:No binding, 1:OCP,2:TAP
      ITR   = 0                 ! 0:
      INR   = 0                 ! 0:
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
      ASB = ASQCD (MUB,INF)
      GS = DSQRT (4D0*PI*AS0)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C
*      GAMT = 1.4911D0
      GF = gfermi
      GAMT = RGAMT (MT,MW,GF,0D0,MB)
*      GAMT = RGAMT (MT,MW,GF,0D0,0D0)
      WRITE (6,*) "GAM_T, GAM_W =", GAMT, GAMW
      TWIDTH = GAMT
C
      NDIM = 12
C
      MWCUT = 5D0
      MIN = 2D0*MB + MW + MWCUT
      MAX = RS
C
*      MIN = 2D0*MT - 15D0
*      MAX = 360D0
      MIN = 330D0
*      MIN = 400D0
      MAX = 3550D0
      DM = .5D0
*      MIN = 2D0*MT
*     MAX = 1000D0
*      DM = 1D0
      MTT = MIN
C
      IF ( IBSS.EQ.2 ) THEN
         CALL BSINIT
         CALL USERIN
      ENDIF
C
      MTT = 344D0
*      DO 99 I = 1, 500
         IF (IDYN.EQ.1) THEN
            MUR = MTT
            MUF = MUR
            AS0 = ASQCD (MUR,INF)
            GS = DSQRT (4D0*PI*AS0)
            GG(1) = DCMPLX (GS,0D0)
            GG(2) = GG(1)
            G = -GS
         ENDIF
         WRITE (2,*) "MTT =", MTT
         CALL SIGBWBLV (RS,MTT,DSDM,ERR,IBSS,ISPR,NDIM)
         WRITE (6,*) MTT, DSDM, ERR
         WRITE (7,*) MTT, DSDM
         IF (MTT.GE.2D0*MT-15D0) DM = .5D0
         IF (MTT.GE.2D0*MT-10D0) DM = .2D0
         IF (MTT.GE.2D0*MT-3D0) DM = .1D0
         IF (MTT.GE.2D0*MT+4D0) DM = 2D0
         IF (MTT.GE.400D0) DM = 5D0
         IF (MTT.GE.500D0) DM = 20D0
         IF (MTT.GE.1000D0) DM = 50D0
         IF ( MTT.GE.MAX ) GOTO 100
         MTT = MTT + DM
 99   CONTINUE
 100  CONTINUE
      STOP
      END
C
C
      SUBROUTINE SIGBWBLV (RS,MTT,DSDM,ERR,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER I,IVAL,JVAL
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS,MTT,B0,DSDM,ERR
      DOUBLE PRECISION GF,GFN,INT12
      EXTERNAL GFN,INT12
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS,VMTT
      INTEGER IBSS,ISPR,JBSS,JJ
      COMMON/V12/VRS,VMTT,JVAL,JBSS,JJ
      INTEGER IH
      COMMON/HST/IH
      include 'parameter.inc'
      IH = 12
      DSDM = 0D0
      ERR = 0D0
      VRS = RS
      VMTT = MTT
      JVAL = IVAL
      JBSS = IBSS
      JJ = 9
      NVP = 10000
      IF (IBSS.EQ.1) THEN
         CALL VEGAS (INT12,1D-3,IVAL,NVP,8,0,0)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BASES (INT12,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (2)
         CALL BHPLOT (2)
*         CALL XHSAVE (3,24)
         IF (ISPR.EQ.1) THEN
            DO I = 1,1000
               CALL SPRING (INT12,50)
            ENDDO
            CALL SPINFO (2)
            CALL SHPLOT (2)
         ENDIF
      ENDIF
      DSDM = S1
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
      nwild = ndim
      do idim = 1, ndim
         xl(idim) = 0d0
         xu(idim) = 1d0
         ig(idim) = 1
      enddo
      call bssetd (ndim,nwild,xl,xu,ig)
      ncall = 25000
      itmx1 = 20
      itmx2 = 20
      acc1  = 1D0
      acc2  = 1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
      do i = 1,ndim
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      return
      end
c
      SUBROUTINE MKHIST (IS,IVAL,IBSS,SIG,X)
      IMPLICIT NONE
      INTEGER I,IS,IVAL,IBSS
      DOUBLE PRECISION SIG,X(12)
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MT1,MT2,MW1
      COMMON/BWM/MT1,MT2,MW1
      DOUBLE PRECISION P(0:3,7),PH(0:3,7)
      COMMON/MOM/P
      DOUBLE PRECISION PT(0:3),PTB(0:3),PTH(0:3)
      DOUBLE PRECISION COST,PTPT,PABS,ETA
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      INTEGER IH
      COMMON/HST/IH
      IF (IH.EQ.IS .AND. IBSS.EQ.2) THEN
         CALL PCM2HCM (X1,X2,P,PH)
         DO I = 0,3
            PT (I) = P(I,3) + P(I,4)
            PTB(I) = P(I,5) + P(I,6)
            PTH (I) = PH(I,3) + PH(I,4)
         ENDDO
         COST = PT(3) / DSQRT(PT(1)**2+PT(2)**2+PT(3)**2)
         PTPT = DSQRT(PTH(1)**2+PTH(2)**2)
         PABS = DSQRT(PTH(1)**2+PTH(2)**2+PTH(3)**2)
         ETA = .5D0 * DLOG((PTH(0)+PTH(3))/(PTH(0)-PTH(3)))
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
      ENDIF
      RETURN
      END
