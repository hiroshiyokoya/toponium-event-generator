      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I,IVAL
      DOUBLE PRECISION MT,MB,MW,GAMT
      COMMON/MASS/MT,MB,MW,GAMT
      DOUBLE PRECISION RS
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      DOUBLE PRECISION PI,SQ2
      PARAMETER (PI=3.141593D0,SQ2=1.41421356D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,GS,ASR,GF,RGWF,RGAMT
      EXTERNAL ASQCD, RGAMT
      INTEGER IPDF,IPORD,ISUB,IERR
      INTEGER IPRC,IEXP
      COMMON/IP/IPRC,IEXP
      DOUBLE PRECISION MUR,MUF
      COMMON/SCL/MUR,MUF
      DOUBLE PRECISION MUB,ASB
      COMMON/BOHR/MUB,ASB
      INTEGER IBSS,ISPR
      INTEGER IDYN
      COMMON/DYN/IDYN
      INTEGER ICUT
      COMMON/CUT/ICUT
      INTEGER IPART,ICLR,FPROC(4)
      DATA FPROC/1,1,1,1/
      INTEGER NPART,OPROC(4)
      COMMON/PRC/NPART,OPROC
      INTEGER IKF
      COMMON/KFAC/KF,IKF
      include 'kfac.inc'
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
      IBSS  = 1                 ! 0:VEGAS, 1:BASES, 2:BSREAD
      ISPR  = 1                 ! SPRING Event Generation 0:OFF, 1:ON
      WRITE (50+IBSS,*) "IBSS, ISPR =", IBSS, ISPR
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
*     RS = 1.96D3               ! Tevatron
      RS = 7D3                  ! LHChalf
*     RS = 14D3                 ! LHC
      WRITE (50+IBSS,*) "IPRC, RS =", IPRC, RS
      WRITE (6,*) "IPRC, RS =", IPRC, RS
C...  'run.inc'
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      lpp(1) = 1
      lpp(2) = 1
C.....
      IPART = 0                 ! 0:Sum,1:qq,2:gg
      ICLR  = 0                 ! 0:Sum,1:singlet,2:octet
      IGRN  = 0                 ! 0:No binding, 1:Green function
      ITR   = 0                 ! 0:Naive, 1:Time Retardation
      INR   = 0                 ! 0:NO NON-RES DIGRAMS, 1:INCLUDE
      IDYN  = 0                 ! 0:Fixed scale,1:Dynamical scale
      IKF   = 0                 ! 0:No, 1:Yes
      ICUT  = 0                 ! 0: no cut, 1: knematical cut
      CALL STEER (IBSS,ISPR,RS,IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT)  ! optional ./topbs.nml
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      WRITE (50+IBSS,*) "IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT =",
     -     IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT
      WRITE (6,*) "IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT =",
     -     IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT
C.....
      IF (IPART.EQ.1) THEN
         FPROC(3) = 0
         FPROC(4) = 0
      ELSEIF (IPART.EQ.2) THEN
         FPROC(1) = 0
         FPROC(2) = 0
      ENDIF
      IF (ICLR.EQ.1) THEN
         FPROC(1) = 0
         FPROC(2) = 0
         FPROC(4) = 0
      ELSEIF (ICLR.EQ.2) THEN
         FPROC(3) = 0
      ENDIF
      NPART = 0
      DO I = 1, 4
         IF (FPROC(I).EQ.1) THEN
            NPART = NPART + 1
            OPROC(NPART) = I
         ENDIF
      ENDDO
      IF (NPART.EQ.0) THEN
         WRITE (6,*) "No process"
         STOP
      ENDIF
C.....
      call setpara ('param_card_sm.dat',.true.)
      MT = TMASS
      MB = BMASS
      MW = WMASS
C.....
      WRITE (50+IBSS,*) "MT, MB, MW =", MT, MB, MW
      WRITE (6,*) "MT, MB, MW =", MT, MB, MW
      GAMT = TWIDTH
C.....
      MUR = MT
      MUF = MUR
      MUB = 20D0
      SCALE = MUR
C.....
      INF = 5
      NF = DBLE (INF)
      ASR = ASQCD (MUR,INF)
      IF (IDYN.EQ.0) THEN
         FIXED_REN_SCALE=.TRUE.
         FIXED_FAC_SCALE=.TRUE.
         SCALE = MUR
         Q2FACT(1) = MUF**2
         Q2FACT(2) = Q2FACT(1)
         ALFAS = ASR
         WRITE (50+IBSS,*) "MUR, ASR =", MUR, ASR
         WRITE (6,*) "MUR, ASR =", MUR, ASR
      ELSE
         FIXED_REN_SCALE=.FALSE.
         FIXED_FAC_SCALE=.FALSE.
         WRITE (50+IBSS,*) "Dynamical scale"
         WRITE (6,*) "Dynamical scale"
      ENDIF
      ASB = 0.1534D0            ! ASQCD (MUB,INF)
      GS = DSQRT (4D0*PI*ASR)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C.....
      GF = gfermi
*     GAMT = RGAMT (MT,MW,GF,0D0,MB)
      WRITE (50+IBSS,*) "GAM_T =", GAMT
      WRITE (6,*) "GAM_T =", GAMT, TWIDTH
*     TWIDTH = GAMT
C.....
      NDIM = 12
C.....
      MIN = 2D0 * (MB + MW) + 1D-4 ! God-Hand
      MAX = RS
*     MIN = 330D0
*     MAX = 380D0
      WRITE (50+IBSS,*) "MIN, MAX =", MIN, MAX
      WRITE (6,*) "MIN, MAX =", MIN, MAX
C.....
      IF ( IBSS.EQ.1 ) THEN
         CALL BSINIT
         CALL USERIN
      ELSEIF (IBSS.EQ.2) THEN
         CALL BSINIT
         CALL USERIN
         OPEN (55,FILE="fort.55",FORM="unformatted")
      ENDIF
C.....
      CALL SIGBWBW (RS,IBSS,ISPR,NDIM)
C.....
      STOP
      END
C
C
      SUBROUTINE SIGBWBW (RS,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER NEVSTR
      EXTERNAL NEVSTR
      INTEGER I,IVAL,JVAL
      DOUBLE PRECISION RS,SIG,ERR
      COMMON/BSS/SIG,ERR
      DOUBLE PRECISION INT12
      EXTERNAL INT12
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS,VMTT
      INTEGER IBSS,ISPR,JBSS
      COMMON/V10/VRS,JVAL,JBSS
      DOUBLE PRECISION MT1,MT2
      COMMON/BWM/MT1,MT2
      include 'parameter.inc'
      include 'thr.inc'
      INTEGER IH
      COMMON/HST/IH
      IH = 12
      SIG = 0D0
      ERR = 0D0
      VRS = RS
      JVAL = IVAL
      JBSS = IBSS
      NVP = 25000
      IF (IBSS.EQ.0) THEN
         CALL VEGAS (INT12,1D-4,IVAL,NVP,10,0,0)
      ELSEIF (IBSS.EQ.1) THEN
         CALL BASES (INT12,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (50+IBSS)
         CALL BHPLOT (50+IBSS)
         SIG = S1
         ERR = S2
         WRITE (50+IBSS,*) "SIGMA =", SIG, " (pb) +-", ERR
         WRITE (6,*) "SIGMA =", SIG, " (pb) +-", ERR
*         CALL XHSAVE2 (10,27)
*         CALL XHSAVE2 (11,28)
*         CALL XHSAVE2 (12,25)
*         CALL XHSAVE2 (13,21)
*         CALL XHSAVE2 (11,32)
*         CALL XHSAVE2 (12,33)
*         CALL XHSAVE2 (13,34)
         CALL BSWRIT (55)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BSREAD (55)
         CLOSE (55)
      ENDIF
      IF (ISPR.EQ.1) THEN
         CALL LHEOUT (-1)
         DO I = 1, NEVSTR(1000000)
            CALL SPRING (INT12,50)
            WRITE (31+IGRN,*) MT1,MT2
            CALL LHEOUT (0)
         ENDDO
         CALL LHEOUT (1)
         CALL SPINFO (54)
         CALL SHPLOT (54)
      ENDIF
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
      INTEGER NPART,OPROC(4)
      COMMON/PRC/NPART,OPROC
      nwild = ndim
      do idim = 1, ndim
         xl(idim) = 0d0
         xu(idim) = 1d0
         ig(idim) = 1
      enddo
      call bssetd (ndim,nwild,xl,xu,ig)
      ncall = 1 000 000
      itmx1 = 15
      itmx2 = 15
      call steerb (ncall,itmx1,itmx2)
      acc1  = .1D0
      acc2  = .1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
      do i = 1,10
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      call xhinit (11,0D0,1D0,NPART,'dS/dX')
      call xhinit (12,0D0,1D0,50,'dS/dX')
      CALL XHINIT (21,168D0,178D0,50,'MT1')
      CALL XHINIT (22,168D0,178D0,50,'MT2')
      CALL XHINIT (23,-1D0,1D0,50,'COST')
      CALL XHINIT (24,0D0,500D0,50,'PT_T')
      CALL XHINIT (25,0D0,1000D0,50,'P_T')
      CALL XHINIT (26,-5D0,5D0,50,'ETA_T')
      CALL XHINIT (27,330D0,380D0,50,'MTT')
      CALL XHINIT (28,300D0,800D0,50,'MTT')
      CALL XHINIT (29,MIN,MAX,50,'MTT')
      CALL DHINIT (31,166D0,180D0,50,166D0,180D0,50,'MT1 vs MT2')
      call xhinit (32,0D0,500D0,50,'dS/dMbb')
      call xhinit (33,0D0,5D0,50,'ds/dRbb')
      call xhinit (34,0D0,3.1416D0,50,'dS/dDeltaPhi_bb')
      return
      end
c
      SUBROUTINE MKHIST (IS,IVAL,IBSS,SIG,X)
      IMPLICIT NONE
      INTEGER I,IS,IVAL,IBSS
      DOUBLE PRECISION SIG,X(10)
      DOUBLE PRECISION X1,X2
      COMMON/XBJ/X1,X2
      DOUBLE PRECISION MT1,MT2
      COMMON/BWM/MT1,MT2
      DOUBLE PRECISION MTT
      COMMON/TT/MTT
      DOUBLE PRECISION P(0:3,6),PH(0:3,6)
      COMMON/MOM/P
      DOUBLE PRECISION PTOP(0:3),PTB(0:3),PTH(0:3)
      DOUBLE PRECISION COST,PTPT,PABS,ETAT
      DOUBLE PRECISION MBB,RBB,PHIBB
      INTEGER IH
      COMMON/HST/IH
      include 'kin_functions.inc'
      IF (IH.EQ.IS .AND. IBSS.GE.1) THEN
         CALL PCM2HCM (X1,X2,P,PH)
         DO I = 0,3
            PTOP(I) = P(I,3) + P(I,4)
            PTB (I) = P(I,5) + P(I,6)
            PTH (I) = PH(I,3) + PH(I,4)
         ENDDO
         COST = PTOP(3) / DSQRT(PTOP(1)**2+PTOP(2)**2+PTOP(3)**2)
         PTPT = DSQRT(PTH(1)**2+PTH(2)**2)
         PABS = DSQRT(PTH(1)**2+PTH(2)**2+PTH(3)**2)
         ETAT = .5D0 * DLOG((PTH(0)+PTH(3))/(PTH(0)-PTH(3)))
         MBB = DSQRT(SumDot(P(0,3),P(0,5),+1D0))
         RBB = DSQRT(    R2(P(0,3),P(0,5)))
         PHIBB = DELTA_PHI(P(0,3),P(0,5))
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
         CALL XHFILL (21,MT1,SIG)
         CALL XHFILL (22,MT2,SIG)
         CALL XHFILL (23,COST,SIG)
         CALL XHFILL (24,PTPT,SIG)
         CALL XHFILL (25,PABS,SIG)
         CALL XHFILL (26,ETAT,SIG)
         CALL XHFILL (27,MTT,SIG)
         CALL XHFILL (28,MTT,SIG)
         CALL XHFILL (29,MTT,SIG)
         CALL DHFILL (31,MT1,MT2,SIG)
         CALL XHFILL (32,MBB,SIG)
         CALL XHFILL (33,RBB,SIG)
         CALL XHFILL (34,PHIBB,SIG)
      ENDIF
      RETURN
      END
