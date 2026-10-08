      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION MT,MB,MW,GAMT,GAMW
      COMMON/MASS/MT,MB,MW,GAMT,GAMW
      DOUBLE PRECISION RS
      DOUBLE PRECISION PI,SQ2
      PARAMETER (PI=3.141593D0,SQ2=1.41421356D0)
      INTEGER INF
      DOUBLE PRECISION NF,ASQCD,ASR,GS,GF,RGWF,RGAMT
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
      DOUBLE PRECISION MIN,MAX
      COMMON/INV/MIN,MAX
      INTEGER ICUT
      DOUBLE PRECISION MWCUT1,MWCUT2
      COMMON/CUT/MWCUT1,MWCUT2,ICUT
      INTEGER IPART,ICLR,FPROC(4)
      DATA FPROC/1,1,1,1/
      INTEGER NPART,OPROC(4)
      COMMON/PRC/NPART,OPROC
      INTEGER IKF
      DOUBLE PRECISION KF(3)
      COMMON/KFAC/KF,IKF
      DATA KF/1.14D0,1.39D0,1.16D0/ ! LHC 14 TeV
*     DATA KF/1.31D0,1.60D0,1.18D0/ ! LHC  7 TeV
      include 'thr.inc'
      include 'coupl_sm.inc'
      include 'parameter.inc'
      include 'run.inc'
C.....
      IBSS  = 2                 ! 0:VEGAS, 1:BASES, 2:BSREAD
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
*     RS = 7D3                  ! LHChalf
      RS = 14D3                 ! LHC
      WRITE (6,*) "IPRC, RS =", IPRC, RS
C...  'run.inc'
      ebeam(1) = rs/2D0
      ebeam(2) = rs/2D0
      lpp(1) = 1
      lpp(2) = 1
C.....
      IPART = 0                 ! 0:Sum,1:qq,2:gg
      ICLR  = 0                 ! 0:Sum,1:singlet,2:octet
      IGRN  = 1                 ! 0:No binding, 1:OCP,2:TAP
      ITR   = 1                 ! 0:Naive, 1:TR scheme
                                ! WWIDTH is not included in running TWIDTH
      INR   = 1                 ! 0:No,1:FWS
      IDYN  = 0                 ! 0:Fixed scale,1:Dynamical scale
      IKF   = 1                 ! 0:No,1:Yes
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
      WRITE (50+IBSS,*) "MT, MB, MW =", MT, MB, MW
      WRITE (6,*) "MT, MB, MW =", MT, MB, MW
      GAMT = TWIDTH
      GAMW = WWIDTH
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
      ASB = ASQCD (MUB,INF)
      GS = DSQRT (4D0*PI*ASR)
      GG(1) = DCMPLX (GS,0D0)
      GG(2) = GG(1)
      G = -GS
C.....
      GF = gfermi
*      GAMT = RGAMT (MT,MW,GF,0D0,MB)
      WRITE (50+IBSS,*) "GAM_T, GAM_W =", GAMT, GAMW
      WRITE (6,*) "GAM_T, GAM_W =", GAMT, GAMW
*      TWIDTH = GAMT
C.....
      NDIM = 18
C.....
      MWCUT1 = MW + 0.1D0       ! to avoid vanishing Gamma_t
      MWCUT2 = 0.1D0
C.....
      MIN = 2D0*MB + 2D0*MWCUT1
      MAX = RS
*      MIN = 340D0
*      MAX = 346D0
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
      CALL SIGBJJBJJ (RS,IBSS,ISPR,NDIM)
C.....
      STOP
      END
C
      SUBROUTINE SIGBJJBJJ (RS,IBSS,ISPR,IVAL)
      IMPLICIT NONE
      INTEGER NEVSTR
      EXTERNAL NEVSTR
      INTEGER I,J,K,IVAL,JVAL
      DOUBLE PRECISION RS,SIG,ERR
      COMMON/BSS/SIG,ERR
      DOUBLE PRECISION INT18
      EXTERNAL INT18
      DOUBLE PRECISION S1,S2,S3,S4
      COMMON/RESULT/S1,S2,S3,S4
      INTEGER NVP
      DOUBLE PRECISION VRS
      INTEGER IBSS,ISPR,JBSS
      COMMON/V16/VRS,JVAL,JBSS
      include 'parameter.inc'
      INTEGER IH
      COMMON/HST/IH
      IH = 18
      SIG = 0D0
      ERR = 0D0
      VRS = RS
      JVAL = IVAL
      JBSS = IBSS
      NVP = 25000
      IF (IBSS.EQ.0) THEN
         CALL VEGAS (INT18,1D-4,IVAL,NVP,10,0,0)
      ELSEIF (IBSS.EQ.1) THEN
         CALL BASES (INT18,S1,S2,CTIME,IT1,IT2)
         CALL BSINFO (50+IBSS)
         CALL BHPLOT (50+IBSS)
         SIG = S1
         ERR = S2
         WRITE (50+IBSS,*) "SIGMA =", SIG, " (pb) +-", ERR
         WRITE (6,*) "SIGMA =", SIG, " (pb) +-", ERR
         CALL BSWRIT (55)
      ELSEIF (IBSS.EQ.2) THEN
         CALL BSREAD (55)
         CLOSE (55)
      ENDIF
      IF (ISPR.EQ.1) THEN
         CALL LHEOUT (-1)
         DO I = 1, NEVSTR(100000)
            CALL SPRING (INT18,50)
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
      ncall = 100 000 000
      itmx1 = 10
      itmx2 = 10
      call steerb (ncall,itmx1,itmx2)
      acc1  = .1D0
      acc2  = .1D0
      call bssetp (ncall,itmx1,itmx2,acc1,acc2)
      do i = 1,16
         call xhinit (i,0D0,1D0,50,'dS/dX')
      enddo
      call xhinit (17,0D0,1D0,NPART,'dS/dX')
      call xhinit (18,0D0,1D0,16,'dS/dX')
      call xhinit (21,300D0,800D0,50,'dS/dMtt')
      call xhinit (22,330D0,360D0,50,'dS/dMtt')
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
      IF (IBSS.GE.1 .AND. IS.EQ.IH) THEN
         CALL PCM2HCM (X1,X2,P,PH)
         DO I = 0,3
            PT (I) = P(I,3) + P(I,4) + P(I,5)
            PTB(I) = P(I,6) + P(I,7) + P(I,8)
         ENDDO
         DO I = 1,IVAL
            CALL XHFILL (I,X(I),SIG)
         ENDDO
         CALL XHFILL (21,MTT,SIG)
         CALL XHFILL (22,MTT,SIG)
      ENDIF
      RETURN
      END
