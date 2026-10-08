C=======================================================================
C     Optional run-time steering (added for the repository release).
C
C     If the file ./topbs.nml exists, the namelist /TOPBS/ in it overrides
C     the defaults hard-coded in the EG_*.f main programs. Entries that are
C     not given keep those defaults, so without the file the programs
C     behave exactly as the original 2015 code.
C
C       &TOPBS  IGRN=1, ITR=1, INR=0, IKF=1, RS=14000D0,
C               NCALL=100000, ITMX1=10, ITMX2=10, NEVENT=1000 /
C
C     IBSS : 0 VEGAS, 1 BASES (integrate + write grid to fort.55),
C            2 read the BASES grid from fort.55
C     ISPR : 1 generate unweighted events with SPRING
C     RS   : collider energy sqrt(s) [GeV] (pp)
C     IPART: 0 sum, 1 qq, 2 gg        ICLR: 0 sum, 1 singlet, 2 octet
C     IGRN, ITR, INR : see thr.inc    IDYN: 0 fixed, 1 dynamical scale
C     IKF  : 1 normalise to NLO (K factors of kfac / DATA KF, 14 TeV)
C     ICUT : 1 kinematical cut
C     NCALL, ITMX1, ITMX2 : BASES sampling points and iterations
C     NEVENT : number of unweighted events
C=======================================================================
      SUBROUTINE STEER (IBSS,ISPR,RS,IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,
     -                  ICUT)
      IMPLICIT NONE
      INTEGER IBSS,ISPR,IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,ICUT
      DOUBLE PRECISION RS
      INTEGER NCALL,ITMX1,ITMX2,NEVENT
      COMMON/STEERC/NCALL,ITMX1,ITMX2,NEVENT
      NAMELIST /TOPBS/ IBSS,ISPR,RS,IPART,ICLR,IGRN,ITR,INR,IDYN,IKF,
     -                 ICUT,NCALL,ITMX1,ITMX2,NEVENT
      LOGICAL EX
      INTEGER IU
      PARAMETER (IU=97)
      NCALL  = -1
      ITMX1  = -1
      ITMX2  = -1
      NEVENT = -1
      INQUIRE (FILE='topbs.nml',EXIST=EX)
      IF (.NOT.EX) RETURN
      OPEN (IU,FILE='topbs.nml',STATUS='OLD')
      READ (IU,NML=TOPBS)
      CLOSE (IU)
      WRITE (6,*) "STEER: settings read from topbs.nml"
      WRITE (6,NML=TOPBS)
      RETURN
      END
C
C     BASES parameters: override the defaults set in USERIN
      SUBROUTINE STEERB (NC,IT1,IT2)
      IMPLICIT NONE
      INTEGER NC,IT1,IT2
      INTEGER NCALL,ITMX1,ITMX2,NEVENT
      COMMON/STEERC/NCALL,ITMX1,ITMX2,NEVENT
      IF (NCALL.GT.0) NC  = NCALL
      IF (ITMX1.GT.0) IT1 = ITMX1
      IF (ITMX2.GT.0) IT2 = ITMX2
      RETURN
      END
C
C     number of events: NDEF unless NEVENT is given
      INTEGER FUNCTION NEVSTR (NDEF)
      IMPLICIT NONE
      INTEGER NDEF
      INTEGER NCALL,ITMX1,ITMX2,NEVENT
      COMMON/STEERC/NCALL,ITMX1,ITMX2,NEVENT
      NEVSTR = NDEF
      IF (NEVENT.GT.0) NEVSTR = NEVENT
      RETURN
      END
