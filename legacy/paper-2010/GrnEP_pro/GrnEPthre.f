      PROGRAM MAIN
      IMPLICIT NONE
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION MUB,ASB
      COMMON/GrnPOT/MUB,ASB
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      DOUBLE PRECISION XIN,XUV,XIR
      COMMON/GrnX/XIN,XUV,XIR
      DOUBLE PRECISION X,R,V,VQCD
      EXTERNAL VQCD
C...  External program
      INTEGER NF,NORD
      DOUBLE PRECISION ALPHAS
      DOUBLE PRECISION RGAMT
      DOUBLE PRECISION GF,MW,MB,VTB
      PARAMETER (GF=1.16639D-5,MW=80.4D0,MB=5D0,VTB=1D0)
*     PARAMETER (GF=1.17801D-5,MW=80D0,MB=4.7D0,VTB=0.998D0)
      DOUBLE PRECISION MZ,ASZ
      PARAMETER (MZ=91.2D0,ASZ=0.118D0)
*     PARAMETER (MZ=91.2D0,ASZ=0.12D0)
      DOUBLE PRECISION ALFST,ALFSP,RQCD
C...  Top Mass, Width
      MU = 173D0
      GAMT = 1.4911D0
*     MU = 176.17D0
*     GAMT = 1.602D0
*     GAMT = RGAMT (MU,GF,MW,MB,VTB) ! LO top-quark width
*     WRITE (6,*) MU,GAMT
C...  CALL alphas_pegasus.f
      NF = 5                    ! num of flavor
      NORD = 1                  ! 1:NLO,2:NNLO
      MUB = 20D0
      ASB = ALPHAS (MUB,MZ,ASZ,NF,NORD)
      WRITE (6,*) MUB,ASB
C...  or give them directly
*     MUB = 20D0
*     ASB = 0.1534D0
C...  Potential and Color setting
*     ICLR = 1                  ! Color of ttbar pair
      DO 100 ICLR = 1, 2        ! 1:Singlet,2:Octet
      WRITE (6,*) ICLR,MU,GAMT,ASB
*     CALL SETVPR (ASZ,MU,MB,MZ,GAMT,ALFST,ALFSP,RQCD)
C...  Setting Initial Conditions
      XIN = -1D0
      XUV = -7D0
      XIR =  2D0
C...  Check QCD Potential
      DO X = XUV, XIR, (XIR-XUV)/100.
         R = 10D0**X
         V = VQCD(R)
         WRITE (ICLR,*) R,V
      ENDDO
C.....
      CALL CALCGRNEP
C.....
 100  CONTINUE
      STOP
      END
C
C
      SUBROUTINE CALCGRNEP
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INCLUDE 'incgrnep.inc'
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER IE,IP,EINT,PINT
      DOUBLE PRECISION EE,PP
      EXTERNAL EE,PP,EINT,PINT
      DOUBLE PRECISION DIMB
      COMMON/GrnIMG/DIMB
C.....
      GA0(1) = DCMPLX( 1D0,0D0)
      GA0(2) = DCMPLX( 1D0,0D0)
      GB0(1) = DCMPLX( 1D0,0D0)
      GB0(2) = DCMPLX(-1D0,0D0)
C...  Energy Loop
      DO 50 IE = 0, NE1
         E = EE (IE)
*     WRITE (6,*) IE,E
C...
         CALL SOLVE (IE)
         CALL FOURIER (IE)
         WRITE (6,*) IE, E, " done, ", DIMB
C...
 50   CONTINUE
 51   CONTINUE
C...  WRITE TABLE TO THE FILE
      CALL SAVETABLE
C...
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE SOLVE (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      INTEGER I,J
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      INTEGER IE
      DOUBLE PRECISION EPS,ERR
      DOUBLE PRECISION H,X,R
      DOUBLE PRECISION X1,R0,R1
      DOUBLE PRECISION XIN,XUV,XIR
      COMMON/GrnX/XIN,XUV,XIR
      DOUBLE COMPLEX GA(N),GAN(N),DGA
      DOUBLE COMPLEX GB(N),GBN(N),DGB
      DOUBLE COMPLEX GP(N),GPN(N)
      DOUBLE PRECISION GPMN
      DOUBLE COMPLEX A1,A2,A5,A6,A7,A8, A1P,A2P,A7P,A8P,BP
      DOUBLE COMPLEX B,B1,B2,GP0(N)
      COMMON/GrnSOL/B1,B2,GP0
      DOUBLE COMPLEX G0N,G1N
      DOUBLE PRECISION DIMB
      DOUBLE PRECISION DEL,D0
      INTEGER I1,I2,I7,I8
      COMMON/GrnIMG/DIMB
C...  Define Accuracy
      EPS = 1D-4
      ERR = EPS/DBLE(NLOOPT)
C...  Take R->0 Limit
C...  Set Final Point of Evaluation
 90   CONTINUE
      X0 = XIN
      X1 = XUV
 91   CONTINUE
C...  Initialize
      X = X0
      R0 = 10D0**X0
      DO I=1,N
         GA(I) = GA0(I)
         GB(I) = GB0(I)
      ENDDO
C...  Reflesh Flag
      I1=0
      I2=0
      I7=0
      I8=1
C...  Start to solve
      H = (X1-X0)/NLOOPT        ! Step of Evaluations
      DO I=1,NLOOPT,1
         CALL RUNGE_KUTTA (X,H,GA,GAN)
         CALL RUNGE_KUTTA (X,H,GB,GBN)
         X = X+H
         R = 10D0**X
         DGA = 1D0/L10/R * GAN(2)
         DGB = 1D0/L10/R * GBN(2)
         A1 = CONE/(DGA - DGB*GAN(1)/GBN(1))
         A2 = -A1*GAN(1)/GBN(1)
         A5 = CONE/GAN(1)
         A6 = DCMPLX(0D0,-DIMAG(A5*DGA))
         A7 = A5
         IF ( I1.EQ.0 .AND. ABS((A1-A1P)/A1P).LE.ERR ) I1 = 1
         IF ( I2.EQ.0 .AND. ABS((A2-A2P)/A2P).LE.ERR ) I2 = 1
         IF ( I7.EQ.0 .AND. ABS((A7-A7P)/A7P).LE.ERR ) I7 = 1
         IF ( I1*I2*I7 .EQ. 1 ) THEN
            RUV = R
            WRITE (11,*) IE, X, "; R->0 converge"
            GOTO 100
         ENDIF
         DO J=1,N
            GA(J) = GAN(J)
            GB(J) = GBN(J)
         ENDDO
         A1P = A1
         A2P = A2
         A7P = A7
      ENDDO
      WRITE (11,*) IE, X1, "; R->0 NOT CONVERGE"
      X1 = X1 - 1D0
      GOTO 91
 100  CONTINUE
C...  Take R->Infinity Limit
C...  Set Final Point of Evaluation
      X1 = XIR
 101  CONTINUE
C...  Initialize
      X = X0
      DO I=1,N
         GA(I) = GA0(I)
         GB(I) = GB0(I)
      ENDDO
C...  Start to solve
      H = (X1-X0)/NLOOPT        ! Step of Evaluations
      DO I=1,NLOOPT,1
         CALL RUNGE_KUTTA (X,H,GA,GAN)
         CALL RUNGE_KUTTA (X,H,GB,GBN)
         X = X+H
         R = 10D0**X
         DGA = 1D0/L10/R * GAN(2)
         DGB = 1D0/L10/R * GBN(2)
         G0N = A1*GAN(1) + A2*GBN(1)
         G1N = A5*GAN(1) + A6*G0N
         A8 = - A7*GAN(1)/G0N
         B = -G1N/G0N
         DEL = ABS((A8-A8P)/A8P)
         IF ( DEL.LE.D0 ) THEN
            I8 = I8*2
            D0 = DEL
         ELSE
            I8 = 1
            D0 = ERR
         ENDIF
         IF ( I8.GT.1000 ) THEN
            RIR = R
            WRITE (11,*) IE, X, "; R->oo converge"
            GOTO 110
         ENDIF
         DO J=1,N
            GA(J) = GAN(J)      ! Set GG for the next step
            GB(J) = GBN(J)      ! Set GG for the next step
         ENDDO
         A8P = A8
         BP = B
      ENDDO
      WRITE (11,*) IE, X1, "; R->oo, NOT CONVERGE"
      NLOOPT = NLOOPT + 5000
      X1 = X1 + .2D0
      GOTO 101
 110  CONTINUE
      WRITE (11,*) IE, E, ": Runge-Kutta Solved ", NLOOPT
      WRITE (11,*) "R_0 =", R0
      B1 = A7 + A1*A8
      B2 =      A2*A8
      GP0(1) = B1*GA0(1) + B2*GB0(1)
      GP0(2) = B1*GA0(2) + B2*GB0(2)
      DIMB = DIMAG(B)
      WRITE (12,*) E, R0, GP0
      WRITE (30+ICLR,*) E, DIMB ! PRINT Im{G}
C...  Check End Points
C...  Re-calc R_UV where the function sufficiently reach to (1.,0.)
C...  Set Final Point of Evaluation
      R1 = RUV
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  R -> R_IR
      H = (X1-X0)/NLOOP1        ! Step of Evaluations
      DO I=1,NLOOP1,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
         IF ( ABS(DREAL(GPN(1))-1D0).LT.EPS .AND.
     -        ABS(DIMAG(GPN(1))-0D0).LT.EPS ) THEN
            WRITE (11,*) "R_UV", RUV, " ->", R
            RUV = R
            GOTO 10
         ENDIF
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
 10   CONTINUE
C...  Re-calc R_IR as the minimum point of ABS(GP)
C...  Set Final Point of Evaluation
      R1 = RIR * 1.5D0
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
      GPMN = ABS(GP0(1))
C...  R -> R_IR
      H = (X1-X0)/NLOOP1        ! Step of Evaluations
      DO I=1,NLOOP1,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
         IF ( ABS(GPN(1)) .LT. GPMN ) THEN
            R0 = R
            GPMN = DBLE(ABS(GPN(1)))
         ENDIF
         IF ( ABS(GPN(1)) .LT. 1D-4 ) THEN
            R0 = R
            GOTO 95
         ENDIF
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
 95   CONTINUE
      WRITE (11,*) "R_IR:", RIR, " ->", R0
      RIR = R0
      R0 = 10D0**X0
      WRITE (13,*) E, RUV, R0, RIR
      CALL REFINE (IE)
      RETURN
      END
C
      SUBROUTINE REFINE (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      INTEGER I,J
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER IE
      DOUBLE PRECISION H,X,R
      DOUBLE PRECISION R0,X1,R1
      DOUBLE COMPLEX GP(N),GPN(N)
      DOUBLE COMPLEX B1,B2,GP0(N)
      COMMON/GrnSOL/B1,B2,GP0
C...  1: R->0 direction
C...  Set Final Point of Evaluation
      R1 = RUV
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  Store starting point data
      R0 = 10D0**X0
      REGP1 (0) = DREAL(GP0(1))
      IMGP1 (0) = DIMAG(GP0(1))
      RP1 (0) = R0
      REGP2 (0) = DREAL(GP0(1))
      IMGP2 (0) = DIMAG(GP0(1))
      RP2 (0) = R0
*     WRITE (14,*) R0,DREAL(GP0(1)),DIMAG(GP0(1)),ABS(GP0(1))
C...  R -> R_UV
      H = (X1-X0)/NLOOP1        ! Step of Evaluations
      DO I=1,NLOOP1,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
C...
         REGP1 (I) = DREAL(GPN(1))
         IMGP1 (I) = DIMAG(GPN(1))
         RP1 (I) = R
*         WRITE (14,*) R,DREAL(GPN(1)),DIMAG(GPN(1)),ABS(GPN(1))
C...
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
C...  2: R->oo direction
C...  Set Final Point of Evaluation
      R1 = RIR
      X1 = LOG10(R1)
C...  Initialize
      X = X0
      DO I=1,N
         GP(I) = GP0(I)
      ENDDO
C...  R -> R_IR
      H = (X1-X0)/NLOOP1        ! Step of Evaluations
      DO I=1,NLOOP1,1
         CALL RUNGE_KUTTA (X,H,GP,GPN)
         X = X+H
         R = 10D0**X
C...
         REGP2 (I) = DREAL(GPN(1))
         IMGP2 (I) = DIMAG(GPN(1))
         RP2 (I) = R
*         WRITE (15,*) R, DREAL(GPN(1)),DIMAG(GPN(1)),ABS(GPN(1))
C...
         DO J=1,N
            GP(J) = GPN(J)
         ENDDO
      ENDDO
C...
      RETURN
      END
C
      SUBROUTINE FOURIER (IE)
      IMPLICIT NONE
      INCLUDE 'incgrnep.inc'
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE COMPLEX CZERO,CONE,CIMG
      PARAMETER (CZERO=(0D0,0D0),CONE=(1D0,0D0), CIMG=(0D0,1D0))
      INTEGER I,IE,IP
      DOUBLE PRECISION P,PP
      EXTERNAL PP
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION DL1,DL2,SIGRE,SIGIM,SIGUV
      DOUBLE PRECISION SIGRE1,SIGIM1,SIGRE2,SIGIM2
      DOUBLE PRECISION S1
      EXTERNAL S1
      DOUBLE PRECISION SMP(2)
      DATA SMP/2D0,4D0/         ! /Even, Odd/
      DOUBLE COMPLEX GAMMA
      DOUBLE COMPLEX GRNEP(0:NE1,0:NP1)
      COMMON/GrnMOM/GRNEP
      DL1 = DABS(DLOG(RP1(1))-DLOG(RP1(0)))
      DL2 = DABS(DLOG(RP2(1))-DLOG(RP2(0)))
      DO 100 IP = 0,NP1,1
         P = PP (IP)
C...  Daikei Koushiki
*         SIGRE1 = 0D0
*         SIGRE2 = 0D0
*         SIGIM1 = 0D0
*         SIGIM2 = 0D0
*         DO I = 0,NLOOP1-1,1
*            SIGRE1 = SIGRE1 + DL1 * ( RP1(I)*REGP1(I)*S1(P,RP1(I),MU)
*     -           + RP1(I+1)*REGP1(I+1)*S1(P,RP1(I+1),MU) ) / 2D0
*            SIGRE2 = SIGRE2 + DL2 * ( RP2(I)*REGP2(I)*S1(P,RP2(I),MU)
*     -           + RP2(I+1)*REGP2(I+1)*S1(P,RP2(I+1),MU) ) / 2D0
*            SIGIM1 = SIGIM1 + DL1 * ( RP1(I)*IMGP1(I)*S1(P,RP1(I),MU)
*     -           +  RP1(I+1)*IMGP1(I+1)*S1(P,RP1(I+1),MU) ) / 2D0
*            SIGIM2 = SIGIM2 + DL2 * ( RP2(I)*IMGP2(I)*S1(P,RP2(I),MU)
*     -           + RP2(I+1)*IMGP2(I+1)*S1(P,RP2(I+1),MU) ) / 2D0
*         ENDDO
*         WRITE (6,*) SIGRE1,SIGRE2,SIGIM1,SIGIM2
C...  Simpthon Integral Formula
         SIGRE1 = DL1/3D0 * ( RP1(0)*REGP1(0)*S1(P,RP1(0),MU)
     -        + RP1(NLOOP1)*REGP1(NLOOP1)*S1(P,RP1(NLOOP1),MU) )
         SIGRE2 = DL2/3D0 * ( RP2(0)*REGP2(0)*S1(P,RP2(0),MU)
     -        + RP2(NLOOP1)*REGP2(NLOOP1)*S1(P,RP2(NLOOP1),MU) )
         SIGIM1 = DL1/3D0 * ( RP1(0)*IMGP1(0)*S1(P,RP1(0),MU)
     -        + RP1(NLOOP1)*IMGP1(NLOOP1)*S1(P,RP1(NLOOP1),MU) )
         SIGIM2 = DL2/3D0 * ( RP2(0)*IMGP2(0)*S1(P,RP2(0),MU)
     -        + RP2(NLOOP1)*IMGP2(NLOOP1)*S1(P,RP2(NLOOP1),MU) )
         DO I = 1,NLOOP1-1, 1
            SIGRE1 = SIGRE1 + DL1/3D0 *
     -           SMP(MOD(I,2)+1)*RP1(I)*REGP1(I)*S1(P,RP1(I),MU)
            SIGRE2 = SIGRE2 + DL2/3D0 *
     -           SMP(MOD(I,2)+1)*RP2(I)*REGP2(I)*S1(P,RP2(I),MU)
            SIGIM1 = SIGIM1 + DL1/3D0 *
     -           SMP(MOD(I,2)+1)*RP1(I)*IMGP1(I)*S1(P,RP1(I),MU)
            SIGIM2 = SIGIM2 + DL2/3D0 *
     -           SMP(MOD(I,2)+1)*RP2(I)*IMGP2(I)*S1(P,RP2(I),MU)
         ENDDO
*         WRITE (6,*) SIGRE1,SIGRE2,SIGIM1,SIGIM2
C...
         IF (P.EQ.0D0) THEN
            SIGUV = P*MU*RUV**2/2D0
         ELSE
            SIGUV = MU/P**2*(1D0-DCOS(P*RUV))
         ENDIF
         SIGRE = SIGRE1 + SIGRE2 + SIGUV
         SIGIM = SIGIM1 + SIGIM2
*         WRITE (6,*) P,SIGRE1, SIGRE2, SIGUV
C...
         GAMMA = - (E + CIMG*GAMT - P**2/MU) * DCMPLX(SIGRE,SIGIM)
         GRNEP(IE,IP) = GAMMA
*         WRITE (6,*) PP(IP), DREAL(GAMMA),DIMAG(GAMMA), IE
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE SAVETABLE
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      DOUBLE PRECISION MUB,ASB
      COMMON/GrnPOT/MUB,ASB
      INTEGER IE,IP,IOF
      PARAMETER (IOF=15)
      DOUBLE COMPLEX GRNEP(0:NE1,0:NP1)
      COMMON/GrnMOM/GRNEP
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      INCLUDE 'filedata.inc'
      IF (ICLR.EQ.1) THEN
         WRITE (6,*) "   Open file : ", T1
         OPEN (IOF,FILE=T1)
      ELSEIF (ICLR.EQ.2) THEN
         WRITE (6,*) "   Open file : ", T8
         OPEN (IOF,FILE=T8)
      ENDIF
      WRITE (IOF,*) MU, GAMT, ASB
      WRITE (IOF,*) NE1, NP1
      WRITE (IOF,*) EMIN1, EMAX1
      WRITE (IOF,*) PMIN1, PMAX1
      DO IE = 0,NE1,1
         DO IP = 0,NP1,1
            WRITE (IOF,*) DREAL(GRNEP(IE,IP)), DIMAG(GRNEP(IE,IP))
         ENDDO
         WRITE (IOF,*) "--- ", IE, " ---"
      ENDDO
      CLOSE (IOF)
      RETURN
      END
C
      SUBROUTINE RUNGE_KUTTA (X,H,GG,GGN)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H
      DOUBLE COMPLEX GG(N),GGN(N)
      DOUBLE COMPLEX KK1(N),KK2(N),KK3(N),KK4(N)
      CALL K1FUNC (X,H,GG,    KK1)
      CALL K2FUNC (X,H,GG,KK1,KK2)
      CALL K3FUNC (X,H,GG,KK2,KK3)
      CALL K4FUNC (X,H,GG,KK3,KK4)
      DO I=1,N
         GGN(I) = GG(I) + 1D0/6D0 *
     -        (KK1(I) + 2D0*KK2(I) + 2D0*KK3(I) + KK4(I))
      ENDDO
      RETURN
      END
C
      SUBROUTINE RKFUNC (X,H,GG,FF)
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=2)
      DOUBLE COMPLEX CONE,CIMG
      PARAMETER (CONE=(1D0,0D0), CIMG=(0D0,1D0))
      DOUBLE PRECISION L10
      PARAMETER (L10=2.3025851D0)
      DOUBLE PRECISION X,H,R
      DOUBLE COMPLEX GG(N), FF(N)
      DOUBLE COMPLEX FFUNC
      EXTERNAL FFUNC
      R = 10D0**X
      FF(1) = H * GG(2)
      FF(2) = H * (L10*GG(2) - L10**2*R**2 * FFUNC (X,GG(1)))
      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION FFUNC (X,G1)
      IMPLICIT NONE
      DOUBLE COMPLEX CONE,CIMG
      PARAMETER (CONE=(1D0,0D0), CIMG=(0D0,1D0))
      DOUBLE PRECISION X,R
      DOUBLE COMPLEX G1
      DOUBLE PRECISION VQCD
      EXTERNAL VQCD
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      R = 10D0**X
      FFUNC = MU * (E + CIMG*GAMT - VQCD (R)) * G1
      RETURN
      END
C
      SUBROUTINE K1FUNC (X,H,GG,KK1)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X1
      DOUBLE COMPLEX GG(N),GG1(N),KK1(N)
      X1 = X
      DO I=1,N
         GG1(I) = GG(I)
      ENDDO
      CALL RKFUNC (X1,H,GG1,KK1)
      RETURN
      END
C
      SUBROUTINE K2FUNC (X,H,GG,KK1,KK2)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X2
      DOUBLE COMPLEX GG(N),GG2(N),KK1(N),KK2(N)
      X2 = X + H/2D0
      DO I=1,N
         GG2(I) = GG(I) + KK1(I)/2D0
      ENDDO
      CALL RKFUNC (X2,H,GG2,KK2)
      RETURN
      END
C
      SUBROUTINE K3FUNC (X,H,GG,KK2,KK3)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X3
      DOUBLE COMPLEX GG(N),GG3(N),KK2(N),KK3(N)
      X3 = X + H/2D0
      DO I=1,N
         GG3(I) = GG(I) + KK2(I)/2D0
      ENDDO
      CALL RKFUNC (X3,H,GG3,KK3)
      RETURN
      END
C
      SUBROUTINE K4FUNC (X,H,GG,KK3,KK4)
      IMPLICIT NONE
      INTEGER I,N
      PARAMETER (N=2)
      DOUBLE PRECISION X,H,X4
      DOUBLE COMPLEX GG(N),GG4(N),KK3(N),KK4(N)
      X4 = X + H
      DO I=1,N
         GG4(I) = GG(I) + KK3(I)
      ENDDO
      CALL RKFUNC (X4,H,GG4,KK4)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION S1 (P,R,MU)
      IMPLICIT NONE
      DOUBLE PRECISION P,R,MU
      IF (P .EQ. 0D0) THEN
         S1 = MU * R
      ELSEIF (P*R .LT. 1D-3) THEN
         S1 = MU * ( R - R**3*P**2/6D0 + R**5*P**4/120D0 )
      ELSE
         S1 = MU/P * DSIN(P*R)
      ENDIF
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION EE (IE)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IE
      EE = EMIN1 + DBLE(IE)/DBLE(NE1) * (EMAX1 - EMIN1)
      RETURN
      END
C
      INTEGER FUNCTION EINT (E)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION E
      EINT = DINT((E-EMIN1)/(EMAX1-EMIN1)*NE1)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PP (IP)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      INTEGER IP
      PP = PMIN1 + DBLE(IP)/DBLE(NP1) * (PMAX1 - PMIN1)
      RETURN
      END
C
      INTEGER FUNCTION PINT (P)
      IMPLICIT NONE
      INCLUDE 'tblgrnep.inc'
      DOUBLE PRECISION P
      PINT = DINT((P-PMIN1)/(PMAX1-PMIN1)*NP1)
      RETURN
      END
