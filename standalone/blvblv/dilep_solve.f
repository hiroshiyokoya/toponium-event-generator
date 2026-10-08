      SUBROUTINE DILEPSOL (PBT,PLP,PBX,PLM,PXMISS,PYMISS,MT,MW,MB)
      IMPLICIT NONE
      INTEGER I,J
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592654D0)
      DOUBLE PRECISION MT,MW,MB
      DOUBLE PRECISION MT1,MT2,MW1,MW2,MB1,MB2,ML1,ML2
      DOUBLE PRECISION PBT(3),PBX(3),PLP(3),PLM(3)
      DOUBLE PRECISION PXMISS,PYMISS
      DOUBLE PRECISION PN1(3),PN2(3)
      DOUBLE PRECISION SOL(16,6)
      INTEGER NSOL
      COMMON/SOL/SOL,NSOL
      DOUBLE PRECISION PB1(3),PB2(3),PL1(3),PL2(3)
      DOUBLE PRECISION EX,EY
      DOUBLE PRECISION EB1,EB2,EL1,EL2
      DOUBLE PRECISION A1,A2,A3,A4,B1,B2,B3,B4
      DOUBLE PRECISION C22,C21,C20,C11,C10,C00
      DOUBLE PRECISION D22,D21,D20,D11,D10,D00
      DOUBLE PRECISION E22,E21,E20,E11,E10,E00
      DOUBLE PRECISION H0,H1,H2,H3,H4
      DOUBLE PRECISION K1,K2,K3
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION S1,S2,S3,Q,R,THETA
      DOUBLE PRECISION U,V
      DOUBLE PRECISION Z(3),POS(3),X(12),Y(8)
      INTEGER IZ,IPOS,IROT,ISOL
      DOUBLE PRECISION C0,C1,C2,D0,D1,D2
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-3)
      DOUBLE PRECISION ENERGY,DOTPRD
      EXTERNAL ENERGY,DOTPRD
C.... SCALE THE MASS DIMENSION VARIABLES
C.... OBSERVED MOMENTUM
      DO I = 1, 3
         PB1(I) = PBT(I)/MT
         PB2(I) = PBX(I)/MT
         PL1(I) = PLP(I)/MT
         PL2(I) = PLM(I)/MT
      ENDDO
C.... MISSING TRANSVERSE MOMENTUM
      EX = PXMISS/MT
      EY = PYMISS/MT
C.... MASS
      MT1 = MT/MT
      MT2 = MT/MT
      MW1 = MW/MT
      MW2 = MW/MT
      MB1 = MB/MT
      MB2 = MB/MT
      ML1 = 0D0
      ML2 = 0D0
C.... ENERGY
      EB1 = ENERGY (PB1,MB1)
      EB2 = ENERGY (PB2,MB2)
      EL1 = ENERGY (PL1,ML1)
      EL2 = ENERGY (PL2,ML2)
*     WRITE (6,*) "Ene", EB1,EB2,EL1,EL2
C.....
      A1 = (EB1 + EL1)*(MW1**2 - ML1**2)
     -     - EL1*(MT1**2 - MB1**2 - ML1**2)
     -     + 2D0*EB1*EL1**2 - 2D0*EL1*DOTPRD(PB1,PL1)
      A2 = 2D0*(EB1*PL1(1) - EL1*PB1(1))
      A3 = 2D0*(EB1*PL1(2) - EL1*PB1(2))
      A4 = 2D0*(EB1*PL1(3) - EL1*PB1(3))
*     WRITE (6,*) "A",A1,A2,A3,A4
C.....
      B1 = (EB2 + EL2)*(MW2**2 - ML2**2)
     -     - EL2*(MT2**2 - MB2**2 - ML2**2)
     -     + 2D0*EB2*EL2**2 - 2D0*EL2*DOTPRD(PB2,PL2)
      B2 = 2D0*(EB2*PL2(1) - EL2*PB2(1))
      B3 = 2D0*(EB2*PL2(2) - EL2*PB2(2))
      B4 = 2D0*(EB2*PL2(3) - EL2*PB2(3))
*     WRITE (6,*) "B",B1,B2,B3,B4
C.....
      C22 = (MW1**2-ML1**2)**2 - 4D0*(EL1**2-PL1(3)**2)*(A1/A4)**2
     -     - 4D0*(MW1**2-ML1**2)*PL1(3)*A1/A4
      C21 = 4D0*(MW1**2-ML1**2)*(PL1(1)-PL1(3)*A2/A4) -
     -     8D0*(EL1**2-PL1(3)**2)*A1*A2/A4**2 - 8D0*PL1(1)*PL1(3)*A1/A4
      C20 = - 4D0*(EL1**2-PL1(1)**2) - 4D0*(EL1**2-PL1(3)**2)*(A2/A4)**2
     -     - 8D0*PL1(1)*PL1(3)*A2/A4
      C11 = 4D0*(MW1**2-ML1**2)*(PL1(2)-PL1(3)*A3/A4) -
     -     8D0*(EL1**2-PL1(3)**2)*A1*A3/A4**2 - 8D0*PL1(2)*PL1(3)*A1/A4
      C10 = - 8D0*(EL1**2-PL1(3)**2)*A2*A3/A4**2 + 8D0*PL1(1)*PL1(2) -
     -     8D0*PL1(1)*PL1(3)*A3/A4 - 8D0*PL1(2)*PL1(3)*A2/A4
      C00 = - 4D0*(EL1**2-PL1(2)**2) - 4D0*(EL1**2-PL1(3)**2)*(A3/A4)**2
     -     - 8D0*PL1(2)*PL1(3)*A3/A4
*     WRITE (6,*) "C",C22,C21,C20,C11,C10,C00
C.....
      E22 = (MW2**2-ML2**2)**2 - 4D0*(EL2**2-PL2(3)**2)*(B1/B4)**2
     -     - 4D0*(MW2**2-ML2**2)*PL2(3)*B1/B4
      E21 = 4D0*(MW2**2-ML2**2)*(PL2(1)-PL2(3)*B2/B4) -
     -     8D0*(EL2**2-PL2(3)**2)*B1*B2/B4**2 - 8D0*PL2(1)*PL2(3)*B1/B4
      E20 = - 4D0*(EL2**2-PL2(1)**2) - 4D0*(EL2**2-PL2(3)**2)*(B2/B4)**2
     -     - 8D0*PL2(1)*PL2(3)*B2/B4
      E11 = 4D0*(MW2**2-ML2**2)*(PL2(2)-PL2(3)*B3/B4) -
     -     8D0*(EL2**2-PL2(3)**2)*B1*B3/B4**2 - 8D0*PL2(2)*PL2(3)*B1/B4
      E10 = - 8D0*(EL2**2-PL2(3)**2)*B2*B3/B4**2 + 8D0*PL2(1)*PL2(2) -
     -     8D0*PL2(1)*PL2(3)*B3/B4 - 8D0*PL2(2)*PL2(3)*B2/B4
      E00 = - 4D0*(EL2**2-PL2(2)**2) - 4D0*(EL2**2-PL2(3)**2)*(B3/B4)**2
     -     - 8D0*PL2(2)*PL2(3)*B3/B4
*     WRITE (6,*) "E",E22,E21,E20,E11,E10,E00
C.....
      D22 = E22 + EX**2*E20 + EY**2*E00 + EX*EY*E10 + EX*E21 + EY*E11
      D21 = - E21 - 2D0*EX*E20 - EY*E10
      D20 = E20
      D11 = - E11 - 2D0*EY*E00 - EX*E10
      D10 = E10
      D00 = E00
*     WRITE (6,*) "D",D22,D21,D20,D11,D10,D00
C.....
      H4 = C00**2*D22**2 + C11*D22*(C11*D00-C00*D11)
     -     + C00*C22*(D11**2-2D0*D00*D22) + C22*D00*(C22*D00-C11*D11)
      H3 = C00*D21*(2D0*C00*D22-C11*D11) + C00*D11*(2D0*C22*D10+C21*D11)
     -     + C22*D00*(2D0*C21*D00-C11*D10) - C00*D22*(C11*D10+C10*D11)
     -     - 2D0*C00*D00*(C22*D21+C21*D22) - D00*D11*(C11*C21+C10*C22)
     -     + C11*D00*(C11*D21+2D0*C10*D22)
      H2 = C00**2*(2D0*D22*D20+D21**2) - C00*D21*(C11*D10+C10*D11)
     -     + C11*D20*(C11*D00-C00*D11) + C00*D10*(C22*D10-C10*D22)
     -     + C00*D11*(2D0*C21*D10+C20*D11) + (2D0*C22*C20+C21**2)*D00**2
     -     - 2D0*C00*D00*(C22*D20+C21*D21+C20*D22)
     -     + C10*D00*(2D0*C11*D21+C10*D22) - D00*D10*(C11*C21+C10*C22)
     -     - D00*D11*(C11*C20+C10*C21)
*      H1 = C10**2*D00*D21 + C20*(2*C21*D00**2-C11*D00*D10-C10*D00*D11 +
*     -     2*C00*D10*D11 - 2*C00*D00*D21) -
*     -  C10*(C21*D00*D10 - 2*C11*D00*D20 + C00*D11*D20 + C00*D10*D21) +
*     -  C00*(C21*D10**2 - 2*C21*D00*D20 - C11*D10*D20 + 2*C00*D20*D21)
      H1 = C00*D21*(2D0*C00*D20-C10*D10) - C00*D20*(C11*D10+C10*D11)
     -     + C00*D10*(C21*D10+2D0*C20*D11)
     -     - 2D0*C00*D00*(C21*D20+C20*D21)
     -     + C10*D00*(2D0*C11*D20+C10*D21)
     -     + C20*D00*(2D0*C21*D00-C10*D11) - D00*D10*(C11*C20+C10*C21)
      H0 = C00**2*D20**2 + C10*D20*(C10*D00-C00*D10)
     -     + C20*D10*(C00*D10-C10*D00) + C20*D00*(C20*D00-2D0*C00*D20)
*     WRITE (6,*) "H",H4,H3,H2,H1,H0
C.....
      H1 = H1/H0
      H2 = H2/H0
      H3 = H3/H0
      H4 = H4/H0
C.....
      K1 = H2 - 3D0*H1**2/8D0
      K2 = H3 + H1**3/8D0 - H1*H2/2D0
      K3 = H4 - 3D0*H1**4/256D0 + H1**2*H2/16D0 - H1*H3/4D0
*     WRITE (6,*) "K",K1,K2,K3
C.....
      IZ = 0
      IPOS = 0
      IROT = 0
      ISOL = 0
C.....
      S1 = 2D0*K1
      S2 = K1**2 - 4D0*K3
      S3 = - K2**2
      Q = (S1**2 - 3D0*S2) / 9D0
      R = (2D0*S1**3 - 9D0*S1*S2 + 27D0*S3)/54D0
      IF (R**2.LT.Q**3) THEN
         IZ = 3
         THETA = DACOS(R/DSQRT(Q**3))
         Z(1) = -2D0*DSQRT(Q)*DCOS(THETA/3D0) - S1/3D0
         Z(2) = -2D0*DSQRT(Q)*DCOS((THETA+2D0*PI)/3D0) - S1/3D0
         Z(3) = -2D0*DSQRT(Q)*DCOS((THETA-2D0*PI)/3D0) - S1/3D0
      ELSE
         IZ = 1
         U = -R/DABS(R)*(DABS(R)+DSQRT(R**2-Q**3))**(1D0/3D0)
         V = Q/U
*     U = (- R + DSQRT(R**2-Q**3))**(1./3.)
*     V = (- R - DSQRT(R**2-Q**3))**(1./3.)
         Z(1) = U + V - S1/3D0
      ENDIF
*     WRITE (6,*) "IZ", IZ, (Z(I),I=1,IZ)
      IF (IZ.EQ.0) RETURN
      DO I = 1,IZ
*     WRITE (6,*) "CUBIC SOL", I, Z(I)**3 + S1*Z(I)**2 + S2*Z(I) + S3
         IF (Z(I).GT.0D0) THEN
            IPOS = IPOS + 1
            POS(IPOS) = Z(I)
         ENDIF
      ENDDO
      IF (IPOS.EQ.0) RETURN
C.....
      DO I = 1, IPOS
         T1 = DSQRT(POS(I))
         T2 = (K1 + T1**2 - K2/T1) * 0.5D0
         IF (T1**2-4D0*T2.GT.0D0) THEN
            Q = -0.5D0 * (T1+T1/DABS(T1)*DSQRT(T1**2-4D0*T2))
            IROT = IROT + 1
            X(IROT) = Q
            IROT = IROT + 1
            X(IROT) = T2/Q
         ENDIF
         IF (T1**2-4D0*K3/T2.GT.0D0) THEN
            Q = -0.5D0 * (-T1-T1/DABS(T1)*DSQRT(T1**2-4D0*K3/T2))
            IROT = IROT + 1
            X(IROT) = Q
            IROT = IROT + 1
            X(IROT) = K3/T2/Q
         ENDIF
         T1 = -DSQRT(POS(I))
         T2 = (K1 + T1**2 - K2/T1) * 0.5D0
         IF (T1**2-4D0*T2.GT.0D0) THEN
            Q = -0.5D0 * (T1+T1/DABS(T1)*DSQRT(T1**2-4D0*T2))
            IROT = IROT + 1
            X(IROT) = Q
            IROT = IROT + 1
            X(IROT) = T2/Q
         ENDIF
         IF (T1**2-4D0*K3/T2.GT.0D0) THEN
            Q = -0.5D0 * (-T1-T1/DABS(T1)*DSQRT(T1**2-4D0*K3/T2))
            IROT = IROT + 1
            X(IROT) = Q
            IROT = IROT + 1
            X(IROT) = K3/T2/Q
         ENDIF
      ENDDO
      IF (IROT.EQ.0) RETURN     ! NO SOLUTION
*      DO I = 1, IROT
*         WRITE (6,*) "QUARTIC SOL", I, X(I)**4+K1*X(I)**2+K2*X(I)+K3
*      ENDDO
*      WRITE (6,*) IROT,(X(I),I=1,IROT)
C.....
      ISOL = 1
      Y(1) = X(1)
      DO I = 2, IROT
         DO J = 1, ISOL
            IF (DABS(X(I)-Y(J))/DABS(X(I)).LT.EPS) THEN
               GOTO 99
            ENDIF
         ENDDO
         ISOL = ISOL + 1
         Y(ISOL) = X(I)
 99      CONTINUE
      ENDDO
C.....
      DO I = NSOL+1, NSOL+ISOL
         PN1(1) = Y(I) - H1/4D0 ! PX1
         PN2(1) = EX - PN1(1)   ! PX2
         C0 = C00
         C1 = C11 + C10*PN1(1)
         C2 = C22 + C21*PN1(1) + C20*PN1(1)**2
         D0 = D00
         D1 = D11 + D10*PN1(1)
         D2 = D22 + D21*PN1(1) + D20*PN1(1)**2
         PN1(2) = (C0*D2-C2*D0)/(C1*D0-C0*D1) ! PY1
         PN2(2) = EY - PN1(2)   ! PY2
         PN1(3) = - (A1+A2*PN1(1)+A3*PN1(2))/A4 ! PZ1
         PN2(3) = - (B1+B2*PN2(1)+B3*PN2(2))/B4 ! PZ2
C.....
         SOL(I,1) = PN1(1) * MT
         SOL(I,2) = PN1(2) * MT
         SOL(I,3) = PN1(3) * MT
         SOL(I,4) = PN2(1) * MT
         SOL(I,5) = PN2(2) * MT
         SOL(I,6) = PN2(3) * MT
C.....
      ENDDO
C.....
      NSOL = NSOL + ISOL
C.....
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION ENERGY (P,M)
      IMPLICIT NONE
      DOUBLE PRECISION P(3),M
      ENERGY = DSQRT(P(1)**2+P(2)**2+P(3)**2+M**2)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DOTPRD (P1,P2)
      IMPLICIT NONE
      DOUBLE PRECISION P1(3),P2(3)
      DOTPRD = P1(1)*P2(1) + P1(2)*P2(2) + P1(3)*P2(3)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION MASCAL (P)
      IMPLICIT NONE
      DOUBLE PRECISION P(0:3)
      MASCAL = DSQRT(P(0)**2-P(1)**2-P(2)**2-P(3)**2)
      RETURN
      END
