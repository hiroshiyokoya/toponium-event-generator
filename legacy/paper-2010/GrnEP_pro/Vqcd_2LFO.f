C...  Two-Loop Fixed-Order QCD potential
C...  Position space representaion; M.Peter, NPB501,471(1997)
C...  A2 and A2O coefficient; B.Kniehl et.al., PLB607,96(2005)
C     and references threin.
      DOUBLE PRECISION FUNCTION VQCD (R)
      IMPLICIT NONE
      DOUBLE PRECISION R
      DOUBLE PRECISION MUB,ASB
      COMMON/GrnPOT/MUB,ASB
      DOUBLE PRECISION MU,GAMT,E
      COMMON/GrnPAR/MU,GAMT,E
      INTEGER NF
      PARAMETER (NF=5)
      DOUBLE PRECISION PI,SQ2,EGAM,ZET3,CF,CA,TR,C
      PARAMETER (PI=3.141592654D0,SQ2=1.41421356D0,EGAM=0.577216D0,
     -     ZET3=1.202057D0,CF=1.33333333D0,CA=3D0,TR=0.5D0)
      DOUBLE PRECISION A1,A2,A2O,B0,B1
      INTEGER ICLR
      COMMON/GrnCOL/ICLR
      DOUBLE PRECISION VH,GF,MH
      PARAMETER (GF=1.16637D-5,MH=120D0)
      B0 = 11D0/3D0*CA - 4D0/3D0*TR*NF
      B1 = 34D0/3D0*CA**2 - 20D0/3D0*CA*TR*NF - 4D0*CF*TR*NF
      A1 = 31D0/9D0*CA - 20D0/9D0*TR*NF
      A2 = (4343D0/162D0 + 4D0*PI**2 - PI**4/4D0 + 22D0/3D0*ZET3)*CA**2
     -     - (1798D0/81D0 + 56D0/3D0*ZET3)*CA*TR*NF
     -     - (55D0/3D0-16D0*ZET3)*CF*TR*NF
     -     + (20D0/9D0*TR*NF)**2
      A2O = (PI**4-12D0*PI**2)*CA**2
      IF (ICLR.EQ.1) THEN
         C = -CF                ! Singlet
      ELSEIF (ICLR.EQ.2) THEN
         C = 1D0/6D0            ! Octet
         A2 = A2 + A2O
      ELSE
         STOP
      ENDIF
      VQCD = C * ASB/R
     -     * ( 1D0 + ASB/(4D0*PI) * (2D0*B0*(DLOG(MUB*R) + EGAM) + A1)
     -     + (ASB/(4D0*PI))**2 * (B0**2*(4D0*(DLOG(MUB*R)+EGAM)**2
     -     + PI**2/3D0)+2D0*(B1+2D0*B0*A1)*(DLOG(MUB*R)+EGAM)+A2) )
*     VH = -GF*MU**2/2D0/SQ2/PI * DEXP(-MH*R) / R
*     VQCD = VQCD + VH
      RETURN
      END
