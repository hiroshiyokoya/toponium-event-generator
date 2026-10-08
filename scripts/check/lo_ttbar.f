C=======================================================================
C     Reference LO cross sections for on-shell pp -> t tbar, used to check
C     the channel decomposition of the generators (issue #9).
C
C     Partonic cross sections (Combridge 1979; Glueck, Owens, Reya 1978):
C       qq~ -> tt~ : 8 pi as^2/(27 s) beta (1 + rho/2)
C       gg  -> tt~ : pi as^2/(3 s) [ (1 + rho + rho^2/16) ln((1+b)/(1-b))
C                                    - b (7/4 + 31 rho/16) ]
C       rho = 4 m^2/s, beta = sqrt(1 - rho)
C     convoluted with LHAPDF (LHAGlue) cteq6l1, alpha_s from the set,
C     mu_R = mu_F = m_t, 4 light flavours (d,u,s,c) in qq~ as in the
C     generators.
C
C     usage: lo_ttbar [sqrt(s) GeV] [m_t GeV]     (defaults 14000, 173)
C=======================================================================
      PROGRAM LOTT
      IMPLICIT NONE
      DOUBLE PRECISION RS,MT,S,MU,AS,PI,GEVPB
      PARAMETER (PI=3.14159265358979D0, GEVPB=0.389379D9)
      DOUBLE PRECISION TAU0,TAU,Y,X1,X2,SH,RHO,BETA,SQQ,SGG
      DOUBLE PRECISION F1(-6:6),F2(-6:6),LQQ,LGG
      DOUBLE PRECISION SIGQQ,SIGGG,WT,U,V,YMAX,ALPHASPDF
      EXTERNAL ALPHASPDF
      INTEGER NU,NV,I,J,K
      PARAMETER (NU=2000, NV=400)
      CHARACTER*32 ARG
      RS = 14000D0
      MT = 173D0
      IF (COMMAND_ARGUMENT_COUNT().GE.1) THEN
         CALL GET_COMMAND_ARGUMENT(1,ARG)
         READ (ARG,*) RS
      ENDIF
      IF (COMMAND_ARGUMENT_COUNT().GE.2) THEN
         CALL GET_COMMAND_ARGUMENT(2,ARG)
         READ (ARG,*) MT
      ENDIF
      CALL INITPDFSETBYNAME ('cteq6l1')
      CALL INITPDF (0)
      S = RS**2
      MU = MT
      AS = ALPHASPDF (MU)
      TAU0 = 4D0*MT**2/S
      SIGQQ = 0D0
      SIGGG = 0D0
C     midpoint rule in u = ln(tau) and v = rapidity fraction
      DO I = 1, NU
         U = DLOG(TAU0) * (1D0 - (DBLE(I)-0.5D0)/DBLE(NU))
         TAU = DEXP(U)
         SH = TAU*S
         RHO = 4D0*MT**2/SH
         BETA = DSQRT(1D0-RHO)
         SQQ = 8D0*PI*AS**2/(27D0*SH) * BETA*(1D0+RHO/2D0)
         SGG = PI*AS**2/(3D0*SH) * ((1D0+RHO+RHO**2/16D0)
     .        * DLOG((1D0+BETA)/(1D0-BETA))
     .        - BETA*(7D0/4D0+31D0*RHO/16D0))
         YMAX = -0.5D0*U
         DO J = 1, NV
            V = (DBLE(J)-0.5D0)/DBLE(NV)
            Y = YMAX*(2D0*V-1D0)
            X1 = DSQRT(TAU)*DEXP(Y)
            X2 = DSQRT(TAU)*DEXP(-Y)
            CALL EVOLVEPDF (X1,MU,F1)
            CALL EVOLVEPDF (X2,MU,F2)
            LQQ = 0D0
            DO K = 1, 4
               LQQ = LQQ + F1(K)*F2(-K) + F1(-K)*F2(K)
            ENDDO
            LGG = F1(0)*F2(0)
C           dsigma = dtau dy f(x1) f(x2) sigma^, with f = xf/x
C           dtau = tau du, dy = 2 ymax dv, 1/(x1 x2) = 1/tau
            WT = (-DLOG(TAU0)/DBLE(NU)) * (2D0*YMAX/DBLE(NV))
            SIGQQ = SIGQQ + WT*LQQ*SQQ
            SIGGG = SIGGG + WT*LGG*SGG
         ENDDO
      ENDDO
      SIGQQ = SIGQQ*GEVPB
      SIGGG = SIGGG*GEVPB
      WRITE (6,'(A,F9.1,A,F7.2,A,F8.5)') ' sqrt(s) =',RS,'  m_t =',MT,
     .     '  alpha_s(m_t) =',AS
      WRITE (6,'(A,F10.3,A)') ' sigma(qq~ -> tt~) =',SIGQQ,' pb'
      WRITE (6,'(A,F10.3,A)') ' sigma(gg  -> tt~) =',SIGGG,' pb'
      WRITE (6,'(A,F8.4)')    ' qq~/gg            =',SIGQQ/SIGGG
      END
