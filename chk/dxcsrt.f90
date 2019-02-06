!*==DXCSRT.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DXCSRT
      SUBROUTINE DXCSRT(Dnu1,Nudiff,Mu1,Mu2,Theta,P,Q,R,Ip,Iq,Ir,C1,Ic1,C2,Ic2,
     &                  Ierror)
      IMPLICIT NONE
!*--DXCSRT6
!*** Start of declarations inserted by SPAG
      INTEGER i , Ic1 , Ic2 , Ierror , Ip , Iq , Ir , ix1 , ix2 , j , k , l , 
     &        lm1 , mu , Mu1 , Mu2 , Nudiff
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DXCSRT
!***PURPOSE  TO COMPUTE CHECK VALUES FOR LEGENDRE FUNCTIONS
!***LIBRARY   SLATEC
!***CATEGORY  C3A2, C9
!***TYPE      SINGLE PRECISION (XCRST-S, DXCRST-D)
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  SMITH, JOHN M., (NBS AND GEORGE MASON UNIVERSITY)
!***DESCRIPTION
!
!        SUBROUTINE DXCSRT CALCULATES CASORATI (CROSS PRODUCT)
!        CHECK VALUES AND STORES THEM IN ARRAYS C1 AND C2 WITH
!        EXPONENTS IN ARRAYS IC1 AND IC2.  CALCULATIONS ARE BASED
!        ON PREVIOUSLY CALCULATED LEGENDRE FUNCTIONS OF THE
!        FIRST KIND (NEGATIVE ORDER) IN ARRAY P, THE SECOND KIND
!        IN ARRAY Q, THE FIRST KIND (POSITIVE ORDER) IN ARRAY R.
!        RESULTS SHOULD BE 1.0 TO WITHIN ROUNDOFF ERROR.
!
!***SEE ALSO  FCNQX2
!***REFERENCES  OLVER AND SMITH,J.COMPUT.PHYSICS,51(1983),NO.3,502-518.
!***ROUTINES CALLED  DXADD, DXADJ, DXRED
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!***END PROLOGUE  DXCSRT
      DOUBLE PRECISION C1 , C2 , dmu , dmu1 , nu , Dnu1 , P , Q , R , Theta , 
     &                 sx , x1 , x2
      DIMENSION P(*) , Ip(*) , Q(*) , Iq(*) , R(*) , Ir(*)
      DIMENSION C1(*) , Ic1(*) , C2(*) , Ic2(*)
!
!         PLACE ALL INPUT IN ADJUSTED FORM.
!
!***FIRST EXECUTABLE STATEMENT  DXCSRT
      Ierror = 0
      l = Nudiff + (Mu2-Mu1) + 1
      lm1 = l - 1
      DO i = 1 , l
        CALL DXADJ(P(i),Ip(i),Ierror)
        IF ( Ierror/=0 ) RETURN
        CALL DXADJ(Q(i),Iq(i),Ierror)
        IF ( Ierror/=0 ) RETURN
        CALL DXADJ(R(i),Ir(i),Ierror)
        IF ( Ierror/=0 ) RETURN
      ENDDO
!
!         CHECKS FOR FIXED MU, VARIABLE NU
!
      IF ( Mu2>Mu1 ) THEN
!
!         CHECKS FOR FIXED NU, VARIABLE MU
!
        sx = SIN(Theta)
        DO i = 1 , lm1
          C1(i) = 0.D0
          C2(i) = 0.D0
!
!         CASORATI 4
!
!         (MU+NU+1)*(MU-NU)*P(-(MU+1),NU,X)*Q(MU,NU,X)
!              -P(-MU,NU,X)*Q(MU+1,NU,X)=COS(MU*PI)/SQRT(1-X**2)
!
          mu = Mu1 + i - 1
          dmu = mu
          x1 = P(i+1)*Q(i)
          ix1 = Ip(i+1) + Iq(i)
          CALL DXADJ(x1,ix1,Ierror)
          IF ( Ierror/=0 ) RETURN
          x2 = P(i)*Q(i+1)
          ix2 = Ip(i) + Iq(i+1)
          CALL DXADJ(x2,ix2,Ierror)
          IF ( Ierror/=0 ) RETURN
          x1 = (dmu+Dnu1+1.D0)*(dmu-Dnu1)*x1
!
!         MULTIPLY BY SQRT(1-X**2)*(-1)**MU SO THAT CHECK VALUE IS 1.
!
          CALL DXADD(x1,ix1,-x2,ix2,C1(i),Ic1(i),Ierror)
          IF ( Ierror/=0 ) RETURN
          C1(i) = sx*C1(i)*(-1)**mu
          CALL DXADJ(C1(i),Ic1(i),Ierror)
          IF ( Ierror/=0 ) RETURN
!
!         CASORATI 3
!
!         P(MU+1,NU,X)*Q(MU,NU,X)-P(MU,NU,X)*Q(MU+1,NU,X)=
!               GAMMA(NU+MU+1)/(GAMMA(NU-MU+1)*SQRT(1-X**2))
!
          IF ( dmu<Dnu1+1.D0.OR.MOD(Dnu1,1.D0)/=0.D0 ) THEN
            x1 = R(i+1)*Q(i)
            ix1 = Ir(i+1) + Iq(i)
            CALL DXADJ(x1,ix1,Ierror)
            IF ( Ierror/=0 ) RETURN
            x2 = R(i)*Q(i+1)
            ix2 = Ir(i) + Iq(i+1)
            CALL DXADJ(x2,ix2,Ierror)
            IF ( Ierror/=0 ) RETURN
            CALL DXADD(x1,ix1,-x2,ix2,C2(i),Ic2(i),Ierror)
            IF ( Ierror/=0 ) RETURN
!
!         MULTIPLY BY SQRT(1-X**2) AND THEN DIVIDE BY
!         (NU+MU),(NU+MU-1),(NU+MU-2),...,(NU-MU+1) SO THAT
!         CHECK VALUE IS 1.
!
            C2(i) = C2(i)*sx
            k = 2*mu
            IF ( k>0 ) THEN
              DO j = 1 , k
                C2(i) = C2(i)/(Dnu1+dmu+1.D0-j)
                CALL DXADJ(C2(i),Ic2(i),Ierror)
              ENDDO
              IF ( Ierror/=0 ) RETURN
            ENDIF
          ENDIF
        ENDDO
      ELSE
        dmu1 = Mu1
        DO i = 1 , lm1
          C1(i) = 0.D0
          C2(i) = 0.D0
          nu = Dnu1 + i - 1.D0
!
!         CASORATI 2
!
!         (MU+NU+1)*P(-MU,NU+1,X)*Q(MU,NU,X)
!               +(MU-NU-1)*P(-MU,NU,X)*Q(MU,NU+1,X)=COS(MU*PI)
!
          x1 = P(i+1)*Q(i)
          ix1 = Ip(i+1) + Iq(i)
          CALL DXADJ(x1,ix1,Ierror)
          IF ( Ierror/=0 ) RETURN
          x2 = P(i)*Q(i+1)
          ix2 = Ip(i) + Iq(i+1)
          CALL DXADJ(x2,ix2,Ierror)
          IF ( Ierror/=0 ) RETURN
          x1 = (dmu1+nu+1.D0)*x1
          x2 = (dmu1-nu-1.D0)*x2
          CALL DXADD(x1,ix1,x2,ix2,C1(i),Ic1(i),Ierror)
          IF ( Ierror/=0 ) RETURN
          CALL DXADJ(C1(i),Ic1(i),Ierror)
          IF ( Ierror/=0 ) RETURN
!
!         MULTIPLY BY (-1)**MU SO THAT CHECK VALUE IS 1.
!
          C1(i) = C1(i)*(-1)**Mu1
!
!         CASORATI 1
!
!         P(MU,NU+1,X)*Q(MU,NU,X)-P(MU,NU,X)*Q(MU,NU+1,X)=
!               GAMMA(NU+MU+1)/GAMMA(NU-MU+2)
!
          IF ( dmu1>=nu+1.D0.AND.MOD(nu,1.D0)==0.D0 ) THEN
            C2(i) = 0.D0
            Ic2(i) = 0
          ELSE
            x1 = R(i+1)*Q(i)
            ix1 = Ir(i+1) + Iq(i)
            CALL DXADJ(x1,ix1,Ierror)
            IF ( Ierror/=0 ) RETURN
            x2 = R(i)*Q(i+1)
            ix2 = Ir(i) + Iq(i+1)
            CALL DXADJ(x2,ix2,Ierror)
            IF ( Ierror/=0 ) RETURN
            CALL DXADD(x1,ix1,-x2,ix2,C2(i),Ic2(i),Ierror)
            IF ( Ierror/=0 ) RETURN
!
!         DIVIDE BY (NU+MU),(NU+MU-1),(NU+MU-2),....(NU-MU+2),
!         SO THAT CHECK VALUE IS 1.
!
            k = 2*Mu1 - 1
            DO j = 1 , k
              IF ( k>0 ) C2(i) = C2(i)/(nu+dmu1+1.D0-j)
              CALL DXADJ(C2(i),Ic2(i),Ierror)
            ENDDO
            IF ( Ierror/=0 ) RETURN
            IF ( k<=0 ) C2(i) = (nu+1.D0)*C2(i)
          ENDIF
        ENDDO
      ENDIF
!
!         PLACE RESULTS IN REDUCED FORM.
!
      DO i = 1 , lm1
        CALL DXRED(C1(i),Ic1(i),Ierror)
        IF ( Ierror/=0 ) RETURN
        CALL DXRED(C2(i),Ic2(i),Ierror)
        IF ( Ierror/=0 ) RETURN
      ENDDO
      END SUBROUTINE DXCSRT
