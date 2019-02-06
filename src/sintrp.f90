!*==SINTRP.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SINTRP
      SUBROUTINE SINTRP(X,Y,Xout,Yout,Ypout,Neqn,Kold,Phi,Ivc,Iv,Kgi,Gi,Alpha,
     &                  Og,Ow,Ox,Oy)
      IMPLICIT NONE
!*--SINTRP6
!*** Start of declarations inserted by SPAG
      REAL alp , Alpha , c , g , gamma , gdi , gdif , Gi , h , hi , hmu , Og , 
     &     Ow , Ox , Oy , Phi , rmu , sigma , temp1 , temp2
      REAL temp3 , w , X , xi , xim1 , xiq , Xout , Y , Yout , Ypout
      INTEGER i , iq , Iv , Ivc , iw , j , jq , Kgi , Kold , kp1 , kp2 , l , m , 
     &        Neqn
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SINTRP
!***PURPOSE  Approximate the solution at XOUT by evaluating the
!            polynomial computed in STEPS at XOUT.  Must be used in
!            conjunction with STEPS.
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      SINGLE PRECISION (SINTRP-S, DINTP-D)
!***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
!             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
!             SMOOTH INTERPOLANT
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   The methods in subroutine  STEPS  approximate the solution near  X
!   by a polynomial.  Subroutine  SINTRP  approximates the solution at
!   XOUT  by evaluating the polynomial there.  Information defining this
!   polynomial is passed from  STEPS  so  SINTRP  cannot be used alone.
!
!   Subroutine STEPS is completely explained and documented in the text,
!   "Computer Solution of Ordinary Differential Equations, the Initial
!   Value Problem"  by L. F. Shampine and M. K. Gordon.
!
!   Input to SINTRP --
!
!   The user provides storage in the calling program for the arrays in
!   the call list
!      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
!                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
!   and defines
!      XOUT -- point at which solution is desired.
!   The remaining parameters are defined in  STEPS  and passed to
!   SINTRP  from that subroutine
!
!   Output from  SINTRP --
!
!      YOUT(*) -- solution at  XOUT
!      YPOUT(*) -- derivative of solution at  XOUT
!   The remaining parameters are returned unaltered from their input
!   values.  Integration with  STEPS  may be continued.
!
!***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
!                 II, Report SAND84-0293, Sandia Laboratories, 1984.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   840201  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINTRP
!
      DIMENSION Y(*) , Yout(*) , Ypout(*) , Phi(Neqn,16) , Oy(*)
      DIMENSION g(13) , c(13) , w(13) , Og(13) , Ow(12) , Alpha(12) , Gi(11) , 
     &          Iv(10)
!
!***FIRST EXECUTABLE STATEMENT  SINTRP
      kp1 = Kold + 1
      kp2 = Kold + 2
!
      hi = Xout - Ox
      h = X - Ox
      xi = hi/h
      xim1 = xi - 1.
!
!   INITIALIZE W(*) FOR COMPUTING G(*)
!
      xiq = xi
      DO iq = 1 , kp1
        xiq = xi*xiq
        temp1 = iq*(iq+1)
        w(iq) = xiq/temp1
      ENDDO
!
!   COMPUTE THE DOUBLE INTEGRAL TERM GDI
!
      IF ( Kold<=Kgi ) THEN
        gdi = Gi(Kold)
      ELSE
        IF ( Ivc>0 ) THEN
          iw = Iv(Ivc)
          gdi = Ow(iw)
          m = Kold - iw + 3
        ELSE
          gdi = 1.0/temp1
          m = 2
        ENDIF
        IF ( m<=Kold ) THEN
          DO i = m , Kold
            gdi = Ow(kp2-i) - Alpha(i)*gdi
          ENDDO
        ENDIF
      ENDIF
!
!   COMPUTE G(*) AND C(*)
!
      g(1) = xi
      g(2) = 0.5*xi*xi
      c(1) = 1.0
      c(2) = xi
      IF ( Kold>=2 ) THEN
        DO i = 2 , Kold
          alp = Alpha(i)
          gamma = 1.0 + xim1*alp
          l = kp2 - i
          DO jq = 1 , l
            w(jq) = gamma*w(jq) - alp*w(jq+1)
          ENDDO
          g(i+1) = w(1)
          c(i+1) = gamma*c(i)
        ENDDO
      ENDIF
!
!   DEFINE INTERPOLATION PARAMETERS
!
      sigma = (w(2)-xim1*w(1))/gdi
      rmu = xim1*c(kp1)/gdi
      hmu = rmu/h
!
!   INTERPOLATE FOR THE SOLUTION -- YOUT
!   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
!
      DO l = 1 , Neqn
        Yout(l) = 0.0
        Ypout(l) = 0.0
      ENDDO
      DO j = 1 , Kold
        i = kp2 - j
        gdif = Og(i) - Og(i-1)
        temp2 = (g(i)-g(i-1)) - sigma*gdif
        temp3 = (c(i)-c(i-1)) + rmu*gdif
        DO l = 1 , Neqn
          Yout(l) = Yout(l) + temp2*Phi(l,i)
          Ypout(l) = Ypout(l) + temp3*Phi(l,i)
        ENDDO
      ENDDO
      DO l = 1 , Neqn
        Yout(l) = ((1.0-sigma)*Oy(l)+sigma*Y(l))
     &            + h*(Yout(l)+(g(1)-sigma*Og(1))*Phi(l,1))
        Ypout(l) = hmu*(Oy(l)-Y(l)) + (Ypout(l)+(c(1)+rmu*Og(1))*Phi(l,1))
      ENDDO
!
      END SUBROUTINE SINTRP
