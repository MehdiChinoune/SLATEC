!*==SDATRP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK SDATRP
SUBROUTINE SDATRP(X,Xout,Yout,Ypout,Neq,Kold,Phi,Psi)
  IMPLICIT NONE
  !*--SDATRP5
  !***BEGIN PROLOGUE  SDATRP
  !***SUBSIDIARY
  !***PURPOSE  Interpolation routine for SDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDATRP-S, DDATRP-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------------
  !     THE METHODS IN SUBROUTINE SDASTP USE POLYNOMIALS
  !     TO APPROXIMATE THE SOLUTION. SDATRP APPROXIMATES THE
  !     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING
  !     ONE OF THESE POLYNOMIALS, AND ITS DERIVATIVE,THERE.
  !     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM
  !     SDASTP, SO SDATRP CANNOT BE USED ALONE.
  !
  !     THE PARAMETERS ARE:
  !     X     THE CURRENT TIME IN THE INTEGRATION.
  !     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED
  !     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT
  !           (THIS IS OUTPUT)
  !     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
  !           (THIS IS OUTPUT)
  !     NEQ   NUMBER OF EQUATIONS
  !     KOLD  ORDER USED ON LAST SUCCESSFUL STEP
  !     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y
  !     PSI   ARRAY OF PAST STEPSIZE HISTORY
  !-----------------------------------------------------------------------
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !***END PROLOGUE  SDATRP
  !
  INTEGER Neq, Kold
  REAL X, Xout, Yout(*), Ypout(*), Phi(Neq,*), Psi(*)
  !
  INTEGER i, j, koldp1
  REAL c, d, gamma, temp1
  !
  !***FIRST EXECUTABLE STATEMENT  SDATRP
  koldp1 = Kold + 1
  temp1 = Xout - X
  DO i = 1, Neq
    Yout(i) = Phi(i,1)
    Ypout(i) = 0.0E0
  ENDDO
  c = 1.0E0
  d = 0.0E0
  gamma = temp1/Psi(1)
  DO j = 2, koldp1
    d = d*gamma + c/Psi(j-1)
    c = c*gamma
    gamma = (temp1+Psi(j-1))/Psi(j)
    DO i = 1, Neq
      Yout(i) = Yout(i) + c*Phi(i,j)
      Ypout(i) = Ypout(i) + d*Phi(i,j)
    ENDDO
  ENDDO
  !
  !------END OF SUBROUTINE SDATRP------
END SUBROUTINE SDATRP
