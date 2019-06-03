!** C9LN2R
COMPLEX(SP) FUNCTION C9LN2R(Z)
  !>
  !  Evaluate LOG(1+Z) from second order relative accuracy so
  !            that  LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      COMPLEX (R9LN2R-S, D9LN2R-D, C9LN2R-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  LOG(1+Z)  from 2-nd order with relative error accuracy so
  ! that     LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
  !
  ! Now  LOG(1+Z) = 0.5*LOG(1+2*X+ABS(Z)**2) + I*CARG(1+Z),
  ! where X = REAL(Z)  and  Y = AIMAG(Z).
  ! We find
  !     Z**3 * C9LN2R(Z) = -X*ABS(Z)**2 - 0.25*ABS(Z)**4
  !        + (2*X+ABS(Z)**2)**3 * R9LN2R(2*X+ABS(Z)**2)
  !        + I * (CARG(1+Z) + (X-1)*Y)
  ! The imaginary part must be evaluated carefully as
  !     (ATAN(Y/(1+X)) - Y/(1+X)) + Y/(1+X) - (1-X)*Y
  !       = (Y/(1+X))**3 * R9ATN1(Y/(1+X)) + X**2*Y/(1+X)
  !
  ! Now we divide through by Z**3 carefully.  Write
  !     1/Z**3 = (X-I*Y)/ABS(Z)**3 * (1/ABS(Z)**3)
  ! then   C9LN2R(Z) = ((X-I*Y)/ABS(Z))**3 * (-X/ABS(Z) - ABS(Z)/4
  !        + 0.5*((2*X+ABS(Z)**2)/ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
  !        + I*Y/(ABS(Z)*(1+X)) * ((X/ABS(Z))**2 +
  !          + (Y/(ABS(Z)*(1+X)))**2 * R9ATN1(Y/(1+X)) ) )
  !
  ! If we let  XZ = X/ABS(Z)  and  YZ = Y/ABS(Z)  we may write
  !     C9LN2R(Z) = (XZ-I*YZ)**3 * (-XZ - ABS(Z)/4
  !        + 0.5*(2*XZ+ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
  !        + I*YZ/(1+X) * (XZ**2 + (YZ/(1+X))**2*R9ATN1(Y/(1+X)) ))
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R9ATN1, R9LN2R

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)

  REAL(SP) aipart, arg, cabsz, rpart, x, xz, y, y1x, yz
  COMPLEX(SP) Z
  !* FIRST EXECUTABLE STATEMENT  C9LN2R
  x = REAL(Z)
  y = AIMAG(Z)
  !
  cabsz = ABS(Z)
  IF ( cabsz>0.8125 ) THEN
    !
    C9LN2R = (LOG(1.0+Z)-Z*(1.0-0.5*Z))/Z**3
    RETURN
  END IF
  !
  C9LN2R = CMPLX(1.0/3.0,0.0)
  IF ( cabsz==0.0 ) RETURN
  !
  xz = x/cabsz
  yz = y/cabsz
  !
  arg = 2.0*xz + cabsz
  rpart = 0.5*arg**3*R9LN2R(cabsz*arg) - xz - 0.25*cabsz
  y1x = yz/(1.0+x)
  aipart = y1x*(xz**2+y1x**2*R9ATN1(cabsz*y1x))
  !
  C9LN2R = CMPLX(xz,-yz)**3*CMPLX(rpart,aipart)
  RETURN
END FUNCTION C9LN2R
