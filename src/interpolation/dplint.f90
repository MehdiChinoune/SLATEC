!** DPLINT
SUBROUTINE DPLINT(N,X,Y,C)
  IMPLICIT NONE
  !>
  !***
  !  Produce the polynomial which interpolates a set of discrete
  !            data points.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E1B
  !***
  ! **Type:**      DOUBLE PRECISION (POLINT-S, DPLINT-D)
  !***
  ! **Keywords:**  POLYNOMIAL INTERPOLATION
  !***
  ! **Author:**  Huddleston, R. E., (SNLL)
  !***
  ! **Description:**
  !
  !     Abstract
  !        Subroutine DPLINT is designed to produce the polynomial which
  !     interpolates the data  (X(I),Y(I)), I=1,...,N.  DPLINT sets up
  !     information in the array C which can be used by subroutine DPOLVL
  !     to evaluate the polynomial and its derivatives and by subroutine
  !     DPOLCF to produce the coefficients.
  !
  !     Formal Parameters
  !     *** All TYPE REAL variables are DOUBLE PRECISION ***
  !     N  - the number of data points  (N .GE. 1)
  !     X  - the array of abscissas (all of which must be distinct)
  !     Y  - the array of ordinates
  !     C  - an array of information used by subroutines
  !     *******  Dimensioning Information  *******
  !     Arrays X,Y, and C must be dimensioned at least N in the calling
  !     program.
  !
  !***
  ! **References:**  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER i, k, km1, N
  REAL(8) :: dif, C(*), X(*), Y(*)
  !* FIRST EXECUTABLE STATEMENT  DPLINT
  IF ( N<=0 ) THEN
    CALL XERMSG('SLATEC','DPLINT','N IS ZERO OR NEGATIVE.',2,1)
    RETURN
  ELSE
    C(1) = Y(1)
    IF ( N==1 ) RETURN
    DO k = 2, N
      C(k) = Y(k)
      km1 = k - 1
      DO i = 1, km1
        !     CHECK FOR DISTINCT X VALUES
        dif = X(i) - X(k)
        IF ( dif==0.0 ) GOTO 100
        C(k) = (C(i)-C(k))/dif
      ENDDO
    ENDDO
    RETURN
  ENDIF
  100  CALL XERMSG('SLATEC','DPLINT','THE ABSCISSAS ARE NOT DISTINCT.',2,1)
END SUBROUTINE DPLINT
