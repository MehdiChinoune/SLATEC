!DECK SCOPYM
SUBROUTINE SCOPYM(N,Sx,Incx,Sy,Incy)
  IMPLICIT NONE
  INTEGER i, Incx, Incy, ix, iy, m, mp1, N, ns
  !***BEGIN PROLOGUE  SCOPYM
  !***PURPOSE  Copy the negative of a vector to a vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A5
  !***TYPE      SINGLE PRECISION (SCOPYM-S, DCOPYM-D)
  !***KEYWORDS  BLAS, COPY, VECTOR
  !***AUTHOR  Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  !       Description of Parameters
  !           The * Flags Output Variables
  !
  !       N   Number of elements in vector(s)
  !      SX   Real vector with N elements
  !    INCX   Storage spacing between elements of SX
  !      SY*  Real negative copy of SX
  !    INCY   Storage spacing between elements of SY
  !
  !      ***  Note that SY = -SX  ***
  !
  !     Copy negative of real SX to real SY.  For I=0 to N-1,
  !     copy  -SX(LX+I*INCX) to SY(LY+I*INCY), where LX=1 if
  !     INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is defined
  !     in a similar way using INCY.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !***END PROLOGUE  SCOPYM
  REAL Sx(*), Sy(*)
  !***FIRST EXECUTABLE STATEMENT  SCOPYM
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy ) THEN
    IF ( Incx<1 ) THEN
    ELSEIF ( Incx==1 ) THEN
      !
      !     Code for both increments equal to 1.
      !
      !     Clean-up loop so remaining vector length is a multiple of 7.
      !
      m = MOD(N,7)
      IF ( m/=0 ) THEN
        DO i = 1, m
          Sy(i) = -Sx(i)
        ENDDO
        IF ( N<7 ) RETURN
      ENDIF
      GOTO 100
    ELSE
      !
      !     Code for equal, positive, non-unit increments.
      !
      ns = N*Incx
      DO i = 1, ns, Incx
        Sy(i) = -Sx(i)
      ENDDO
      RETURN
    ENDIF
  ENDIF
  !
  !     Code for unequal or nonpositive increments.
  !
  ix = 1
  iy = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  IF ( Incy<0 ) iy = (-N+1)*Incy + 1
  DO i = 1, N
    Sy(iy) = -Sx(ix)
    ix = ix + Incx
    iy = iy + Incy
  ENDDO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 7
    Sy(i) = -Sx(i)
    Sy(i+1) = -Sx(i+1)
    Sy(i+2) = -Sx(i+2)
    Sy(i+3) = -Sx(i+3)
    Sy(i+4) = -Sx(i+4)
    Sy(i+5) = -Sx(i+5)
    Sy(i+6) = -Sx(i+6)
  ENDDO
  RETURN
END SUBROUTINE SCOPYM
