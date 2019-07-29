!** ZUCHK
ELEMENTAL SUBROUTINE ZUCHK(Y,Nz,Ascle,Tol)
  !> Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and ZKSCL
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUCHK-A, ZUCHK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
  !      EXP(-ALIM)=ASCLE=1.0E+3*tiny_dp/TOL. THE TEST IS MADE TO SEE
  !      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
  !      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
  !      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
  !      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
  !      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
  !
  !***
  ! **See also:**  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Ascle, Tol
  COMPLEX(DP), INTENT(IN) :: Y
  !
  REAL(DP) :: ss, st, yr, yi
  !* FIRST EXECUTABLE STATEMENT  ZUCHK
  Nz = 0
  yr = REAL(Y,DP)
  yi = AIMAG(Y)
  yr = ABS(yr)
  yi = ABS(yi)
  st = MIN(yr,yi)
  IF( st>Ascle ) RETURN
  ss = MAX(yr,yi)
  st = st/Tol
  IF( ss<st ) Nz = 1
  !
END SUBROUTINE ZUCHK