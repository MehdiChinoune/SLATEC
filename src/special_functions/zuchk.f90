!** ZUCHK
SUBROUTINE ZUCHK(Yr,Yi,Nz,Ascle,Tol)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and
  !            ZKSCL
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
  !      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
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
  
  !
  !     COMPLEX Y
  REAL(8) :: Ascle, ss, st, Tol, wr, wi, Yr, Yi
  INTEGER Nz
  !* FIRST EXECUTABLE STATEMENT  ZUCHK
  Nz = 0
  wr = ABS(Yr)
  wi = ABS(Yi)
  st = MIN(wr,wi)
  IF ( st>Ascle ) RETURN
  ss = MAX(wr,wi)
  st = st/Tol
  IF ( ss<st ) Nz = 1
END SUBROUTINE ZUCHK
