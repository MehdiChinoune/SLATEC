!** D9PAK
REAL(DP) ELEMENTAL FUNCTION D9PAK(Y,N)
  !> Pack a base 2 exponent into a floating point number.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  A6B
  !***
  ! **Type:**      DOUBLE PRECISION (R9PAK-S, D9PAK-D)
  !***
  ! **Keywords:**  FNLIB, PACK
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Pack a base 2 exponent into floating point number X.  This routine is
  ! almost the inverse of D9UPAK.  It is not exactly the inverse, because
  ! ABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N
  ! were known to be in range we could compute
  !               D9PAK = X *2.0d0**N
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9UPAK, I1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891009  Corrected error when XERROR called.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901009  Routine used I1MACH(7) where it should use I1MACH(10), Corrected (RWC)
  USE service, ONLY : D1MACH, I1MACH
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN) :: Y
  INTEGER :: nsum, ny
  REAL(DP), PARAMETER :: a1n2b = D1MACH(5)/LOG(2._DP)
  INTEGER, PARAMETER :: nmin = INT( a1n2b*I1MACH(15) ), nmax = INT( a1n2b*I1MACH(16) )
  !* FIRST EXECUTABLE STATEMENT  D9PAK
  !
  CALL D9UPAK(Y,D9PAK,ny)
  !
  nsum = N + ny
  IF( nsum<nmin ) THEN
    ! CALL XERMSG('D9PAK','PACKED NUMBER UNDERFLOWS',1,1)
    D9PAK = 0._DP
  ELSEIF( nsum>nmax ) THEN
    ERROR STOP 'D9PAK : PACKED NUMBER OVERFLOWS'
  ELSEIF( nsum==0 ) THEN
    RETURN
  ELSEIF( nsum>0 ) THEN
    DO
      D9PAK = 2._DP*D9PAK
      nsum = nsum - 1
      IF( nsum==0 ) EXIT
    END DO
  ELSE
    DO
      D9PAK = 0.5_DP*D9PAK
      nsum = nsum + 1
      IF( nsum==0 ) RETURN
    END DO
  END IF

  RETURN
END FUNCTION D9PAK