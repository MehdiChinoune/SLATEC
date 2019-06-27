!** R9PAK
REAL(SP) ELEMENTAL FUNCTION R9PAK(Y,N)
  !> Pack a base 2 exponent into a floating point number.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  A6B
  !***
  ! **Type:**      SINGLE PRECISION (R9PAK-S, D9PAK-D)
  !***
  ! **Keywords:**  FNLIB, PACK
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Pack a base 2 exponent into floating point number Y.  This
  ! routine is almost the inverse of R9UPAK.  It is not exactly
  ! the inverse, because ABS(X) need not be between 0.5 and
  ! 1.0.  If both R9PAK and 2.0**N were known to be in range, we
  ! could compute
  !       R9PAK = Y * 2.0**N .
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  I1MACH, R1MACH, R9UPAK, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901009  Routine used I1MACH(7) where it should use I1MACH(10), Corrected (RWC)
  USE service, ONLY : XERMSG, R1MACH, I1MACH
  INTEGER, INTENT(IN) :: N
  REAL(SP), INTENT(IN) :: Y
  INTEGER :: nsum, ny
  REAL(SP), PARAMETER :: a1n2b = R1MACH(5)/LOG(2._SP)
  INTEGER, PARAMETER :: nmin = INT( a1n2b*I1MACH(12) ), nmax = INT( a1n2b*I1MACH(13) )
  !* FIRST EXECUTABLE STATEMENT  R9PAK
  !
  CALL R9UPAK(Y,R9PAK,ny)
  !
  nsum = N + ny
  IF( nsum<nmin ) THEN
    ! CALL XERMSG('R9PAK : PACKED NUMBER UNDERFLOWS',1,1)
    R9PAK = 0._SP
  ELSEIF( nsum>nmax ) THEN
    ERROR STOP 'R9PAK : PACKED NUMBER OVERFLOWS'
  ELSEIF( nsum==0 ) THEN
    RETURN
  ELSEIF( nsum>0 ) THEN
    DO
      R9PAK = 2._SP*R9PAK
      nsum = nsum - 1
      IF( nsum==0 ) EXIT
    END DO
  ELSE
    DO
      R9PAK = 0.5_SP*R9PAK
      nsum = nsum + 1
      IF( nsum==0 ) RETURN
    END DO
  END IF

  RETURN
END FUNCTION R9PAK