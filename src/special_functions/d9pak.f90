!** D9PAK
REAL(8) FUNCTION D9PAK(Y,N)
  !>
  !  Pack a base 2 exponent into a floating point number.
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
  !   901009  Routine used I1MACH(7) where it should use I1MACH(10),
  !           Corrected (RWC)
  USE service, ONLY : XERMSG, D1MACH, I1MACH
  INTEGER N, nsum, ny
  REAL(8) :: Y, a1n2b
  INTEGER, SAVE :: nmin, nmax
  REAL(8), PARAMETER :: a1n210 = 3.321928094887362347870319429489D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9PAK
  IF ( first ) THEN
    a1n2b = 1.0D0
    IF ( I1MACH(10)/=2 ) a1n2b = D1MACH(5)*a1n210
    nmin = INT( a1n2b*I1MACH(15) )
    nmax = INT( a1n2b*I1MACH(16) )
    first = .FALSE.
  END IF
  !
  CALL D9UPAK(Y,D9PAK,ny)
  !
  nsum = N + ny
  IF ( nsum<nmin ) THEN
    !
    CALL XERMSG('SLATEC','D9PAK','PACKED NUMBER UNDERFLOWS',1,1)
    D9PAK = 0.0D0
    RETURN
  ELSE
    IF ( nsum>nmax ) CALL XERMSG('SLATEC','D9PAK','PACKED NUMBER OVERFLOWS',&
      1,2)
    !
    IF ( nsum==0 ) RETURN
    IF ( nsum>0 ) THEN
      DO
        !
        D9PAK = 2.0D0*D9PAK
        nsum = nsum - 1
        IF ( nsum==0 ) EXIT
      END DO
    ELSE
      DO
        !
        D9PAK = 0.5D0*D9PAK
        nsum = nsum + 1
        IF ( nsum==0 ) RETURN
      END DO
    END IF
  END IF
  RETURN
END FUNCTION D9PAK
