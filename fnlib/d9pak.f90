!*==D9PAK.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9PAK
REAL(8) FUNCTION D9PAK(Y,N)
  IMPLICIT NONE
  !*--D9PAK5
  !*** Start of declarations inserted by SPAG
  INTEGER I1MACH, N, nmax, nmin, nsum, ny
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9PAK
  !***PURPOSE  Pack a base 2 exponent into a floating point number.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  A6B
  !***TYPE      DOUBLE PRECISION (R9PAK-S, D9PAK-D)
  !***KEYWORDS  FNLIB, PACK
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Pack a base 2 exponent into floating point number X.  This routine is
  ! almost the inverse of D9UPAK.  It is not exactly the inverse, because
  ! ABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N
  ! were known to be in range we could compute
  !               D9PAK = X *2.0d0**N
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9UPAK, I1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891009  Corrected error when XERROR called.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901009  Routine used I1MACH(7) where it should use I1MACH(10),
  !           Corrected (RWC)
  !***END PROLOGUE  D9PAK
  REAL(8) :: Y, a1n2b, a1n210, D1MACH
  LOGICAL first
  SAVE nmin, nmax, a1n210, first
  DATA a1n210/3.321928094887362347870319429489D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9PAK
  IF ( first ) THEN
    a1n2b = 1.0D0
    IF ( I1MACH(10)/=2 ) a1n2b = D1MACH(5)*a1n210
    nmin = a1n2b*I1MACH(15)
    nmax = a1n2b*I1MACH(16)
  ENDIF
  first = .FALSE.
  !
  CALL D9UPAK(Y,D9PAK,ny)
  !
  nsum = N + ny
  IF ( nsum<nmin ) THEN
    !
    CALL XERMSG('SLATEC','D9PAK','PACKED NUMBER UNDERFLOWS',1,1)
    D9PAK = 0.0D0
    GOTO 99999
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
      ENDDO
    ELSE
      DO
        !
        D9PAK = 0.5D0*D9PAK
        nsum = nsum + 1
        IF ( nsum==0 ) RETURN
      ENDDO
    ENDIF
  ENDIF
  RETURN
  !
  99999 CONTINUE
  END FUNCTION D9PAK
