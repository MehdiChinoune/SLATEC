!DECK R9PAK
FUNCTION R9PAK(Y,N)
  IMPLICIT NONE
  REAL a1n210, a1n2b, R1MACH, R9PAK, Y
  INTEGER I1MACH, N, nmax, nmin, nsum, ny
  !***BEGIN PROLOGUE  R9PAK
  !***PURPOSE  Pack a base 2 exponent into a floating point number.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  A6B
  !***TYPE      SINGLE PRECISION (R9PAK-S, D9PAK-D)
  !***KEYWORDS  FNLIB, PACK
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Pack a base 2 exponent into floating point number Y.  This
  ! routine is almost the inverse of R9UPAK.  It is not exactly
  ! the inverse, because ABS(X) need not be between 0.5 and
  ! 1.0.  If both R9PAK and 2.0**N were known to be in range, we
  ! could compute
  !       R9PAK = Y * 2.0**N .
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  I1MACH, R1MACH, R9UPAK, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901009  Routine used I1MACH(7) where it should use I1MACH(10),
  !           Corrected (RWC)
  !***END PROLOGUE  R9PAK
  LOGICAL first
  SAVE nmin, nmax, a1n210, first
  DATA a1n210/3.321928094887362E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  R9PAK
  IF ( first ) THEN
    a1n2b = 1.0
    IF ( I1MACH(10)/=2 ) a1n2b = R1MACH(5)*a1n210
    nmin = a1n2b*I1MACH(12)
    nmax = a1n2b*I1MACH(13)
  ENDIF
  first = .FALSE.
  !
  CALL R9UPAK(Y,R9PAK,ny)
  !
  nsum = N + ny
  IF ( nsum<nmin ) THEN
    !
    CALL XERMSG('SLATEC','R9PAK','PACKED NUMBER UNDERFLOWS',1,1)
    R9PAK = 0.0
    GOTO 99999
  ELSE
    IF ( nsum>nmax ) CALL XERMSG('SLATEC','R9PAK','PACKED NUMBER OVERFLOWS',&
      2,2)
    !
    IF ( nsum==0 ) RETURN
    IF ( nsum>0 ) THEN
      DO
        !
        R9PAK = 2.0*R9PAK
        nsum = nsum - 1
        IF ( nsum==0 ) EXIT
      ENDDO
    ELSE
      DO
        !
        R9PAK = 0.5*R9PAK
        nsum = nsum + 1
        IF ( nsum==0 ) RETURN
      ENDDO
    ENDIF
  ENDIF
  RETURN
  !
  99999 CONTINUE
  END FUNCTION R9PAK
