!*==DTEST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DTEST
SUBROUTINE DTEST(Len,Dcomp,Dtrue,Dsize,Dfac,Kprint)
  IMPLICIT NONE
  !*--DTEST5
  !*** Start of declarations inserted by SPAG
  INTEGER i, ICAse, INCx, INCy, Kprint, Len, MODe, N, NPRint
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DTEST
  !***PURPOSE  Compare arrays DCOMP and DTRUE.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (STEST-S, DTEST-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !   This subroutine compares arrays DCOMP and DTRUE of length LEN to
  !   see if the term by term differences, multiplied by DFAC, are
  !   negligible.  In the case of a significant difference, appropriate
  !   messages are written.
  !
  !***ROUTINES CALLED  D1MACH
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   741210  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900820  Modified IF test to use function DDIFF and made cosmetic
  !           changes to routine.  (WRB)
  !   901005  Removed usage of DDIFF in favour of D1MACH.  (RWC)
  !   910501  Added TYPE record.  (WRB)
  !   920211  Code restructured and information added to the DESCRIPTION
  !           section.  (WRB)
  !***END PROLOGUE  DTEST
  REAL(8) :: Dcomp(*), Dtrue(*), Dsize(*), Dfac, dd, releps, &
    D1MACH
  LOGICAL PASs
  COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
  SAVE releps
  DATA releps/0.0D0/
  !***FIRST EXECUTABLE STATEMENT  DTEST
  IF ( releps==0.0D0 ) releps = D1MACH(4)
  DO i = 1, Len
    dd = ABS(Dcomp(i)-Dtrue(i))
    IF ( Dfac*dd>ABS(Dsize(i))*releps ) THEN
      !
      !         Here DCOMP(I) is not close to DTRUE(I).
      !
      IF ( PASs ) THEN
        !
        !           Print FAIL message and header.
        !
        PASs = .FALSE.
        IF ( Kprint>=3 ) THEN
          WRITE (NPRint,99001)
          99001         FORMAT ('+',39X,'FAIL')
          WRITE (NPRint,99002)
          99002         FORMAT ('0CASE  N INCX INCY MODE  I',29X,'COMP(I)',29X,'TRUE(I)',&
            2X,'DIFFERENCE',5X,'SIZE(I)'/1X)
        ENDIF
      ENDIF
      IF ( Kprint>=3 ) WRITE (NPRint,99003) ICAse, N, INCx, INCy, MODe, &
        i, Dcomp(i), Dtrue(i), dd, Dsize(i)
      99003     FORMAT (1X,I4,I3,3I5,I3,2D36.18,2D12.4)
    ENDIF
  ENDDO
  RETURN
END SUBROUTINE DTEST
