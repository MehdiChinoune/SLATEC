!*==STEST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK STEST
SUBROUTINE STEST(Len,Scomp,Strue,Ssize,Sfac,Kprint)
  IMPLICIT NONE
  !*--STEST5
  !*** Start of declarations inserted by SPAG
  INTEGER i, ICAse, INCx, INCy, Kprint, Len, MODe, N, NPRint
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  STEST
  !***PURPOSE  Compare arrays SCOMP and STRUE.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (STEST-S, DTEST-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !   This subroutine compares arrays SCOMP and STRUE of length LEN to
  !   see if the term by term differences, multiplied by SFAC, are
  !   negligible.  In the case of a significant difference, appropriate
  !   messages are written.
  !
  !***ROUTINES CALLED  R1MACH
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   741210  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900820  Modified IF test to use function DIFF and made cosmetic
  !           changes to routine.  (WRB)
  !   901005  Removed usage of DIFF in favour of R1MACH.  (RWC)
  !   910501  Added TYPE record.  (WRB)
  !   920211  Code restructured and information added to the DESCRIPTION
  !           section.  (WRB)
  !***END PROLOGUE  STEST
  REAL Scomp(*), Strue(*), Ssize(*), Sfac, sd, releps, R1MACH
  LOGICAL PASs
  COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
  SAVE releps
  DATA releps/0.0E0/
  !***FIRST EXECUTABLE STATEMENT  STEST
  IF ( releps==0.0E0 ) releps = R1MACH(4)
  DO i = 1, Len
    sd = ABS(Scomp(i)-Strue(i))
    IF ( Sfac*sd>ABS(Ssize(i))*releps ) THEN
      !
      !         Here SCOMP(I) is not close to STRUE(I).
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
        i, Scomp(i), Strue(i), sd, Ssize(i)
      99003     FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
    ENDIF
  ENDDO
  RETURN
END SUBROUTINE STEST
