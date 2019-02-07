!*==ITEST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK ITEST
SUBROUTINE ITEST(Len,Icomp,Itrue,Kprint)
  IMPLICIT NONE
  !*--ITEST5
  !*** Start of declarations inserted by SPAG
  INTEGER i , ICAse , id , INCx , INCy , Kprint , Len , MODe , N , NPRint
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  ITEST
  !***PURPOSE  Compare arrays ICOMP and ITRUE.
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (ITEST-I)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !   This subroutine compares the arrays ICOMP and ITRUE of length LEN
  !   for equality.  In the case of an unequal compare, appropriate
  !   messages are written.
  !
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   741210  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920211  Code restructured and information added to the DESCRIPTION
  !           section.  (WRB)
  !***END PROLOGUE  ITEST
  INTEGER Icomp(*) , Itrue(*)
  LOGICAL PASs
  COMMON /COMBLA/ NPRint , ICAse , N , INCx , INCy , MODe , PASs
  !***FIRST EXECUTABLE STATEMENT  ITEST
  DO i = 1 , Len
    IF ( Icomp(i)/=Itrue(i) ) THEN
      !
      !         Here ICOMP(I) is not equal to ITRUE(I).
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
            2X,'DIFFERENCE'/1X)
        ENDIF
      ENDIF
      IF ( Kprint>=3 ) THEN
        id = Icomp(i) - Itrue(i)
        WRITE (NPRint,99003) ICAse , N , INCx , INCy , MODe , i , Icomp(i) , &
          Itrue(i) , id
        99003       FORMAT (1X,I4,I3,3I5,I3,2I36,I12)
      ENDIF
    ENDIF
  ENDDO
  RETURN
END SUBROUTINE ITEST
