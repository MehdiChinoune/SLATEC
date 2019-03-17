!DECK DBESKS
SUBROUTINE DBESKS(Xnu,X,Nin,Bk)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DBESKS
  !***PURPOSE  Compute a sequence of modified Bessel functions of the
  !            third kind of fractional order.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B3
  !***TYPE      DOUBLE PRECISION (BESKS-S, DBESKS-D)
  !***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
  !             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBESKS computes a sequence of modified Bessel functions of the third
  ! kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1),
  ! and I = 0, 1, ..., NIN - 1, if NIN is positive and I = 0, 1, ... ,
  ! NIN + 1, if NIN is negative.  On return, the vector BK(.) contains
  ! the results at X for order starting at XNU.  XNU, X, and BK are
  ! double precision.  NIN is an integer.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DBSKES, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBESKS
  INTEGER i, n, Nin
  REAL(8) :: Xnu, X, Bk(*), expxi, xmax, D1MACH
  SAVE xmax
  DATA xmax/0.D0/
  !***FIRST EXECUTABLE STATEMENT  DBESKS
  IF ( xmax==0.D0 ) xmax = -LOG(D1MACH(1))
  !
  IF ( X>xmax ) CALL XERMSG('SLATEC','DBESKS','X SO BIG BESSEL K UNDERFLOWS'&
    ,1,2)
  !
  CALL DBSKES(Xnu,X,Nin,Bk)
  !
  expxi = EXP(-X)
  n = ABS(Nin)
  DO i = 1, n
    Bk(i) = expxi*Bk(i)
  ENDDO
  !
END SUBROUTINE DBESKS
