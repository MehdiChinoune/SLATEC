!** DBESKS
PURE SUBROUTINE DBESKS(Xnu,X,Nin,Bk)
  !> Compute a sequence of modified Bessel functions of the third kind of fractional order.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      DOUBLE PRECISION (BESKS-S, DBESKS-D)
  !***
  ! **Keywords:**  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION,
  !             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS, THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESKS computes a sequence of modified Bessel functions of the third
  ! kind of order XNU + I at X, where X > 0, XNU lies in (-1,1),
  ! and I = 0, 1, ..., NIN - 1, if NIN is positive and I = 0, 1, ... ,
  ! NIN + 1, if NIN is negative.  On return, the vector BK(.) contains
  ! the results at X for order starting at XNU.  XNU, X, and BK are
  ! double precision.  NIN is an integer.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DBSKES, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : D1MACH
  INTEGER, INTENT(IN) :: Nin
  REAL(DP), INTENT(IN) :: Xnu, X
  REAL(DP), INTENT(OUT) :: Bk(Nin)
  INTEGER :: i, n
  REAL(DP) :: expxi
  REAL(DP), PARAMETER :: xmax = -LOG(D1MACH(1))
  !* FIRST EXECUTABLE STATEMENT  DBESKS
  !
  IF( X>xmax ) ERROR STOP 'DBESKS : X SO BIG BESSEL K UNDERFLOWS'
  !
  CALL DBSKES(Xnu,X,Nin,Bk)
  !
  expxi = EXP(-X)
  n = ABS(Nin)
  DO i = 1, n
    Bk(i) = expxi*Bk(i)
  END DO
  !
END SUBROUTINE DBESKS