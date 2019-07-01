!** DPCHST
REAL(DP) ELEMENTAL FUNCTION DPCHST(Arg1,Arg2)
  !> DPCHIP Sign-Testing Routine
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      DOUBLE PRECISION (PCHST-S, DPCHST-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !         DPCHST:  DPCHIP Sign-Testing Routine.
  !
  !
  !     Returns:
  !        -1. if ARG1 and ARG2 are of opposite sign.
  !         0. if either argument is zero.
  !        +1. if ARG1 and ARG2 are of the same sign.
  !
  !     The object is to do this without multiplying ARG1*ARG2, to avoid
  !     possible over/underflow problems.
  !
  !  Fortran intrinsics used:  SIGN.
  !
  !***
  ! **See also:**  DPCHCE, DPCHCI, DPCHCS, DPCHIM
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811103  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   930503  Improved purpose.  (FNF)

  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  REAL(DP), INTENT(IN) :: Arg1, Arg2
  !
  !  DECLARE LOCAL VARIABLES.
  !
  REAL(DP), PARAMETER :: zero = 0._DP, one = 1._DP
  !
  !  PERFORM THE TEST.
  !
  !* FIRST EXECUTABLE STATEMENT  DPCHST
  DPCHST = SIGN(one,Arg1)*SIGN(one,Arg2)
  IF( (Arg1==zero) .OR. (Arg2==zero) ) DPCHST = zero
  !
  !------------- LAST LINE OF DPCHST FOLLOWS -----------------------------
END FUNCTION DPCHST