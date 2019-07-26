!** PCHST
REAL(SP) ELEMENTAL FUNCTION PCHST(Arg1,Arg2)
  !> PCHIP Sign-Testing Routine
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      SINGLE PRECISION (PCHST-S, DPCHST-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !         PCHST:  PCHIP Sign-Testing Routine.
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
  ! **See also:**  PCHCE, PCHCI, PCHCS, PCHIM
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811103  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890411  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   930503  Improved purpose.  (FNF)

  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  REAL(SP), INTENT(IN) :: Arg1, Arg2
  !
  !  PERFORM THE TEST.
  !
  !* FIRST EXECUTABLE STATEMENT  PCHST
  PCHST = SIGN(1._SP,Arg1)*SIGN(1._SP,Arg2)
  IF( (Arg1==0._SP) .OR. (Arg2==0._SP) ) PCHST = 0._SP
  !
  !------------- LAST LINE OF PCHST FOLLOWS ------------------------------
END FUNCTION PCHST