!DECK PCHST
REAL FUNCTION PCHST(Arg1,Arg2)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PCHST
  !***SUBSIDIARY
  !***PURPOSE  PCHIP Sign-Testing Routine
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      SINGLE PRECISION (PCHST-S, DPCHST-D)
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
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
  !***SEE ALSO  PCHCE, PCHCI, PCHCS, PCHIM
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   811103  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870813  Minor cosmetic changes.
  !   890411  Added SAVE statements (Vers. 3.2).
  !   890411  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   930503  Improved purpose.  (FNF)
  !***END PROLOGUE  PCHST
  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  REAL Arg1, Arg2
  !
  !  DECLARE LOCAL VARIABLES.
  !
  REAL one, zero
  SAVE zero, one
  DATA zero/0./, one/1./
  !
  !  PERFORM THE TEST.
  !
  !***FIRST EXECUTABLE STATEMENT  PCHST
  PCHST = SIGN(one,Arg1)*SIGN(one,Arg2)
  IF ( (Arg1==zero).OR.(Arg2==zero) ) PCHST = zero
  !
  !------------- LAST LINE OF PCHST FOLLOWS ------------------------------
END FUNCTION PCHST
