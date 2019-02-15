!DECK XGETF
SUBROUTINE XGETF(Kontrl)
  IMPLICIT NONE
  INTEGER J4SAVE, Kontrl
  !***BEGIN PROLOGUE  XGETF
  !***PURPOSE  Return the current value of the error control flag.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XGETF-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !   Abstract
  !        XGETF returns the current value of the error control flag
  !        in KONTRL.  See subroutine XSETF for flag value meanings.
  !        (KONTRL is an output parameter only.)
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  J4SAVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XGETF
  !***FIRST EXECUTABLE STATEMENT  XGETF
  Kontrl = J4SAVE(2,0,.FALSE.)
END SUBROUTINE XGETF
