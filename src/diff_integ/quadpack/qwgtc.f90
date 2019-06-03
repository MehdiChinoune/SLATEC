!** QWGTC
REAL(SP) FUNCTION QWGTC(X,C,P2,P3,P4,Kp)
  !>
  !  This function subprogram is used together with the
  !            routine QAWC and defines the WEIGHT function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (QWGTC-S, DQWGTC-D)
  !***
  ! **Keywords:**  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **See also:**  QK15W
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   830518  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  REAL(SP) C, P2, P3, P4, X
  INTEGER Kp
  !* FIRST EXECUTABLE STATEMENT  QWGTC
  QWGTC = 0.1E+01/(X-C)
END FUNCTION QWGTC
