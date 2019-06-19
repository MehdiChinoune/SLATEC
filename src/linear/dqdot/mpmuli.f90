!** MPMULI
SUBROUTINE MPMULI(X,Iy,Z)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPMULI-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! Multiplies 'mp' X by single-precision integer IY giving 'mp' Z.
  ! This is faster than using MPMUL.  Result is ROUNDED.
  ! Multiplication by 1 may be used to normalize a number
  ! even if the last digit is B.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI
  !***
  ! **Routines called:**  MPMUL2

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Iy
  INTEGER :: X(*), Z(*)
  !* FIRST EXECUTABLE STATEMENT  MPMULI
  CALL MPMUL2(X,Iy,Z,0)
END SUBROUTINE MPMULI
