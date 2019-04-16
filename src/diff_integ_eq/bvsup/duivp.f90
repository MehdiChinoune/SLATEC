!** DUIVP
SUBROUTINE DUIVP(X,Y,Yp)
  !>
  !***
  !  Dummy routine for DBVSUP quick check.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (UIVP-S, DUIVP-D)
  !***
  ! **Keywords:**  QUICK CHECK
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   This routine is never called;  it is here to prevent loaders from
  !   complaining about undefined externals while testing DBVSUP.
  !
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920401  Variables declaration and TYPE sections added.  (WRB)
  
  !     .. Scalar Arguments ..
  REAL(8) :: X
  !     .. Array Arguments ..
  REAL(8) :: Y(*), Yp(*)
  !* FIRST EXECUTABLE STATEMENT  DUIVP
  STOP
END SUBROUTINE DUIVP
