!** UIVP
SUBROUTINE UIVP(X,Y,Yp)
  !>
  !***
  !  Dummy routine for BVSUP quick check.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (UIVP-S, DUIVP-D)
  !***
  ! **Keywords:**  QUICK CHECK
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   This routine is never called;  it is here to prevent loaders from
  !   complaining about undefined externals while testing BVSUP.
  !
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920401  Variables declaration and TYPE sections added.  (WRB)
  
  !     .. Scalar Arguments ..
  REAL X
  !     .. Array Arguments ..
  REAL Y(*), Yp(*)
  !* FIRST EXECUTABLE STATEMENT  UIVP
  STOP
END SUBROUTINE UIVP
