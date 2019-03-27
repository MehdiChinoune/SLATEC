!** RSCO
SUBROUTINE RSCO(Rsav,Isav)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (RSCO-S, DRSCO-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   RSCO transfers data from arrays to a common block within the
  !   integrator package DEBDF.
  !
  !***
  ! **See also:**  DEBDF
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DEBDF1

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  
  !
  !
  !-----------------------------------------------------------------------
  ! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
  ! BLOCK DEBDF1 , WHICH IS USED INTERNALLY IN THE DEBDF
  ! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
  ! OF SUBROUTINE SVCO OR THE EQUIVALENT.
  !-----------------------------------------------------------------------
  INTEGER Isav(*), i, ILS(33), lenils, lenrls
  REAL Rsav(*), RLS(218)
  COMMON /DEBDF1/ RLS, ILS
  SAVE lenrls, lenils
  DATA lenrls/218/, lenils/33/
  !
  !* FIRST EXECUTABLE STATEMENT  RSCO
  DO i = 1, lenrls
    RLS(i) = Rsav(i)
  ENDDO
  DO i = 1, lenils
    ILS(i) = Isav(i)
  ENDDO
  !----------------------- END OF SUBROUTINE RSCO -----------------------
END SUBROUTINE RSCO
