!** SVCO
SUBROUTINE SVCO(Rsav,Isav)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SVCO-S, DSVCO-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   SVCO transfers data from a common block to arrays within the
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
  
  !
  !
  !-----------------------------------------------------------------------
  ! THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
  ! DEBDF1 , WHICH IS USED INTERNALLY IN THE DEBDF PACKAGE.
  !
  ! RSAV = REAL ARRAY OF LENGTH 218 OR MORE.
  ! ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
  !-----------------------------------------------------------------------
  INTEGER Isav(*), i, ILS(33), lenils, lenrls
  REAL Rsav(*), RLS(218)
  COMMON /DEBDF1/ RLS, ILS
  SAVE lenrls, lenils
  DATA lenrls/218/, lenils/33/
  !
  !* FIRST EXECUTABLE STATEMENT  SVCO
  DO i = 1, lenrls
    Rsav(i) = RLS(i)
  ENDDO
  DO i = 1, lenils
    Isav(i) = ILS(i)
  ENDDO
  !----------------------- END OF SUBROUTINE SVCO -----------------------
END SUBROUTINE SVCO
