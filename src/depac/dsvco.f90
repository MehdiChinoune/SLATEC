!** DSVCO
SUBROUTINE DSVCO(Rsav,Isav)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SVCO-S, DSVCO-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   DSVCO transfers data from a common block to arrays within the
  !   integrator package DDEBDF.
  !
  !***
  ! **See also:**  DDEBDF
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DDEBD1

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !-----------------------------------------------------------------------
  ! THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
  ! DDEBD1 , WHICH IS USED INTERNALLY IN THE DDEBDF PACKAGE.
  !
  ! RSAV = DOUBLE PRECISION ARRAY OF LENGTH 218 OR MORE.
  ! ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
  !-----------------------------------------------------------------------
  INTEGER i, ILS, Isav(*)
  REAL(8) :: RLS, Rsav(*)
  COMMON /DDEBD1/ RLS(218), ILS(33)
  INTEGER, PARAMETER :: lenrls = 218, lenils = 33
  !
  !* FIRST EXECUTABLE STATEMENT  DSVCO
  DO i = 1, lenrls
    Rsav(i) = RLS(i)
  ENDDO
  DO i = 1, lenils
    Isav(i) = ILS(i)
  ENDDO
  !     ----------------------- END OF SUBROUTINE DSVCO
  !     -----------------------
END SUBROUTINE DSVCO
