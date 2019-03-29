!** DRSCO
SUBROUTINE DRSCO(Rsav,Isav)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (RSCO-S, DRSCO-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   DRSCO transfers data from arrays to a common block within the
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
  !   910722  Updated AUTHOR section.  (ALS)

  !-----------------------------------------------------------------------
  ! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
  ! BLOCK DDEBD1 , WHICH IS USED INTERNALLY IN THE DDEBDF
  ! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
  ! OF SUBROUTINE DSVCO OR THE EQUIVALENT.
  !-----------------------------------------------------------------------
  !
  INTEGER i, ILS, Isav(*)
  REAL(8) :: RLS, Rsav(*)
  COMMON /DDEBD1/ RLS(218), ILS(33)
  INTEGER, PARAMETER :: lenrls = 218, lenils = 33
  !
  !* FIRST EXECUTABLE STATEMENT  DRSCO
  DO i = 1, lenrls
    RLS(i) = Rsav(i)
  ENDDO
  DO i = 1, lenils
    ILS(i) = Isav(i)
  ENDDO
  !     ----------------------- END OF SUBROUTINE DRSCO
  !     -----------------------
END SUBROUTINE DRSCO
