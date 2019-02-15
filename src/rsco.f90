!DECK RSCO
SUBROUTINE RSCO(Rsav,Isav)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  RSCO
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DEBDF
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (RSCO-S, DRSCO-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   RSCO transfers data from arrays to a common block within the
  !   integrator package DEBDF.
  !
  !***SEE ALSO  DEBDF
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DEBDF1
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  RSCO
  !
  !
  !-----------------------------------------------------------------------
  ! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
  ! BLOCK DEBDF1 , WHICH IS USED INTERNALLY IN THE DEBDF
  ! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
  ! OF SUBROUTINE SVCO OR THE EQUIVALENT.
  !-----------------------------------------------------------------------
  INTEGER Isav, i, ILS, lenils, lenrls
  REAL Rsav, RLS
  DIMENSION Rsav(*), Isav(*)
  COMMON /DEBDF1/ RLS(218), ILS(33)
  SAVE lenrls, lenils
  DATA lenrls/218/, lenils/33/
  !
  !***FIRST EXECUTABLE STATEMENT  RSCO
  DO i = 1, lenrls
    RLS(i) = Rsav(i)
  ENDDO
  DO i = 1, lenils
    ILS(i) = Isav(i)
  ENDDO
  !----------------------- END OF SUBROUTINE RSCO -----------------------
END SUBROUTINE RSCO
