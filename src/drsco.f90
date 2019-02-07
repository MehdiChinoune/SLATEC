!*==DRSCO.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DRSCO
SUBROUTINE DRSCO(Rsav,Isav)
  IMPLICIT NONE
  !*--DRSCO5
  !***BEGIN PROLOGUE  DRSCO
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDEBDF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (RSCO-S, DRSCO-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   DRSCO transfers data from arrays to a common block within the
  !   integrator package DDEBDF.
  !
  !***SEE ALSO  DDEBDF
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DDEBD1
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DRSCO
  !-----------------------------------------------------------------------
  ! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
  ! BLOCK DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF
  ! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
  ! OF SUBROUTINE DSVCO OR THE EQUIVALENT.
  !-----------------------------------------------------------------------
  !
  INTEGER i , ILS , Isav , lenils , lenrls
  REAL(8) :: RLS , Rsav
  DIMENSION Rsav(*) , Isav(*)
  SAVE lenrls , lenils
  COMMON /DDEBD1/ RLS(218) , ILS(33)
  DATA lenrls/218/ , lenils/33/
  !
  !***FIRST EXECUTABLE STATEMENT  DRSCO
  DO i = 1 , lenrls
    RLS(i) = Rsav(i)
  ENDDO
  DO i = 1 , lenils
    ILS(i) = Isav(i)
  ENDDO
  !     ----------------------- END OF SUBROUTINE DRSCO
  !     -----------------------
END SUBROUTINE DRSCO
