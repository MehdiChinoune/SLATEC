!*==DSVCO.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DSVCO
SUBROUTINE DSVCO(Rsav,Isav)
  IMPLICIT NONE
  !*--DSVCO5
  !***BEGIN PROLOGUE  DSVCO
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDEBDF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SVCO-S, DSVCO-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   DSVCO transfers data from a common block to arrays within the
  !   integrator package DDEBDF.
  !
  !***SEE ALSO  DDEBDF
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DDEBD1
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DSVCO
  !-----------------------------------------------------------------------
  ! THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
  ! DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF PACKAGE.
  !
  ! RSAV = DOUBLE PRECISION ARRAY OF LENGTH 218 OR MORE.
  ! ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
  !-----------------------------------------------------------------------
  INTEGER i , ILS , Isav , lenils , lenrls
  DOUBLE PRECISION RLS , Rsav
  DIMENSION Rsav(*) , Isav(*)
  SAVE lenrls , lenils
  COMMON /DDEBD1/ RLS(218) , ILS(33)
  DATA lenrls/218/ , lenils/33/
  !
  !***FIRST EXECUTABLE STATEMENT  DSVCO
  DO i = 1 , lenrls
    Rsav(i) = RLS(i)
  ENDDO
  DO i = 1 , lenils
    Isav(i) = ILS(i)
  ENDDO
  !     ----------------------- END OF SUBROUTINE DSVCO
  !     -----------------------
END SUBROUTINE DSVCO
