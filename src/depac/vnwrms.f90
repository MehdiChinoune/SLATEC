!** VNWRMS
REAL FUNCTION VNWRMS(N,V,W)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (VNWRMS-S, DVNRMS-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   VNWRMS computes a weighted root-mean-square vector norm for the
  !   integrator package DEBDF.
  !
  !***
  ! **See also:**  DEBDF
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  
  !
  !
  !LLL. OPTIMIZE
  !-----------------------------------------------------------------------
  ! THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM
  ! OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
  ! CONTAINED IN THE ARRAY W OF LENGTH N..
  !   VNWRMS = SQRT( (1/N) * SUM( V(I)/W(I) )**2 )
  !-----------------------------------------------------------------------
  INTEGER N, i
  REAL V, W, sum
  DIMENSION V(*), W(*)
  !* FIRST EXECUTABLE STATEMENT  VNWRMS
  sum = 0.0E0
  DO i = 1, N
    sum = sum + (V(i)/W(i))**2
  ENDDO
  VNWRMS = SQRT(sum/N)
  !----------------------- END OF FUNCTION VNWRMS ------------------------
END FUNCTION VNWRMS
