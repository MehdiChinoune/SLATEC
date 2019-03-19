!** DVNRMS
REAL(8) FUNCTION DVNRMS(N,V,W)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   DVNRMS computes a weighted root-mean-square vector norm for the
  !   integrator package DDEBDF.
  !
  !***
  ! **See also:**  DDEBDF
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  
  INTEGER i, N
  REAL(8) :: sum, V, W
  DIMENSION V(*), W(*)
  !* FIRST EXECUTABLE STATEMENT  DVNRMS
  sum = 0.0D0
  DO i = 1, N
    sum = sum + (V(i)/W(i))**2
  ENDDO
  DVNRMS = SQRT(sum/N)
  !     ----------------------- END OF FUNCTION DVNRMS
  !     ------------------------
END FUNCTION DVNRMS
