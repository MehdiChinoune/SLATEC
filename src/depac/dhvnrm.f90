!** DHVNRM
REAL(8) FUNCTION DHVNRM(V,Ncomp)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DDEABM, DDEBDF and DDERKF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !     Compute the maximum norm of the vector V(*) of length NCOMP and
  !     return the result as DHVNRM
  !
  !***
  ! **See also:**  DDEABM, DDEBDF, DDERKF
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891024  Changed references from DVNORM to DHVNRM.  (WRB)
  !   891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  
  !
  INTEGER k, Ncomp
  REAL(8) :: V(*)
  !* FIRST EXECUTABLE STATEMENT  DHVNRM
  DHVNRM = 0.0D0
  DO k = 1, Ncomp
    DHVNRM = MAX(DHVNRM,ABS(V(k)))
  ENDDO
END FUNCTION DHVNRM
