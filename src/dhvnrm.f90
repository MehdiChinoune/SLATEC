!*==DHVNRM.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DHVNRM
REAL(8) FUNCTION DHVNRM(V,Ncomp)
  IMPLICIT NONE
  !*--DHVNRM5
  !***BEGIN PROLOGUE  DHVNRM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     Compute the maximum norm of the vector V(*) of length NCOMP and
  !     return the result as DHVNRM
  !
  !***SEE ALSO  DDEABM, DDEBDF, DDERKF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891024  Changed references from DVNORM to DHVNRM.  (WRB)
  !   891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DHVNRM
  !
  INTEGER k , Ncomp
  REAL(8) :: V
  DIMENSION V(*)
  !***FIRST EXECUTABLE STATEMENT  DHVNRM
  DHVNRM = 0.0D0
  DO k = 1 , Ncomp
    DHVNRM = MAX(DHVNRM,ABS(V(k)))
  ENDDO
END FUNCTION DHVNRM
