!*==DWNLT3.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DWNLT3
SUBROUTINE DWNLT3(I,Imax,M,Mdw,Ipivot,H,W)
  IMPLICIT NONE
  !*--DWNLT35
  !***BEGIN PROLOGUE  DWNLT3
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to WNLIT
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
  !
  !     Perform column interchange.
  !     Exchange elements of permuted index vector and perform column
  !     interchanges.
  !
  !***SEE ALSO  DWNLIT
  !***ROUTINES CALLED  DSWAP
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
  !   900604  DP version created from SP version.  (RWC)
  !***END PROLOGUE  DWNLT3
  INTEGER I , Imax , Ipivot(*) , M , Mdw
  REAL(8) :: H(*) , W(Mdw,*)
  !
  EXTERNAL DSWAP
  !
  REAL(8) :: t
  INTEGER itemp
  !
  !***FIRST EXECUTABLE STATEMENT  DWNLT3
  IF ( Imax/=I ) THEN
    itemp = Ipivot(I)
    Ipivot(I) = Ipivot(Imax)
    Ipivot(Imax) = itemp
    !
    CALL DSWAP(M,W(1,Imax),1,W(1,I),1)
    !
    t = H(Imax)
    H(Imax) = H(I)
    H(I) = t
  ENDIF
END SUBROUTINE DWNLT3
