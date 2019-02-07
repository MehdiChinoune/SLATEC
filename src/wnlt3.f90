!*==WNLT3.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK WNLT3
SUBROUTINE WNLT3(I,Imax,M,Mdw,Ipivot,H,W)
  IMPLICIT NONE
  !*--WNLT35
  !***BEGIN PROLOGUE  WNLT3
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to WNLIT
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (WNLT3-S, DWNLT3-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
  !
  !     Perform column interchange.
  !     Exchange elements of permuted index vector and perform column
  !     interchanges.
  !
  !***SEE ALSO  WNLIT
  !***ROUTINES CALLED  SSWAP
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLT and made a subroutine.  (RWC))
  !***END PROLOGUE  WNLT3
  INTEGER I , Imax , Ipivot(*) , M , Mdw
  REAL H(*) , W(Mdw,*)
  !
  EXTERNAL SSWAP
  !
  REAL t
  INTEGER itemp
  !
  !***FIRST EXECUTABLE STATEMENT  WNLT3
  IF ( Imax/=I ) THEN
    itemp = Ipivot(I)
    Ipivot(I) = Ipivot(Imax)
    Ipivot(Imax) = itemp
    !
    CALL SSWAP(M,W(1,Imax),1,W(1,I),1)
    !
    t = H(Imax)
    H(Imax) = H(I)
    H(I) = t
  ENDIF
END SUBROUTINE WNLT3
