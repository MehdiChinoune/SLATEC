!** DWNLT3
SUBROUTINE DWNLT3(I,Imax,M,Mdw,Ipivot,H,W)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to WNLIT
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !     Perform column interchange.
  !     Exchange elements of permuted index vector and perform column
  !     interchanges.
  !
  !***
  ! **See also:**  DWNLIT
  !***
  ! **Routines called:**  DSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
  !   900604  DP version created from SP version.  (RWC)

  INTEGER I, Imax, Ipivot(*), M, Mdw
  REAL(8) :: H(*), W(Mdw,*)
  !
  EXTERNAL :: DSWAP
  !
  REAL(8) :: t
  INTEGER itemp
  !
  !* FIRST EXECUTABLE STATEMENT  DWNLT3
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
  END IF
END SUBROUTINE DWNLT3
