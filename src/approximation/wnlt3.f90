!** WNLT3
SUBROUTINE WNLT3(I,Imax,M,Mdw,Ipivot,H,W)
  !>
  !  Subsidiary to WNLIT
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (WNLT3-S, DWNLT3-D)
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
  ! **See also:**  WNLIT
  !***
  ! **Routines called:**  SSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLT and made a subroutine.  (RWC))
  USE linear, ONLY : SSWAP
  INTEGER :: I, Imax, M, Mdw, Ipivot(:)
  REAL(SP) :: H(:), W(Mdw,M)
  !
  INTEGER :: itemp
  REAL(SP) :: t
  !
  !* FIRST EXECUTABLE STATEMENT  WNLT3
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
  END IF
END SUBROUTINE WNLT3
