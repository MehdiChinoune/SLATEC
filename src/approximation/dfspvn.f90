!** DFSPVN
SUBROUTINE DFSPVN(T,Jhigh,Indexx,X,Ileft,Vnikx)
  !> Subsidiary to DFC
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BSPLVN-S, DFSPVN-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  **** Double Precision version of BSPLVN ****
  !
  ! Calculates the value of all possibly nonzero B-splines at *X* of
  !  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*.
  !
  !***
  ! **See also:**  DFC
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  INTEGER, INTENT(IN) :: Ileft, Indexx, Jhigh
  REAL(DP), INTENT(IN) :: T(:), X
  REAL(DP), INTENT(INOUT) :: Vnikx(*)
  REAL(DP) :: vm, vmprev
  INTEGER :: imjp1, ipj, jp1, jp1ml, l
  INTEGER, SAVE :: j = 1
  REAL(DP), SAVE :: deltam(20) = 0._DP, deltap(20) = 0._DP
  !* FIRST EXECUTABLE STATEMENT  DFSPVN
  IF( Indexx/=2 ) THEN
    j = 1
    Vnikx(1) = 1._DP
    IF( j>=Jhigh ) RETURN
  END IF
  DO
    !
    ipj = Ileft + j
    deltap(j) = T(ipj) - X
    imjp1 = Ileft - j + 1
    deltam(j) = X - T(imjp1)
    vmprev = 0._DP
    jp1 = j + 1
    DO l = 1, j
      jp1ml = jp1 - l
      vm = Vnikx(l)/(deltap(l)+deltam(jp1ml))
      Vnikx(l) = vm*deltap(l) + vmprev
      vmprev = vm*deltam(jp1ml)
    END DO
    Vnikx(jp1) = vmprev
    j = jp1
    IF( j>=Jhigh ) EXIT
  END DO
  !
  RETURN
END SUBROUTINE DFSPVN