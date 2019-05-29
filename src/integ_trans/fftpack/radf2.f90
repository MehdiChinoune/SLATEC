!** RADF2
SUBROUTINE RADF2(Ido,L1,Cc,Ch,Wa1)
  !>
  !  Calculate the fast Fourier transform of subvectors of
  !            length two.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADF2-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           changing dummy array size declarations (1) to (*).
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Ido, L1
  REAL :: Cc(Ido,L1,2), Ch(Ido,2,L1), Wa1(Ido)
  INTEGER :: i, ic, idp2, k
  REAL :: ti2, tr2
  !* FIRST EXECUTABLE STATEMENT  RADF2
  DO k = 1, L1
    Ch(1,1,k) = Cc(1,k,1) + Cc(1,k,2)
    Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,2)
  END DO
  IF ( Ido<2 ) RETURN
  IF ( Ido/=2 ) THEN
    idp2 = Ido + 2
    IF ( (Ido-1)/2<L1 ) THEN
      DO i = 3, Ido, 2
        ic = idp2 - i
        DO k = 1, L1
          tr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
          ti2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
          Ch(i,1,k) = Cc(i,k,1) + ti2
          Ch(ic,2,k) = ti2 - Cc(i,k,1)
          Ch(i-1,1,k) = Cc(i-1,k,1) + tr2
          Ch(ic-1,2,k) = Cc(i-1,k,1) - tr2
        END DO
      END DO
    ELSE
      DO k = 1, L1
        DO i = 3, Ido, 2
          ic = idp2 - i
          tr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
          ti2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
          Ch(i,1,k) = Cc(i,k,1) + ti2
          Ch(ic,2,k) = ti2 - Cc(i,k,1)
          Ch(i-1,1,k) = Cc(i-1,k,1) + tr2
          Ch(ic-1,2,k) = Cc(i-1,k,1) - tr2
        END DO
      END DO
    END IF
    IF ( MOD(Ido,2)==1 ) RETURN
  END IF
  DO k = 1, L1
    Ch(1,2,k) = -Cc(Ido,k,2)
    Ch(Ido,1,k) = Cc(Ido,k,1)
  END DO
  RETURN
END SUBROUTINE RADF2
