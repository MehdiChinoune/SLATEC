!** RADF4
SUBROUTINE RADF4(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3)
  !> Calculate the fast Fourier transform of subvectors of
  !            length four.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADF4-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*).
  !           (b) changing definition of variable HSQT2 by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Ido, L1
  REAL(SP) :: Cc(Ido,L1,4), Ch(Ido,4,L1), Wa1(Ido), Wa2(Ido), Wa3(ido)
  INTEGER :: i, ic, idp2, k
  REAL(SP) :: ci2, ci3, ci4, cr2, cr3, cr4, hsqt2, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
  !* FIRST EXECUTABLE STATEMENT  RADF4
  hsqt2 = 0.5_SP*SQRT(2._SP)
  DO k = 1, L1
    tr1 = Cc(1,k,2) + Cc(1,k,4)
    tr2 = Cc(1,k,1) + Cc(1,k,3)
    Ch(1,1,k) = tr1 + tr2
    Ch(Ido,4,k) = tr2 - tr1
    Ch(Ido,2,k) = Cc(1,k,1) - Cc(1,k,3)
    Ch(1,3,k) = Cc(1,k,4) - Cc(1,k,2)
  END DO
  IF( Ido<2 ) RETURN
  IF( Ido/=2 ) THEN
    idp2 = Ido + 2
    IF( (Ido-1)/2<L1 ) THEN
      DO i = 3, Ido, 2
        ic = idp2 - i
        DO k = 1, L1
          cr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
          ci2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
          cr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
          ci3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
          cr4 = Wa3(i-2)*Cc(i-1,k,4) + Wa3(i-1)*Cc(i,k,4)
          ci4 = Wa3(i-2)*Cc(i,k,4) - Wa3(i-1)*Cc(i-1,k,4)
          tr1 = cr2 + cr4
          tr4 = cr4 - cr2
          ti1 = ci2 + ci4
          ti4 = ci2 - ci4
          ti2 = Cc(i,k,1) + ci3
          ti3 = Cc(i,k,1) - ci3
          tr2 = Cc(i-1,k,1) + cr3
          tr3 = Cc(i-1,k,1) - cr3
          Ch(i-1,1,k) = tr1 + tr2
          Ch(ic-1,4,k) = tr2 - tr1
          Ch(i,1,k) = ti1 + ti2
          Ch(ic,4,k) = ti1 - ti2
          Ch(i-1,3,k) = ti4 + tr3
          Ch(ic-1,2,k) = tr3 - ti4
          Ch(i,3,k) = tr4 + ti3
          Ch(ic,2,k) = tr4 - ti3
        END DO
      END DO
    ELSE
      DO k = 1, L1
        DO i = 3, Ido, 2
          ic = idp2 - i
          cr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
          ci2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
          cr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
          ci3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
          cr4 = Wa3(i-2)*Cc(i-1,k,4) + Wa3(i-1)*Cc(i,k,4)
          ci4 = Wa3(i-2)*Cc(i,k,4) - Wa3(i-1)*Cc(i-1,k,4)
          tr1 = cr2 + cr4
          tr4 = cr4 - cr2
          ti1 = ci2 + ci4
          ti4 = ci2 - ci4
          ti2 = Cc(i,k,1) + ci3
          ti3 = Cc(i,k,1) - ci3
          tr2 = Cc(i-1,k,1) + cr3
          tr3 = Cc(i-1,k,1) - cr3
          Ch(i-1,1,k) = tr1 + tr2
          Ch(ic-1,4,k) = tr2 - tr1
          Ch(i,1,k) = ti1 + ti2
          Ch(ic,4,k) = ti1 - ti2
          Ch(i-1,3,k) = ti4 + tr3
          Ch(ic-1,2,k) = tr3 - ti4
          Ch(i,3,k) = tr4 + ti3
          Ch(ic,2,k) = tr4 - ti3
        END DO
      END DO
    END IF
    IF( MOD(Ido,2)==1 ) RETURN
  END IF
  DO k = 1, L1
    ti1 = -hsqt2*(Cc(Ido,k,2)+Cc(Ido,k,4))
    tr1 = hsqt2*(Cc(Ido,k,2)-Cc(Ido,k,4))
    Ch(Ido,1,k) = tr1 + Cc(Ido,k,1)
    Ch(Ido,3,k) = Cc(Ido,k,1) - tr1
    Ch(1,2,k) = ti1 - Cc(Ido,k,3)
    Ch(1,4,k) = ti1 + Cc(Ido,k,3)
  END DO
  RETURN
END SUBROUTINE RADF4
