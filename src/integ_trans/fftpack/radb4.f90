!** RADB4
SUBROUTINE RADB4(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3)
  !>
  !  Calculate the fast Fourier transform of subvectors of
  !            length four.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADB4-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variable SQRT2 by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Ido, L1
  REAL(SP) :: Cc(Ido,4,L1), Ch(Ido,L1,4), Wa1(Ido), Wa2(Ido), Wa3(Ido)
  INTEGER :: i, ic, idp2, k
  REAL(SP) :: ci2, ci3, ci4, cr2, cr3, cr4, sqrt2, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
  !* FIRST EXECUTABLE STATEMENT  RADB4
  sqrt2 = SQRT(2.)
  DO k = 1, L1
    tr1 = Cc(1,1,k) - Cc(Ido,4,k)
    tr2 = Cc(1,1,k) + Cc(Ido,4,k)
    tr3 = Cc(Ido,2,k) + Cc(Ido,2,k)
    tr4 = Cc(1,3,k) + Cc(1,3,k)
    Ch(1,k,1) = tr2 + tr3
    Ch(1,k,2) = tr1 - tr4
    Ch(1,k,3) = tr2 - tr3
    Ch(1,k,4) = tr1 + tr4
  END DO
  IF ( Ido<2 ) RETURN
  IF ( Ido/=2 ) THEN
    idp2 = Ido + 2
    IF ( (Ido-1)/2<L1 ) THEN
      DO i = 3, Ido, 2
        ic = idp2 - i
        DO k = 1, L1
          ti1 = Cc(i,1,k) + Cc(ic,4,k)
          ti2 = Cc(i,1,k) - Cc(ic,4,k)
          ti3 = Cc(i,3,k) - Cc(ic,2,k)
          tr4 = Cc(i,3,k) + Cc(ic,2,k)
          tr1 = Cc(i-1,1,k) - Cc(ic-1,4,k)
          tr2 = Cc(i-1,1,k) + Cc(ic-1,4,k)
          ti4 = Cc(i-1,3,k) - Cc(ic-1,2,k)
          tr3 = Cc(i-1,3,k) + Cc(ic-1,2,k)
          Ch(i-1,k,1) = tr2 + tr3
          cr3 = tr2 - tr3
          Ch(i,k,1) = ti2 + ti3
          ci3 = ti2 - ti3
          cr2 = tr1 - tr4
          cr4 = tr1 + tr4
          ci2 = ti1 + ti4
          ci4 = ti1 - ti4
          Ch(i-1,k,2) = Wa1(i-2)*cr2 - Wa1(i-1)*ci2
          Ch(i,k,2) = Wa1(i-2)*ci2 + Wa1(i-1)*cr2
          Ch(i-1,k,3) = Wa2(i-2)*cr3 - Wa2(i-1)*ci3
          Ch(i,k,3) = Wa2(i-2)*ci3 + Wa2(i-1)*cr3
          Ch(i-1,k,4) = Wa3(i-2)*cr4 - Wa3(i-1)*ci4
          Ch(i,k,4) = Wa3(i-2)*ci4 + Wa3(i-1)*cr4
        END DO
      END DO
    ELSE
      DO k = 1, L1
        DO i = 3, Ido, 2
          ic = idp2 - i
          ti1 = Cc(i,1,k) + Cc(ic,4,k)
          ti2 = Cc(i,1,k) - Cc(ic,4,k)
          ti3 = Cc(i,3,k) - Cc(ic,2,k)
          tr4 = Cc(i,3,k) + Cc(ic,2,k)
          tr1 = Cc(i-1,1,k) - Cc(ic-1,4,k)
          tr2 = Cc(i-1,1,k) + Cc(ic-1,4,k)
          ti4 = Cc(i-1,3,k) - Cc(ic-1,2,k)
          tr3 = Cc(i-1,3,k) + Cc(ic-1,2,k)
          Ch(i-1,k,1) = tr2 + tr3
          cr3 = tr2 - tr3
          Ch(i,k,1) = ti2 + ti3
          ci3 = ti2 - ti3
          cr2 = tr1 - tr4
          cr4 = tr1 + tr4
          ci2 = ti1 + ti4
          ci4 = ti1 - ti4
          Ch(i-1,k,2) = Wa1(i-2)*cr2 - Wa1(i-1)*ci2
          Ch(i,k,2) = Wa1(i-2)*ci2 + Wa1(i-1)*cr2
          Ch(i-1,k,3) = Wa2(i-2)*cr3 - Wa2(i-1)*ci3
          Ch(i,k,3) = Wa2(i-2)*ci3 + Wa2(i-1)*cr3
          Ch(i-1,k,4) = Wa3(i-2)*cr4 - Wa3(i-1)*ci4
          Ch(i,k,4) = Wa3(i-2)*ci4 + Wa3(i-1)*cr4
        END DO
      END DO
    END IF
    IF ( MOD(Ido,2)==1 ) RETURN
  END IF
  DO k = 1, L1
    ti1 = Cc(1,2,k) + Cc(1,4,k)
    ti2 = Cc(1,4,k) - Cc(1,2,k)
    tr1 = Cc(Ido,1,k) - Cc(Ido,3,k)
    tr2 = Cc(Ido,1,k) + Cc(Ido,3,k)
    Ch(Ido,k,1) = tr2 + tr2
    Ch(Ido,k,2) = sqrt2*(tr1-ti1)
    Ch(Ido,k,3) = ti2 + ti2
    Ch(Ido,k,4) = -sqrt2*(tr1+ti1)
  END DO
  RETURN
END SUBROUTINE RADB4
