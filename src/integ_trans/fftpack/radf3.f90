!** RADF3
SUBROUTINE RADF3(Ido,L1,Cc,Ch,Wa1,Wa2)
  !>
  !  Calculate the fast Fourier transform of subvectors of
  !            length three.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADF3-S)
  !***
  ! **Author:**  Swarztrauber, P. N., (NCAR)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variable TAUI by using
  !               FORTRAN intrinsic function SQRT instead of a DATA
  !               statement.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER i, ic, Ido, idp2, k, L1
  REAL Cc(Ido,L1,3), Ch(Ido,3,*), ci2, cr2, di2, di3, dr2, dr3, taui, taur, ti2, &
    ti3, tr2, tr3, Wa1(*), Wa2(*)
  !* FIRST EXECUTABLE STATEMENT  RADF3
  taur = -.5
  taui = .5*SQRT(3.)
  DO k = 1, L1
    cr2 = Cc(1,k,2) + Cc(1,k,3)
    Ch(1,1,k) = Cc(1,k,1) + cr2
    Ch(1,3,k) = taui*(Cc(1,k,3)-Cc(1,k,2))
    Ch(Ido,2,k) = Cc(1,k,1) + taur*cr2
  END DO
  IF ( Ido==1 ) RETURN
  idp2 = Ido + 2
  IF ( (Ido-1)/2<L1 ) THEN
    DO i = 3, Ido, 2
      ic = idp2 - i
      !DIR$ IVDEP
      DO k = 1, L1
        dr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
        di2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
        dr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
        di3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
        cr2 = dr2 + dr3
        ci2 = di2 + di3
        Ch(i-1,1,k) = Cc(i-1,k,1) + cr2
        Ch(i,1,k) = Cc(i,k,1) + ci2
        tr2 = Cc(i-1,k,1) + taur*cr2
        ti2 = Cc(i,k,1) + taur*ci2
        tr3 = taui*(di2-di3)
        ti3 = taui*(dr3-dr2)
        Ch(i-1,3,k) = tr2 + tr3
        Ch(ic-1,2,k) = tr2 - tr3
        Ch(i,3,k) = ti2 + ti3
        Ch(ic,2,k) = ti3 - ti2
      END DO
    END DO
    RETURN
  END IF
  DO k = 1, L1
    !DIR$ IVDEP
    DO i = 3, Ido, 2
      ic = idp2 - i
      dr2 = Wa1(i-2)*Cc(i-1,k,2) + Wa1(i-1)*Cc(i,k,2)
      di2 = Wa1(i-2)*Cc(i,k,2) - Wa1(i-1)*Cc(i-1,k,2)
      dr3 = Wa2(i-2)*Cc(i-1,k,3) + Wa2(i-1)*Cc(i,k,3)
      di3 = Wa2(i-2)*Cc(i,k,3) - Wa2(i-1)*Cc(i-1,k,3)
      cr2 = dr2 + dr3
      ci2 = di2 + di3
      Ch(i-1,1,k) = Cc(i-1,k,1) + cr2
      Ch(i,1,k) = Cc(i,k,1) + ci2
      tr2 = Cc(i-1,k,1) + taur*cr2
      ti2 = Cc(i,k,1) + taur*ci2
      tr3 = taui*(di2-di3)
      ti3 = taui*(dr3-dr2)
      Ch(i-1,3,k) = tr2 + tr3
      Ch(ic-1,2,k) = tr2 - tr3
      Ch(i,3,k) = ti2 + ti3
      Ch(ic,2,k) = ti3 - ti2
    END DO
  END DO
  RETURN
END SUBROUTINE RADF3
