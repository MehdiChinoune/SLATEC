!** RADB2
SUBROUTINE RADB2(Ido,L1,Cc,Ch,Wa1)
  IMPLICIT NONE
  !>
  !***
  !  Calculate the fast Fourier transform of subvectors of
  !            length two.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (RADB2-S)
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
  
  REAL Cc, Ch, ti2, tr2, Wa1
  INTEGER i, ic, Ido, idp2, k, L1
  DIMENSION Cc(Ido,2,*), Ch(Ido,L1,2), Wa1(*)
  !* FIRST EXECUTABLE STATEMENT  RADB2
  DO k = 1, L1
    Ch(1,k,1) = Cc(1,1,k) + Cc(Ido,2,k)
    Ch(1,k,2) = Cc(1,1,k) - Cc(Ido,2,k)
  ENDDO
  IF ( Ido<2 ) RETURN
  IF ( Ido/=2 ) THEN
    idp2 = Ido + 2
    IF ( (Ido-1)/2<L1 ) THEN
      DO i = 3, Ido, 2
        ic = idp2 - i
        !DIR$ IVDEP
        DO k = 1, L1
          Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(ic-1,2,k)
          tr2 = Cc(i-1,1,k) - Cc(ic-1,2,k)
          Ch(i,k,1) = Cc(i,1,k) - Cc(ic,2,k)
          ti2 = Cc(i,1,k) + Cc(ic,2,k)
          Ch(i-1,k,2) = Wa1(i-2)*tr2 - Wa1(i-1)*ti2
          Ch(i,k,2) = Wa1(i-2)*ti2 + Wa1(i-1)*tr2
        ENDDO
      ENDDO
    ELSE
      DO k = 1, L1
        !DIR$ IVDEP
        DO i = 3, Ido, 2
          ic = idp2 - i
          Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(ic-1,2,k)
          tr2 = Cc(i-1,1,k) - Cc(ic-1,2,k)
          Ch(i,k,1) = Cc(i,1,k) - Cc(ic,2,k)
          ti2 = Cc(i,1,k) + Cc(ic,2,k)
          Ch(i-1,k,2) = Wa1(i-2)*tr2 - Wa1(i-1)*ti2
          Ch(i,k,2) = Wa1(i-2)*ti2 + Wa1(i-1)*tr2
        ENDDO
      ENDDO
    ENDIF
    IF ( MOD(Ido,2)==1 ) RETURN
  ENDIF
  DO k = 1, L1
    Ch(Ido,k,1) = Cc(Ido,1,k) + Cc(Ido,1,k)
    Ch(Ido,k,2) = -(Cc(1,2,k)+Cc(1,2,k))
  ENDDO
  RETURN
END SUBROUTINE RADB2
