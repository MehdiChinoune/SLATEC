!** PASSB2
SUBROUTINE PASSB2(Ido,L1,Cc,Ch,Wa1)
  IMPLICIT NONE
  !>
  !***
  !  Calculate the fast Fourier transform of subvectors of
  !            length two.
  !***
  ! **Library:**   SLATEC (FFTPACK)
  !***
  ! **Type:**      SINGLE PRECISION (PASSB2-S)
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
  
  INTEGER i, Ido, k, L1
  REAL Cc(Ido,2,*), Ch(Ido,L1,2), ti2, tr2, Wa1(*)
  !* FIRST EXECUTABLE STATEMENT  PASSB2
  IF ( Ido<=2 ) THEN
    DO k = 1, L1
      Ch(1,k,1) = Cc(1,1,k) + Cc(1,2,k)
      Ch(1,k,2) = Cc(1,1,k) - Cc(1,2,k)
      Ch(2,k,1) = Cc(2,1,k) + Cc(2,2,k)
      Ch(2,k,2) = Cc(2,1,k) - Cc(2,2,k)
    ENDDO
    RETURN
  ELSEIF ( Ido/2<L1 ) THEN
    DO i = 2, Ido, 2
      !DIR$ IVDEP
      DO k = 1, L1
        Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(i-1,2,k)
        tr2 = Cc(i-1,1,k) - Cc(i-1,2,k)
        Ch(i,k,1) = Cc(i,1,k) + Cc(i,2,k)
        ti2 = Cc(i,1,k) - Cc(i,2,k)
        Ch(i,k,2) = Wa1(i-1)*ti2 + Wa1(i)*tr2
        Ch(i-1,k,2) = Wa1(i-1)*tr2 - Wa1(i)*ti2
      ENDDO
    ENDDO
    RETURN
  ENDIF
  DO k = 1, L1
    !DIR$ IVDEP
    DO i = 2, Ido, 2
      Ch(i-1,k,1) = Cc(i-1,1,k) + Cc(i-1,2,k)
      tr2 = Cc(i-1,1,k) - Cc(i-1,2,k)
      Ch(i,k,1) = Cc(i,1,k) + Cc(i,2,k)
      ti2 = Cc(i,1,k) - Cc(i,2,k)
      Ch(i,k,2) = Wa1(i-1)*ti2 + Wa1(i)*tr2
      Ch(i-1,k,2) = Wa1(i-1)*tr2 - Wa1(i)*ti2
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE PASSB2
