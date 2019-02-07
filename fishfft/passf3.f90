!*==PASSF3.f90  processed by SPAG 6.72Dc at 10:55 on  6 Feb 2019
!DECK PASSF3
SUBROUTINE PASSF3(Ido,L1,Cc,Ch,Wa1,Wa2)
  IMPLICIT NONE
  !*--PASSF35
  !*** Start of declarations inserted by SPAG
  REAL Cc, Ch, ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, &
    taur, ti2, tr2, Wa1, Wa2
  INTEGER i, Ido, k, L1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PASSF3
  !***SUBSIDIARY
  !***PURPOSE  Calculate the fast Fourier transform of subvectors of
  !            length three.
  !***LIBRARY   SLATEC (FFTPACK)
  !***TYPE      SINGLE PRECISION (PASSF3-S)
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
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
  !***END PROLOGUE  PASSF3
  DIMENSION Cc(Ido,3,*), Ch(Ido,L1,3), Wa1(*), Wa2(*)
  !***FIRST EXECUTABLE STATEMENT  PASSF3
  taur = -.5
  taui = -.5*SQRT(3.)
  IF ( Ido==2 ) THEN
    DO k = 1, L1
      tr2 = Cc(1,2,k) + Cc(1,3,k)
      cr2 = Cc(1,1,k) + taur*tr2
      Ch(1,k,1) = Cc(1,1,k) + tr2
      ti2 = Cc(2,2,k) + Cc(2,3,k)
      ci2 = Cc(2,1,k) + taur*ti2
      Ch(2,k,1) = Cc(2,1,k) + ti2
      cr3 = taui*(Cc(1,2,k)-Cc(1,3,k))
      ci3 = taui*(Cc(2,2,k)-Cc(2,3,k))
      Ch(1,k,2) = cr2 - ci3
      Ch(1,k,3) = cr2 + ci3
      Ch(2,k,2) = ci2 + cr3
      Ch(2,k,3) = ci2 - cr3
    ENDDO
    RETURN
  ELSEIF ( Ido/2<L1 ) THEN
    DO i = 2, Ido, 2
      !DIR$ IVDEP
      DO k = 1, L1
        tr2 = Cc(i-1,2,k) + Cc(i-1,3,k)
        cr2 = Cc(i-1,1,k) + taur*tr2
        Ch(i-1,k,1) = Cc(i-1,1,k) + tr2
        ti2 = Cc(i,2,k) + Cc(i,3,k)
        ci2 = Cc(i,1,k) + taur*ti2
        Ch(i,k,1) = Cc(i,1,k) + ti2
        cr3 = taui*(Cc(i-1,2,k)-Cc(i-1,3,k))
        ci3 = taui*(Cc(i,2,k)-Cc(i,3,k))
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
        Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
        Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
        Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
      ENDDO
    ENDDO
    GOTO 99999
  ENDIF
  DO k = 1, L1
    !DIR$ IVDEP
    DO i = 2, Ido, 2
      tr2 = Cc(i-1,2,k) + Cc(i-1,3,k)
      cr2 = Cc(i-1,1,k) + taur*tr2
      Ch(i-1,k,1) = Cc(i-1,1,k) + tr2
      ti2 = Cc(i,2,k) + Cc(i,3,k)
      ci2 = Cc(i,1,k) + taur*ti2
      Ch(i,k,1) = Cc(i,1,k) + ti2
      cr3 = taui*(Cc(i-1,2,k)-Cc(i-1,3,k))
      ci3 = taui*(Cc(i,2,k)-Cc(i,3,k))
      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3
      Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
      Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
      Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
      Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
    ENDDO
  ENDDO
  RETURN
  99999 CONTINUE
  END SUBROUTINE PASSF3
