!DECK PASSF5
SUBROUTINE PASSF5(Ido,L1,Cc,Ch,Wa1,Wa2,Wa3,Wa4)
  IMPLICIT NONE
  REAL Cc, Ch, ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, di3, &
    di4, di5, dr2, dr3, dr4, dr5, pi, ti11
  REAL ti12, ti2, ti3, ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5, &
    Wa1, Wa2, Wa3, Wa4
  INTEGER i, Ido, k, L1
  !***BEGIN PROLOGUE  PASSF5
  !***SUBSIDIARY
  !***PURPOSE  Calculate the fast Fourier transform of subvectors of
  !            length five.
  !***LIBRARY   SLATEC (FFTPACK)
  !***TYPE      SINGLE PRECISION (PASSF5-S)
  !***AUTHOR  Swarztrauber, P. N., (NCAR)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   830401  Modified to use SLATEC library source file format.
  !   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
  !           (a) changing dummy array size declarations (1) to (*),
  !           (b) changing definition of variables PI, TI11, TI12,
  !               TR11, TR12 by using FORTRAN intrinsic functions ATAN
  !               and SIN instead of DATA statements.
  !   881128  Modified by Dick Valent to meet prologue standards.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PASSF5
  DIMENSION Cc(Ido,5,*), Ch(Ido,L1,5), Wa1(*), Wa2(*), Wa3(*), Wa4(*)
  !***FIRST EXECUTABLE STATEMENT  PASSF5
  pi = 4.*ATAN(1.)
  tr11 = SIN(.1*pi)
  ti11 = -SIN(.4*pi)
  tr12 = -SIN(.3*pi)
  ti12 = -SIN(.2*pi)
  IF ( Ido==2 ) THEN
    DO k = 1, L1
      ti5 = Cc(2,2,k) - Cc(2,5,k)
      ti2 = Cc(2,2,k) + Cc(2,5,k)
      ti4 = Cc(2,3,k) - Cc(2,4,k)
      ti3 = Cc(2,3,k) + Cc(2,4,k)
      tr5 = Cc(1,2,k) - Cc(1,5,k)
      tr2 = Cc(1,2,k) + Cc(1,5,k)
      tr4 = Cc(1,3,k) - Cc(1,4,k)
      tr3 = Cc(1,3,k) + Cc(1,4,k)
      Ch(1,k,1) = Cc(1,1,k) + tr2 + tr3
      Ch(2,k,1) = Cc(2,1,k) + ti2 + ti3
      cr2 = Cc(1,1,k) + tr11*tr2 + tr12*tr3
      ci2 = Cc(2,1,k) + tr11*ti2 + tr12*ti3
      cr3 = Cc(1,1,k) + tr12*tr2 + tr11*tr3
      ci3 = Cc(2,1,k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      Ch(1,k,2) = cr2 - ci5
      Ch(1,k,5) = cr2 + ci5
      Ch(2,k,2) = ci2 + cr5
      Ch(2,k,3) = ci3 + cr4
      Ch(1,k,3) = cr3 - ci4
      Ch(1,k,4) = cr3 + ci4
      Ch(2,k,4) = ci3 - cr4
      Ch(2,k,5) = ci2 - cr5
    ENDDO
    RETURN
  ELSEIF ( Ido/2<L1 ) THEN
    DO i = 2, Ido, 2
      !DIR$ IVDEP
      DO k = 1, L1
        ti5 = Cc(i,2,k) - Cc(i,5,k)
        ti2 = Cc(i,2,k) + Cc(i,5,k)
        ti4 = Cc(i,3,k) - Cc(i,4,k)
        ti3 = Cc(i,3,k) + Cc(i,4,k)
        tr5 = Cc(i-1,2,k) - Cc(i-1,5,k)
        tr2 = Cc(i-1,2,k) + Cc(i-1,5,k)
        tr4 = Cc(i-1,3,k) - Cc(i-1,4,k)
        tr3 = Cc(i-1,3,k) + Cc(i-1,4,k)
        Ch(i-1,k,1) = Cc(i-1,1,k) + tr2 + tr3
        Ch(i,k,1) = Cc(i,1,k) + ti2 + ti3
        cr2 = Cc(i-1,1,k) + tr11*tr2 + tr12*tr3
        ci2 = Cc(i,1,k) + tr11*ti2 + tr12*ti3
        cr3 = Cc(i-1,1,k) + tr12*tr2 + tr11*tr3
        ci3 = Cc(i,1,k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
        Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
        Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
        Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
        Ch(i-1,k,4) = Wa3(i-1)*dr4 + Wa3(i)*di4
        Ch(i,k,4) = Wa3(i-1)*di4 - Wa3(i)*dr4
        Ch(i-1,k,5) = Wa4(i-1)*dr5 + Wa4(i)*di5
        Ch(i,k,5) = Wa4(i-1)*di5 - Wa4(i)*dr5
      ENDDO
    ENDDO
    GOTO 99999
  ENDIF
  DO k = 1, L1
    !DIR$ IVDEP
    DO i = 2, Ido, 2
      ti5 = Cc(i,2,k) - Cc(i,5,k)
      ti2 = Cc(i,2,k) + Cc(i,5,k)
      ti4 = Cc(i,3,k) - Cc(i,4,k)
      ti3 = Cc(i,3,k) + Cc(i,4,k)
      tr5 = Cc(i-1,2,k) - Cc(i-1,5,k)
      tr2 = Cc(i-1,2,k) + Cc(i-1,5,k)
      tr4 = Cc(i-1,3,k) - Cc(i-1,4,k)
      tr3 = Cc(i-1,3,k) + Cc(i-1,4,k)
      Ch(i-1,k,1) = Cc(i-1,1,k) + tr2 + tr3
      Ch(i,k,1) = Cc(i,1,k) + ti2 + ti3
      cr2 = Cc(i-1,1,k) + tr11*tr2 + tr12*tr3
      ci2 = Cc(i,1,k) + tr11*ti2 + tr12*ti3
      cr3 = Cc(i-1,1,k) + tr12*tr2 + tr11*tr3
      ci3 = Cc(i,1,k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5
      Ch(i-1,k,2) = Wa1(i-1)*dr2 + Wa1(i)*di2
      Ch(i,k,2) = Wa1(i-1)*di2 - Wa1(i)*dr2
      Ch(i-1,k,3) = Wa2(i-1)*dr3 + Wa2(i)*di3
      Ch(i,k,3) = Wa2(i-1)*di3 - Wa2(i)*dr3
      Ch(i-1,k,4) = Wa3(i-1)*dr4 + Wa3(i)*di4
      Ch(i,k,4) = Wa3(i-1)*di4 - Wa3(i)*dr4
      Ch(i-1,k,5) = Wa4(i-1)*dr5 + Wa4(i)*di5
      Ch(i,k,5) = Wa4(i-1)*di5 - Wa4(i)*dr5
    ENDDO
  ENDDO
  RETURN
  99999 CONTINUE
  END SUBROUTINE PASSF5
