!DECK TRIX
SUBROUTINE TRIX(Idegbr,Idegcr,M,A,B,C,Y,Tcos,D,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TRIX
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to GENBUN
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Subroutine to solve a system of linear equations where the
  !     coefficient matrix is a rational function in the matrix given by
  !     TRIDIAGONAL  ( . . ., A(I), B(I), C(I), . . . ).
  !
  !***SEE ALSO  GENBUN
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  TRIX
  REAL A, B, C, D, Tcos, W, x, xx, Y, z
  INTEGER i, Idegbr, Idegcr, ip, k, l, lint, M, mm1
  DIMENSION A(*), B(*), C(*), Y(*), Tcos(*), D(*), W(*)
  INTEGER kb, kc
  !***FIRST EXECUTABLE STATEMENT  TRIX
  mm1 = M - 1
  kb = Idegbr + 1
  kc = Idegcr + 1
  l = (Idegbr+1)/(Idegcr+1)
  lint = 1
  DO k = 1, Idegbr
    x = Tcos(k)
    IF ( k==l ) THEN
      i = Idegbr + lint
      xx = x - Tcos(i)
      DO i = 1, M
        W(i) = Y(i)
        Y(i) = xx*Y(i)
      ENDDO
    ENDIF
    z = 1./(B(1)-x)
    D(1) = C(1)*z
    Y(1) = Y(1)*z
    DO i = 2, mm1
      z = 1./(B(i)-x-A(i)*D(i-1))
      D(i) = C(i)*z
      Y(i) = (Y(i)-A(i)*Y(i-1))*z
    ENDDO
    z = B(M) - x - A(M)*D(mm1)
    IF ( z/=0. ) THEN
      Y(M) = (Y(M)-A(M)*Y(mm1))/z
    ELSE
      Y(M) = 0.
    ENDIF
    DO ip = 1, mm1
      i = M - ip
      Y(i) = Y(i) - D(i)*Y(i+1)
    ENDDO
    IF ( k==l ) THEN
      DO i = 1, M
        Y(i) = Y(i) + W(i)
      ENDDO
      lint = lint + 1
      l = (lint*kb)/kc
    ENDIF
  ENDDO
END SUBROUTINE TRIX