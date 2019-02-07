!*==MPMUL2.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK MPMUL2
SUBROUTINE MPMUL2(X,Iy,Z,Trunc)
  IMPLICIT NONE
  !*--MPMUL25
  !*** Start of declarations inserted by SPAG
  INTEGER i, ij, is, ix, Iy, j, j1, j2, LUN, M, MXR
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  MPMUL2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPMUL2-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !  Multiplies 'mp' X by single-precision integer IY giving 'mp' Z.
  !  Multiplication by 1 may be used to normalize a number even if some
  !  digits are greater than B-1. Result is rounded if TRUNC.EQ.0,
  !  otherwise truncated.
  !
  !  The arguments X(*) and Z(*), and the variable R in COMMON are all
  !  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
  !  for the reason for this choice.
  !
  !***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
  !***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPOVFL, MPSTR
  !***COMMON BLOCKS    MPCOM
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  !***END PROLOGUE  MPMUL2
  COMMON /MPCOM / B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Z(*), Trunc, re, rs
  INTEGER c, c1, c2, ri, t1, t3, t4
  !***FIRST EXECUTABLE STATEMENT  MPMUL2
  rs = X(1)
  IF ( rs/=0 ) THEN
    j = Iy
    IF ( j<0 ) THEN
      j = -j
      rs = -rs
      ! CHECK FOR MULTIPLICATION BY B
      IF ( j/=B ) GOTO 200
      IF ( X(2)<M ) THEN
        CALL MPSTR(X,Z)
        Z(1) = rs
        Z(2) = X(2) + 1
        RETURN
      ELSE
        CALL MPCHK(1,4)
        WRITE (LUN,99001)
        99001       FORMAT (' *** OVERFLOW OCCURRED IN MPMUL2 ***')
        CALL MPOVFL(Z)
        RETURN
      ENDIF
    ELSEIF ( j/=0 ) THEN
      GOTO 200
    ENDIF
  ENDIF
  ! RESULT ZERO
  100  Z(1) = 0
  RETURN
  ! SET EXPONENT TO EXPONENT(X) + 4
  200  re = X(2) + 4
  ! FORM PRODUCT IN ACCUMULATOR
  c = 0
  t1 = T + 1
  t3 = T + 3
  t4 = T + 4
  ! IF J*B NOT REPRESENTABLE AS AN INTEGER WE HAVE TO SIMULATE
  ! DOUBLE-PRECISION MULTIPLICATION.
  IF ( j>=MAX(8*B,32767/B) ) THEN
    ! HERE J IS TOO LARGE FOR SINGLE-PRECISION MULTIPLICATION
    j1 = j/B
    j2 = j - j1*B
    ! FORM PRODUCT
    DO ij = 1, t4
      c1 = c/B
      c2 = c - B*c1
      i = t1 - ij
      ix = 0
      IF ( i>0 ) ix = X(i+2)
      ri = j2*ix + c2
      is = ri/B
      c = j1*ix + c1 + is
      R(i+4) = ri - B*is
    ENDDO
    IF ( c<0 ) GOTO 400
    IF ( c==0 ) GOTO 300
  ELSE
    DO ij = 1, T
      i = t1 - ij
      ri = j*X(i+2) + c
      c = ri/B
      R(i+4) = ri - B*c
    ENDDO
    ! CHECK FOR INTEGER OVERFLOW
    IF ( ri<0 ) GOTO 400
    ! HAVE TO TREAT FIRST FOUR WORDS OF R SEPARATELY
    DO ij = 1, 4
      i = 5 - ij
      ri = c
      c = ri/B
      R(i) = ri - B*c
    ENDDO
    IF ( c==0 ) GOTO 300
  ENDIF
  DO
    ! HAVE TO SHIFT RIGHT HERE AS CARRY OFF END
    DO ij = 1, t3
      i = t4 - ij
      R(i+1) = R(i)
    ENDDO
    ri = c
    c = ri/B
    R(1) = ri - B*c
    re = re + 1
    IF ( c<0 ) GOTO 400
    IF ( c==0 ) EXIT
  ENDDO
  ! NORMALIZE AND ROUND OR TRUNCATE RESULT
  300  CALL MPNZR(rs,re,Z,Trunc)
  RETURN
  ! CAN ONLY GET HERE IF INTEGER OVERFLOW OCCURRED
  400  CALL MPCHK(1,4)
  WRITE (LUN,99002)
  99002 FORMAT (' *** INTEGER OVERFLOW IN MPMUL2, B TOO LARGE ***')
  CALL MPERR
  GOTO 100
END SUBROUTINE MPMUL2
