!** MPMUL
SUBROUTINE MPMUL(X,Y,Z)
  USE MPCOM
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPMUL-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Multiplies X and Y, returning result in Z, for 'mp' X, Y and Z.
  !  The simple o(t**2) algorithm is used, with four guard digits and
  !  R*-rounding. Advantage is taken of zero digits in X, but not in Y.
  !  Asymptotically faster algorithms are known (see Knuth, VOL. 2),
  !  but are difficult to implement in FORTRAN in an efficient and
  !  machine-independent manner. In comments to other 'mp' routines,
  !  M(t) is the time to perform t-digit 'mp' multiplication. Thus
  !  M(t) = o(t**2) with the present version of MPMUL, but
  !  M(t) = o(t.log(t).log(log(t))) is theoretically possible.
  !
  !  The arguments X(*), Y(*), and Z(*), and the variable R in COMMON are
  !  all INTEGER arrays of size 30.  See the comments in the routine
  !  MPBLAS for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPCHK, MPERR, MPMLP, MPNZR
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)

  INTEGER i, i2, i2p, j, j1, X(*), Y(*), Z(*), rs, re, xi, c, ri
  !* FIRST EXECUTABLE STATEMENT  MPMUL
  CALL MPCHK(1,4)
  i2 = T + 4
  i2p = i2 + 1
  ! FORM SIGN OF PRODUCT
  rs = X(1)*Y(1)
  IF ( rs/=0 ) THEN
    ! FORM EXPONENT OF PRODUCT
    re = X(2) + Y(2)
    ! CLEAR ACCUMULATOR
    DO i = 1, i2
      R(i) = 0
    END DO
    ! PERFORM MULTIPLICATION
    c = 8
    DO i = 1, T
      xi = X(i+2)
      ! FOR SPEED, PUT THE NUMBER WITH MANY ZEROS FIRST
      IF ( xi/=0 ) THEN
        CALL MPMLP(R(i+1),Y(3),xi,MIN(T,i2-i))
        c = c - 1
        IF ( c<=0 ) THEN
          ! CHECK FOR LEGAL BASE B DIGIT
          IF ( (xi<0).OR.(xi>=B) ) GOTO 200
          ! PROPAGATE CARRIES AT END AND EVERY EIGHTH TIME,
          ! FASTER THAN DOING IT EVERY TIME.
          DO j = 1, i2
            j1 = i2p - j
            ri = R(j1) + c
            IF ( ri<0 ) GOTO 100
            c = ri/B
            R(j1) = ri - B*c
          END DO
          IF ( c/=0 ) GOTO 200
          c = 8
        END IF
      END IF
    END DO
    IF ( c/=8 ) THEN
      IF ( (xi<0).OR.(xi>=B) ) GOTO 200
      c = 0
      DO j = 1, i2
        j1 = i2p - j
        ri = R(j1) + c
        IF ( ri<0 ) GOTO 100
        c = ri/B
        R(j1) = ri - B*c
      END DO
      IF ( c/=0 ) GOTO 200
    END IF
    ! NORMALIZE AND ROUND RESULT
    CALL MPNZR(rs,re,Z,0)
    RETURN
  ELSE
    ! SET RESULT TO ZERO
    Z(1) = 0
    RETURN
  END IF
  100  WRITE (LUN,99001)
  99001 FORMAT (' *** INTEGER OVERFLOW IN MPMUL, B TOO LARGE ***')
  GOTO 300
  200  WRITE (LUN,99002)
  99002 FORMAT (' *** ILLEGAL BASE B DIGIT IN CALL TO MPMUL,',&
    ' POSSIBLE OVERWRITING PROBLEM ***')
  300  CALL MPERR
  Z(1) = 0
END SUBROUTINE MPMUL
