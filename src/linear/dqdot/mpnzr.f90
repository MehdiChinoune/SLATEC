!** MPNZR
SUBROUTINE MPNZR(Rs,Re,Z,Trunc)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPNZR-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Modified for use with BLAS.  Blank COMMON changed to named COMMON.
  !  Assumes long (i.e. (t+4)-DIGIT) fraction in R, sign = RS, exponent
  !  = RE.  Normalizes, and returns 'mp' result in Z. Integer arguments
  !  RS and RE are not preserved. R*-rounding is used if TRUNC=0
  !
  !  The argument Z(*) and the variable R in COMMON are INTEGER arrays
  !  of size 30.  See the comments in the routine MPBLAS for the reason
  !  for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPERR, MPOVFL, MPUNFL
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  USE MPCOM, ONLY : b_com, lun_com, m_com, t_com, r_com
  INTEGER :: Re, Rs, Trunc, Z(30)
  INTEGER :: i, i2, i2m, i2p, is, it, j, k, b2
  !* FIRST EXECUTABLE STATEMENT  MPNZR
  i2 = t_com + 4
  IF( Rs/=0 ) THEN
    ! CHECK THAT SIGN = +-1
    IF( ABS(Rs)<=1 ) THEN
      ! LOOK FOR FIRST NONZERO DIGIT
      DO i = 1, i2
        is = i - 1
        IF( r_com(i)>0 ) GOTO 100
        ! FRACTION ZERO
      END DO
    ELSE
      WRITE (lun_com,99001)
      99001 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPNZR,',&
        ' POSSIBLE OVERWRITING PROBLEM ***')
      CALL MPERR
    END IF
  END IF
  ! STORE ZERO IN Z
  Z(1) = 0
  RETURN
  100 CONTINUE
  IF( is/=0 ) THEN
    ! NORMALIZE
    Re = Re - is
    i2m = i2 - is
    DO j = 1, i2m
      k = j + is
      r_com(j) = r_com(k)
    END DO
    i2p = i2m + 1
    DO j = i2p, i2
      r_com(j) = 0
    END DO
  END IF
  ! CHECK TO SEE IF TRUNCATION IS DESIRED
  IF( Trunc==0 ) THEN
    ! SEE IF ROUNDING NECESSARY
    ! TREAT EVEN AND ODD BASES DIFFERENTLY
    b2 = b_com/2
    IF( (2*b2)/=b_com ) THEN
      ! ODD BASE, ROUND IF R(T+1)... > 1/2
      DO i = 1, 4
        it = t_com + i
        IF( r_com(it)<b2 ) EXIT
        IF( r_com(it)/=b2 ) GOTO 150
      END DO
      GOTO 200
    ELSE
      ! B EVEN.  ROUND IF R(T+1)>=B2 UNLESS R(T) ODD AND ALL ZEROS
      ! AFTER R(T+2).
      IF( r_com(t_com+1)<b2 ) GOTO 200
      IF( r_com(t_com+1)==b2 ) THEN
        IF( MOD(r_com(t_com),2)/=0 ) THEN
          IF( (r_com(t_com+2)+r_com(t_com+3)+r_com(t_com+4))==0 ) GOTO 200
        END IF
      END IF
    END IF
    ! ROUND
    150 CONTINUE
    DO j = 1, t_com
      i = t_com + 1 - j
      r_com(i) = r_com(i) + 1
      IF( r_com(i)<b_com ) GOTO 200
      r_com(i) = 0
    END DO
    ! EXCEPTIONAL CASE, ROUNDED UP TO .10000...
    Re = Re + 1
    r_com(1) = 1
  END IF
  ! CHECK FOR OVERFLOW
  200 CONTINUE
  IF( Re>m_com ) THEN
    WRITE (lun_com,99002)
    99002 FORMAT (' *** OVERFLOW OCCURRED IN MPNZR ***')
    CALL MPOVFL(Z)
    RETURN
    ! CHECK FOR UNDERFLOW
  ELSEIF( Re<(-m_com) ) THEN
    ! UNDERFLOW HERE
    CALL MPUNFL(Z)
    RETURN
  END IF
  ! STORE RESULT IN Z
  Z(1) = Rs
  Z(2) = Re
  DO i = 1, t_com
    Z(i+2) = r_com(i)
  END DO
  RETURN
END SUBROUTINE MPNZR
