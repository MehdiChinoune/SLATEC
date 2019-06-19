!** MPDIVI
SUBROUTINE MPDIVI(X,Iy,Z)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPDIVI-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Divides 'mp' X by the single-precision integer IY giving 'mp' Z.
  !  This is much faster than division by an 'mp' number.
  !
  !  The arguments X(*) and Z(*), and the variable R in COMMON are all
  !  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
  !  for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPCHK, MPERR, MPNZR, MPSTR, MPUNFL
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
  USE MPCOM, ONLY : b_com, lun_com, m_com, t_com, r_com
  INTEGER :: Iy
  INTEGER :: X(30), Z(30)
  INTEGER :: i, i2, iq, iqj, ir, j, j1, j11, j2, k, kh, rs, re, r1, c, c2, b2
  !* FIRST EXECUTABLE STATEMENT  MPDIVI
  rs = X(1)
  j = Iy
  IF( j<0 ) THEN
    j = -j
    rs = -rs
  ELSEIF( j==0 ) THEN
    WRITE (lun_com,99001)
    99001 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MPDIVI ***')
    GOTO 500
  END IF
  re = X(2)
  ! CHECK FOR ZERO DIVIDEND
  IF( rs/=0 ) THEN
    ! CHECK FOR DIVISION BY B
    IF( j==b_com ) THEN
      CALL MPSTR(X,Z)
      IF( re<=(-m_com) ) THEN
        ! UNDERFLOW HERE
        CALL MPUNFL(Z)
        RETURN
      ELSE
        Z(1) = rs
        Z(2) = re - 1
        RETURN
      END IF
      ! CHECK FOR DIVISION BY 1 OR -1
    ELSEIF( j/=1 ) THEN
      c = 0
      i2 = t_com + 4
      i = 0
      ! IF J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE
      ! LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD.
      b2 = MAX(8*b_com,32767/b_com)
      IF( j>=b2 ) THEN
        ! HERE NEED SIMULATED DOUBLE-PRECISION DIVISION
        c2 = 0
        j1 = j/b_com
        j2 = j - j1*b_com
        j11 = j1 + 1
        DO
          ! LOOK FOR FIRST NONZERO DIGIT
          i = i + 1
          c = b_com*c + c2
          c2 = 0
          IF( i<=t_com ) c2 = X(i+2)
          IF( c<j1 ) THEN
          ELSEIF( c==j1 ) THEN
            IF( c2>=j2 ) GOTO 200
          ELSE
            GOTO 200
          END IF
        END DO
      ELSE
        DO
          ! LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT
          i = i + 1
          c = b_com*c
          IF( i<=t_com ) c = c + X(i+2)
          r1 = c/j
          IF( r1<0 ) GOTO 400
          IF( r1/=0 ) THEN
            ! ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT
            re = re + 1 - i
            r_com(1) = r1
            c = b_com*(c-j*r1)
            kh = 2
            IF( i<t_com ) THEN
              kh = 1 + t_com - i
              DO k = 2, kh
                i = i + 1
                c = c + X(i+2)
                r_com(k) = c/j
                c = b_com*(c-j*r_com(k))
              END DO
              IF( c<0 ) GOTO 400
              kh = kh + 1
            END IF
            DO k = kh, i2
              r_com(k) = c/j
              c = b_com*(c-j*r_com(k))
            END DO
            IF( c>=0 ) EXIT
            GOTO 400
          END IF
        END DO
      END IF
    ELSE
      CALL MPSTR(X,Z)
      Z(1) = rs
      RETURN
    END IF
  END IF
  ! NORMALIZE AND ROUND RESULT
  100  CALL MPNZR(rs,re,Z,0)
  RETURN
  ! COMPUTE T+4 QUOTIENT DIGITS
  200  re = re + 1 - i
  k = 1
  ! GET APPROXIMATE QUOTIENT FIRST
  300  ir = c/j11
  ! NOW REDUCE SO OVERFLOW DOES NOT OCCUR
  iq = c - ir*j1
  IF( iq>=b2 ) THEN
    ! HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR
    ir = ir + 1
    iq = iq - j1
  END IF
  iq = iq*b_com - ir*j2
  IF( iq<0 ) THEN
    ! HERE IQ NEGATIVE SO IR WAS TOO LARGE
    ir = ir - 1
    iq = iq + j
  END IF
  IF( i<=t_com ) iq = iq + X(i+2)
  iqj = iq/j
  ! R(K) = QUOTIENT, C = REMAINDER
  r_com(k) = iqj + ir
  c = iq - j*iqj
  IF( c>=0 ) THEN
    ! MAIN LOOP FOR LARGE ABS(IY) CASE
    k = k + 1
    IF( k>i2 ) GOTO 100
    i = i + 1
    GOTO 300
  END IF
  ! CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED
  400  CALL MPCHK(1,4)
  WRITE (lun_com,99002)
  99002 FORMAT (' *** INTEGER OVERFLOW IN MPDIVI, b_com TOO LARGE ***')
  500  CALL MPERR
  Z(1) = 0
  RETURN
END SUBROUTINE MPDIVI
