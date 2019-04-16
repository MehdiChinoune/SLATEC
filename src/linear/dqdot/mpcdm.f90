!** MPCDM
SUBROUTINE MPCDM(Dx,Z)
  USE MPCOM
  !>
  !***
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPCDM-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  ! Converts double-precision number DX to multiple-precision Z.
  ! Some numbers will not convert exactly on machines with base
  ! other than two, four or sixteen. This routine is not called
  ! by any other routine in 'mp', so may be omitted if double-
  ! precision is not available.
  !
  ! The argument Z(*) and the variable R in COMMON are both INTEGER
  ! arrays of size 30.  See the comments in the routine MPBLAS for the
  ! for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPCHK, MPDIVI, MPMULI, MPNZR
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)

  INTEGER i, i2, ib, ie, k, Z(*), rs, re, tp
  REAL(8) :: db, dj, Dx
  !* FIRST EXECUTABLE STATEMENT  MPCDM
  CALL MPCHK(1,4)
  i2 = T + 4
  ! CHECK SIGN
  IF ( Dx<0 ) THEN
    ! DX .LT. 0D0
    rs = -1
    dj = -Dx
    ie = 0
  ELSEIF ( Dx==0 ) THEN
    ! IF DX = 0D0 RETURN 0
    Z(1) = 0
    RETURN
  ELSE
    ! DX .GT. 0D0
    rs = 1
    dj = Dx
    ie = 0
  END IF
  DO WHILE ( dj>=1D0 )
    ! INCREASE IE AND DIVIDE DJ BY 16.
    ie = ie + 1
    dj = 0.0625D0*dj
  END DO
  DO WHILE ( dj<0.0625D0 )
    ie = ie - 1
    dj = 16D0*dj
  END DO
  ! NOW DJ IS DY DIVIDED BY SUITABLE POWER OF 16
  ! SET EXPONENT TO 0
  re = 0
  db = REAL(B, 8)
  ! CONVERSION LOOP (ASSUME DOUBLE-PRECISION OPS. EXACT)
  DO i = 1, i2
    dj = db*dj
    R(i) = INT(dj)
    dj = dj - REAL(R(i), 8)
  END DO
  ! NORMALIZE RESULT
  CALL MPNZR(rs,re,Z,0)
  ib = MAX(7*B*B,32767)/16
  tp = 1
  ! NOW MULTIPLY BY 16**IE
  IF ( ie<0 ) THEN
    k = -ie
    DO i = 1, k
      tp = 16*tp
      IF ( (tp>ib).OR.(tp==B).OR.(i>=k) ) THEN
        CALL MPDIVI(Z,tp,Z)
        tp = 1
      END IF
    END DO
    RETURN
  ELSEIF ( ie/=0 ) THEN
    DO i = 1, ie
      tp = 16*tp
      IF ( (tp>ib).OR.(tp==B).OR.(i>=ie) ) THEN
        CALL MPMULI(Z,tp,Z)
        tp = 1
      END IF
    END DO
  END IF
END SUBROUTINE MPCDM
