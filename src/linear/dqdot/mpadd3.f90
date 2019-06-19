!** MPADD3
SUBROUTINE MPADD3(X,Y,S,Med,Re)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPADD3-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   Called by MPADD2; does inner loops of addition
  !
  !   The arguments X(*) and Y(*) and the variable R in COMMON are all
  !   INTEGER arrays of size 30.  See the comments in the routine MPBLAS
  !   for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  USE MPCOM, ONLY : b_com, t_com, r_com
  INTEGER :: Med, S, Re
  INTEGER :: X(30), Y(30)
  INTEGER :: i, i2, i2p, j, c, ted
  !* FIRST EXECUTABLE STATEMENT  MPADD3
  ted = t_com + Med
  i2 = t_com + 4
  i = i2
  c = 0
  ! CLEAR GUARD DIGITS TO RIGHT OF X DIGITS
  DO WHILE( i>ted )
    r_com(i) = 0
    i = i - 1
  END DO
  IF( S<0 ) THEN
    DO WHILE( i>t_com )
      ! HERE DO SUBTRACTION, ABS(Y) > ABS(X)
      j = i - Med
      r_com(i) = c - X(j+2)
      c = 0
      IF( r_com(i)<0 ) THEN
        ! BORROW GENERATED HERE
        c = -1
        r_com(i) = r_com(i) + b_com
      END IF
      i = i - 1
    END DO
    DO WHILE( i>Med )
      j = i - Med
      c = Y(i+2) + c - X(j+2)
      IF( c>=0 ) THEN
        ! NO BORROW GENERATED HERE
        r_com(i) = c
        c = 0
        i = i - 1
      ELSE
        ! BORROW GENERATED HERE
        r_com(i) = c + b_com
        c = -1
        i = i - 1
      END IF
    END DO
    DO
      IF( i<=0 ) RETURN
      c = Y(i+2) + c
      IF( c>=0 ) EXIT
      r_com(i) = c + b_com
      c = -1
      i = i - 1
    END DO
  ELSE
    ! HERE DO ADDITION, EXPONENT(Y) >= EXPONENT(X)
    IF( i>=t_com ) THEN
      DO
        j = i - Med
        r_com(i) = X(j+2)
        i = i - 1
        IF( i<=t_com ) EXIT
      END DO
    END IF
    DO WHILE( i>Med )
      j = i - Med
      c = Y(i+2) + X(j+2) + c
      IF( c<b_com ) THEN
        ! NO CARRY GENERATED HERE
        r_com(i) = c
        c = 0
        i = i - 1
      ELSE
        ! CARRY GENERATED HERE
        r_com(i) = c - b_com
        c = 1
        i = i - 1
      END IF
    END DO
    DO WHILE( i>0 )
      c = Y(i+2) + c
      IF( c<b_com ) GOTO 100
      r_com(i) = 0
      c = 1
      i = i - 1
    END DO
    IF( c==0 ) RETURN
    ! MUST SHIFT RIGHT HERE AS CARRY OFF END
    i2p = i2 + 1
    DO j = 2, i2
      i = i2p - j
      r_com(i+1) = r_com(i)
    END DO
    r_com(1) = 1
    Re = Re + 1
    RETURN
  END IF
  100  r_com(i) = c
  i = i - 1
  DO
    ! NO CARRY POSSIBLE HERE
    IF( i<=0 ) RETURN
    r_com(i) = Y(i+2)
    i = i - 1
  END DO
END SUBROUTINE MPADD3
