!** DX4
SUBROUTINE DX4(U,Idmn,I,J,Uxxx,Uxxxx)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (DX4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This program computes second order finite difference
  !     approximations to the third and fourth X
  !     partial derivatives of U at the (I,J) mesh point.
  !
  !***
  ! **See also:**  SEPX4
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPL4

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL AIT, BIT, CIT, DIT, DLX, DLX4, DLY, DLY4, TDLx3, TDLy3, U, Uxxx, Uxxxx
  INTEGER I, Idmn, IS, J, JS, K, KSWx, KSWy, L, MIT, MS, NIT, NS
  COMMON /SPL4  / KSWx, KSWy, K, L, AIT, BIT, CIT, DIT, MIT, NIT, &
    IS, MS, JS, NS, DLX, DLY, TDLx3, TDLy3, DLX4, DLY4
  DIMENSION U(Idmn,*)
  !* FIRST EXECUTABLE STATEMENT  DX4
  IF ( I>2.AND.I<(K-1) ) THEN
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
    !
    Uxxx = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLx3
    Uxxxx = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/DLX4
    RETURN
  ELSE
    IF ( I/=1 ) THEN
      IF ( I==2 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
        !
        IF ( KSWx==1 ) THEN
          !
          !     PERIODIC AT X=A+DLX
          !
          Uxxx = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(TDLx3)
          Uxxxx = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/DLX4
          RETURN
        ELSE
          Uxxx = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))&
            /TDLx3
          Uxxxx = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)&
            -U(6,J))/DLX4
          RETURN
        ENDIF
      ELSEIF ( I==K-1 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
        !
        IF ( KSWx==1 ) THEN
          !
          !     PERIODIC AT X=B-DLX
          !
          Uxxx = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLx3
          Uxxxx = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J))&
            /DLX4
          RETURN
        ELSE
          Uxxx = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)&
            +3.0*U(K,J))/TDLx3
          Uxxxx = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)&
            -9.0*U(K-1,J)+2.0*U(K,J))/DLX4
          RETURN
        ENDIF
      ELSEIF ( I==K ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
        !
        Uxxx = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)&
          +5.0*U(K,J))/TDLx3
        Uxxxx = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)&
          -14.0*U(K-1,J)+3.0*U(K,J))/DLX4
        RETURN
      ENDIF
    ENDIF
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
    !
    IF ( KSWx/=1 ) THEN
      Uxxx = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-3.0*U(5,J))&
        /(TDLx3)
      Uxxxx = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+11.0*U(5,J)&
        -2.0*U(6,J))/DLX4
      RETURN
    ENDIF
  ENDIF
  !
  !     PERIODIC AT X=A
  !
  Uxxx = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/(TDLx3)
  Uxxxx = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/DLX4
  RETURN
END SUBROUTINE DX4
