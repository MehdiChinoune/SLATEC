!** DX
SUBROUTINE DX(U,Idmn,I,J,Uxxx,Uxxxx)
  !> Subsidiary to SEPELI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (DX-S)
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
  ! **See also:**  SEPELI
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPLPCM

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPLPCM, ONLY : dlx4_com, kswx_com, tdlx3_com, k_com, l_com
  INTEGER :: I, Idmn, J
  REAL(SP) :: U(Idmn,l_com), Uxxx, Uxxxx
  !* FIRST EXECUTABLE STATEMENT  DX
  IF( I>2 .AND. I<(k_com-1) ) THEN
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
    !
    Uxxx = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/tdlx3_com
    Uxxxx = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/dlx4_com
    RETURN
  ELSE
    IF( I/=1 ) THEN
      IF( I==2 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+dlx_com
        !
        IF( kswx_com==1 ) THEN
          !
          !     PERIODIC AT X=A+DLX
          !
          Uxxx = (-U(k_com-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(tdlx3_com)
          Uxxxx = (U(k_com-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/dlx4_com
          RETURN
        ELSE
          Uxxx = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))/tdlx3_com
          Uxxxx = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)&
            -U(6,J))/dlx4_com
          RETURN
        END IF
      ELSEIF( I==k_com-1 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
        !
        IF( kswx_com==1 ) THEN
          !
          !     PERIODIC AT X=B-DLX
          !
          Uxxx = (-U(k_com-3,J)+2.0*U(k_com-2,J)-2.0*U(1,J)+U(2,J))/tdlx3_com
          Uxxxx = (U(k_com-3,J)-4.0*U(k_com-2,J)+6.0*U(k_com-1,J)-4.0*U(1,J)+U(2,J))/dlx4_com
          RETURN
        ELSE
          Uxxx = (U(k_com-4,J)-6.0*U(k_com-3,J)+12.0*U(k_com-2,J)-10.0*U(k_com-1,J)&
            +3.0*U(k_com,J))/tdlx3_com
          Uxxxx = (-U(k_com-5,J)+6.0*U(k_com-4,J)-14.0*U(k_com-3,J)+16.0*U(k_com-2,J)&
            -9.0*U(k_com-1,J)+2.0*U(k_com,J))/dlx4_com
          RETURN
        END IF
      ELSEIF( I==k_com ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
        !
        Uxxx = -(3.0*U(k_com-4,J)-14.0*U(k_com-3,J)+24.0*U(k_com-2,J)-18.0*U(k_com-1,J)&
          +5.0*U(k_com,J))/tdlx3_com
        Uxxxx = (-2.0*U(k_com-5,J)+11.0*U(k_com-4,J)-24.0*U(k_com-3,J)+26.0*U(k_com-2,J)&
          -14.0*U(k_com-1,J)+3.0*U(k_com,J))/dlx4_com
        RETURN
      END IF
    END IF
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
    !
    IF( kswx_com/=1 ) THEN
      Uxxx = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-3.0*U(5,J))/tdlx3_com
      Uxxxx = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+11.0*U(5,J)&
        -2.0*U(6,J))/dlx4_com
      RETURN
    END IF
  END IF
  !
  !     PERIODIC AT X=A
  !
  Uxxx = (-U(k_com-2,J)+2.0*U(k_com-1,J)-2.0*U(2,J)+U(3,J))/(tdlx3_com)
  Uxxxx = (U(k_com-2,J)-4.0*U(k_com-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/dlx4_com
  RETURN
END SUBROUTINE DX
