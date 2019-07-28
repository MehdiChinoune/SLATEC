!** DY4
PURE SUBROUTINE DY4(U,Idmn,I,J,Uyyy,Uyyyy)
  !> Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (DY4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This program computes second order finite difference
  !     approximations to the third and fourth Y
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
  USE SPL4, ONLY : l_com, dly4_com, kswy_com, tdly3_com
  !
  INTEGER, INTENT(IN) :: I, Idmn, J
  REAL(SP), INTENT(IN) :: U(Idmn,l_com)
  REAL(SP), INTENT(OUT) :: Uyyy, Uyyyy
  !* FIRST EXECUTABLE STATEMENT  DY4
  IF( J>2 .AND. J<(l_com-1) ) THEN
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
    !
    Uyyy = (-U(I,J-2)+2._SP*U(I,J-1)-2._SP*U(I,J+1)+U(I,J+2))/tdly3_com
    Uyyyy = (U(I,J-2)-4._SP*U(I,J-1)+6._SP*U(I,J)-4._SP*U(I,J+1)+U(I,J+2))/dly4_com
    RETURN
  ELSE
    IF( J/=1 ) THEN
      IF( J==2 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
        !
        IF( kswy_com==1 ) THEN
          !
          !     PERIODIC AT Y=C+DLY
          !
          Uyyy = (-U(I,l_com-1)+2._SP*U(I,1)-2._SP*U(I,3)+U(I,4))/tdly3_com
          Uyyyy = (U(I,l_com-1)-4._SP*U(I,1)+6._SP*U(I,2)-4._SP*U(I,3)+U(I,4))/dly4_com
          RETURN
        ELSE
          Uyyy = (-3._SP*U(I,1)+10._SP*U(I,2)-12._SP*U(I,3)+6._SP*U(I,4)-U(I,5))/tdly3_com
          Uyyyy = (2._SP*U(I,1)-9._SP*U(I,2)+16._SP*U(I,3)-14._SP*U(I,4)+6._SP*U(I,5)&
            -U(I,6))/dly4_com
          RETURN
        END IF
      ELSEIF( J==l_com-1 ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
        !
        IF( kswy_com==1 ) THEN
          !
          !     PERIODIC AT Y=D-DLY
          !
          Uyyy = (-U(I,l_com-3)+2._SP*U(I,l_com-2)-2._SP*U(I,1)+U(I,2))/tdly3_com
          Uyyyy = (U(I,l_com-3)-4._SP*U(I,l_com-2)+6._SP*U(I,l_com-1)&
            -4._SP*U(I,1)+U(I,2))/dly4_com
          RETURN
        ELSE
          Uyyy = (U(I,l_com-4)-6._SP*U(I,l_com-3)+12._SP*U(I,l_com-2)-10._SP*U(I,l_com-1)&
            +3._SP*U(I,l_com))/tdly3_com
          Uyyyy = (-U(I,l_com-5)+6._SP*U(I,l_com-4)-14._SP*U(I,l_com-3)+16._SP*U(I,l_com-2)&
            -9._SP*U(I,l_com-1)+2._SP*U(I,l_com))/dly4_com
          RETURN
        END IF
      ELSEIF( J==l_com ) THEN
        !
        !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
        !
        Uyyy = -(3._SP*U(I,l_com-4)-14._SP*U(I,l_com-3)+24._SP*U(I,l_com-2)&
          -18._SP*U(I,l_com-1)+5._SP*U(I,l_com))/tdly3_com
        Uyyyy = (-2._SP*U(I,l_com-5)+11._SP*U(I,l_com-4)-24._SP*U(I,l_com-3)&
          +26._SP*U(I,l_com-2)-14._SP*U(I,l_com-1)+3._SP*U(I,l_com))/dly4_com
        RETURN
      END IF
    END IF
    !
    !     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
    !
    IF( kswy_com/=1 ) THEN
      Uyyy = (-5._SP*U(I,1)+18._SP*U(I,2)-24._SP*U(I,3)+14._SP*U(I,4)-3._SP*U(I,5))/tdly3_com
      Uyyyy = (3._SP*U(I,1)-14._SP*U(I,2)+26._SP*U(I,3)-24._SP*U(I,4)+11._SP*U(I,5)&
        -2._SP*U(I,6))/dly4_com
      RETURN
    END IF
  END IF
  !
  !     PERIODIC AT X=A
  !
  Uyyy = (-U(I,l_com-2)+2._SP*U(I,l_com-1)-2._SP*U(I,2)+U(I,3))/tdly3_com
  Uyyyy = (U(I,l_com-2)-4._SP*U(I,l_com-1)+6._SP*U(I,1)-4._SP*U(I,2)+U(I,3))/dly4_com
  !
  RETURN
END SUBROUTINE DY4