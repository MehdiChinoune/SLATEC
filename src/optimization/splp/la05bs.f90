!** LA05BS
PURE SUBROUTINE LA05BS(A,Ind,Ia,N,Ip,Iw,W,G,B,Trans)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (LA05BS-S, LA05BD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
  !     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
  !     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
  !     THE FINAL LETTER =S= IN THE NAMES USED HERE.
  !     REVISED SEP. 13, 1979.
  !
  !     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
  !     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
  !     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
  !     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
  !     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
  !
  ! IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U.
  ! IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U.
  ! IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  XERMSG, XSETUN
  !***
  ! COMMON BLOCKS    LA05DS

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900402  Added TYPE section.  (WRB)
  !   920410  Corrected second dimension on IW declaration.  (WRB)
  USE LA05DS, ONLY : lp_com, lenl_com

  INTEGER, INTENT(IN) :: Ia, N
  LOGICAL, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Ind(Ia,2), Iw(N,8)
  INTEGER, INTENT(INOUT) :: Ip(N,2)
  REAL(SP), INTENT(IN) :: G, A(Ia)
  REAL(SP), INTENT(INOUT) :: B(:)
  REAL(SP), INTENT(OUT) :: W(:)
  INTEGER :: i, ii, j, k, k2, kk, kl, kll, kp, kpc, l1, n1, nz
  REAL(SP) :: am
  !* FIRST EXECUTABLE STATEMENT  LA05BS
  IF( G<0. ) THEN
    !
    IF( lp_com>0 ) ERROR STOP 'LA05BS : EARLIER ENTRY GAVE ERROR RETURN.'
  ELSE
    kll = Ia - lenl_com + 1
    IF( Trans ) THEN
      !
      !     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U
      DO i = 1, N
        W(i) = B(i)
        B(i) = 0._SP
      END DO
      DO ii = 1, N
        i = Iw(ii,4)
        am = W(i)
        IF( am/=0. ) THEN
          j = Iw(ii,3)
          kp = Ip(j,1)
          am = am/A(kp)
          B(j) = am
          kl = Iw(j,1) + kp - 1
          IF( kp/=kl ) THEN
            k2 = kp + 1
            DO k = k2, kl
              i = Ind(k,2)
              W(i) = W(i) - am*A(k)
            END DO
          END IF
        END IF
      END DO
      !
      !     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L
      IF( kll>Ia ) RETURN
      DO k = kll, Ia
        j = Ind(k,2)
        IF( B(j)/=0. ) THEN
          i = Ind(k,1)
          B(i) = B(i) + A(k)*B(j)
        END IF
      END DO
    ELSE
      !
      !     MULTIPLY VECTOR BY INVERSE OF L
      IF( lenl_com>0 ) THEN
        l1 = Ia + 1
        DO kk = 1, lenl_com
          k = l1 - kk
          i = Ind(k,1)
          IF( B(i)/=0. ) THEN
            j = Ind(k,2)
            B(j) = B(j) + A(k)*B(i)
          END IF
        END DO
      END IF
      DO i = 1, N
        W(i) = B(i)
        B(i) = 0._SP
      END DO
      !
      !     MULTIPLY VECTOR BY INVERSE OF U
      n1 = N + 1
      DO ii = 1, N
        i = n1 - ii
        i = Iw(i,3)
        am = W(i)
        kp = Ip(i,1)
        IF( kp<=0 ) THEN
          kp = -kp
          Ip(i,1) = kp
          nz = Iw(i,1)
          kl = kp - 1 + nz
          k2 = kp + 1
          DO k = k2, kl
            j = Ind(k,2)
            am = am - A(k)*B(j)
          END DO
        END IF
        IF( am/=0. ) THEN
          j = Ind(kp,2)
          B(j) = am/A(kp)
          kpc = Ip(j,2)
          kl = Iw(j,2) + kpc - 1
          IF( kl/=kpc ) THEN
            k2 = kpc + 1
            DO k = k2, kl
              i = Ind(k,1)
              Ip(i,1) = -ABS(Ip(i,1))
            END DO
          END IF
        END IF
      END DO
    END IF
  END IF

END SUBROUTINE LA05BS