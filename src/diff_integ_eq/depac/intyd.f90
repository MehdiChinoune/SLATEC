!** INTYD
SUBROUTINE INTYD(T,K,Yh,Nyh,Dky,Iflag)
  USE DEBDF1, ONLY : H, HU, TN, UROund, L, N, NQ
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (INTYD-S, DINTYD-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   INTYD approximates the solution and derivatives at T by polynomial
  !   interpolation. Must be used in conjunction with the integrator
  !   package DEBDF.
  ! ----------------------------------------------------------------------
  ! INTYD computes interpolated values of the K-th derivative of the
  ! dependent variable vector Y, and stores it in DKY.
  ! This routine is called by DEBDF with K = 0,1 and T = TOUT, but may
  ! also be called by the user for any K up to the current order.
  ! (see detailed instructions in LSODE usage documentation.)
  ! ----------------------------------------------------------------------
  ! The computed values in DKY are gotten by interpolation using the
  ! Nordsieck history array YH.  This array corresponds uniquely to a
  ! vector-valued polynomial of degree NQCUR or less, and DKY is set
  ! to the K-th derivative of this polynomial at T.
  ! The formula for DKY is..
  !              Q
  !  DKY(I)  =  sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
  !             J=K
  ! where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
  ! The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
  ! communicated by common.  The above sum is done in reverse order.
  ! IFLAG is returned negative if either K or T is out of bounds.
  ! ----------------------------------------------------------------------
  !
  !***
  ! **See also:**  DEBDF
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DEBDF1

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)

  !
  !LLL. OPTIMIZE
  INTEGER K, Nyh, Iflag, i, ic, j, jb, jb2, jj, jj1, jp1
  REAL T, Yh(Nyh,*), Dky(*), c, r, s, tp
  !
  !* FIRST EXECUTABLE STATEMENT  INTYD
  Iflag = 0
  IF ( K<0.OR.K>NQ ) THEN
    !
    Iflag = -1
    RETURN
  ELSE
    tp = TN - HU*(1.0E0+100.0E0*UROund)
    IF ( (T-tp)*(T-TN)>0.0E0 ) THEN
      Iflag = -2
      RETURN
    ELSE
      !
      s = (T-TN)/H
      ic = 1
      IF ( K/=0 ) THEN
        jj1 = L - K
        DO jj = jj1, NQ
          ic = ic*jj
        END DO
      END IF
      c = ic
      DO i = 1, N
        Dky(i) = c*Yh(i,L)
      END DO
      IF ( K/=NQ ) THEN
        jb2 = NQ - K
        DO jb = 1, jb2
          j = NQ - jb
          jp1 = j + 1
          ic = 1
          IF ( K/=0 ) THEN
            jj1 = jp1 - K
            DO jj = jj1, j
              ic = ic*jj
            END DO
          END IF
          c = ic
          DO i = 1, N
            Dky(i) = c*Yh(i,jp1) + s*Dky(i)
          END DO
        END DO
        IF ( K==0 ) RETURN
      END IF
    END IF
  END IF
  r = H**(-K)
  DO i = 1, N
    Dky(i) = r*Dky(i)
  END DO
  RETURN
  !----------------------- END OF SUBROUTINE INTYD -----------------------
  RETURN
END SUBROUTINE INTYD
