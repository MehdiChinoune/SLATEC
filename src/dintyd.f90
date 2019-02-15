!DECK DINTYD
SUBROUTINE DINTYD(T,K,Yh,Nyh,Dky,Iflag)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DINTYD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDEBDF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   DINTYD approximates the solution and derivatives at T by polynomial
  !   interpolation. Must be used in conjunction with the integrator
  !   package DDEBDF.
  ! ----------------------------------------------------------------------
  ! DINTYD computes interpolated values of the K-th derivative of the
  ! dependent variable vector Y, and stores it in DKY.
  ! This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may
  ! also be called by the user for any K up to the current order.
  ! (see detailed instructions in LSODE usage documentation.)
  ! ----------------------------------------------------------------------
  ! The computed values in DKY are gotten by interpolation using the
  ! Nordsieck history array YH.  This array corresponds uniquely to a
  ! vector-valued polynomial of degree NQCUR or less, and DKY is set
  ! to the K-th derivative of this polynomial at T.
  ! The formula for DKY is..
  !              Q
  !  DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
  !             J=K
  ! where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
  ! The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
  ! communicated by common.  The above sum is done in reverse order.
  ! IFLAG is returned negative if either K or T is out of bounds.
  ! ----------------------------------------------------------------------
  !
  !***SEE ALSO  DDEBDF
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DDEBD1
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DINTYD
  !
  INTEGER i, ic, IER, Iflag, IOWnd, IOWns, j, jb, jb2, jj, jj1, &
    jp1, JSTart, K, KFLag, L, MAXord, METh, MITer, N, NFE, &
    NJE, NQ, NQU, NST, Nyh
  REAL(8) :: c, Dky, EL0, H, HMIn, HMXi, HU, r, ROWnd, &
    ROWns, s, T, TN, tp, UROund, Yh
  DIMENSION Yh(Nyh,*), Dky(*)
  COMMON /DDEBD1/ ROWnd, ROWns(210), EL0, H, HMIn, HMXi, HU, TN, &
    UROund, IOWnd(14), IOWns(6), IER, JSTart, KFLag, L, &
    METh, MITer, MAXord, N, NQ, NST, NFE, NJE, NQU
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 130
  !***FIRST EXECUTABLE STATEMENT  DINTYD
  Iflag = 0
  IF ( K<0.OR.K>NQ ) THEN
    !
    Iflag = -1
  ELSE
    tp = TN - HU*(1.0D0+100.0D0*UROund)
    IF ( (T-tp)*(T-TN)<=0.0D0 ) THEN
      !
      s = (T-TN)/H
      ic = 1
      IF ( K/=0 ) THEN
        jj1 = L - K
        DO jj = jj1, NQ
          ic = ic*jj
        ENDDO
      ENDIF
      c = ic
      DO i = 1, N
        Dky(i) = c*Yh(i,L)
      ENDDO
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
            ENDDO
          ENDIF
          c = ic
          DO i = 1, N
            Dky(i) = c*Yh(i,jp1) + s*Dky(i)
          ENDDO
        ENDDO
        !     .........EXIT
        IF ( K==0 ) GOTO 99999
      ENDIF
      r = H**(-K)
      DO i = 1, N
        Dky(i) = r*Dky(i)
      ENDDO
    ELSE
      !     .........EXIT
      Iflag = -2
    ENDIF
  ENDIF
  !     ----------------------- END OF SUBROUTINE DINTYD
  !     -----------------------
  99999 CONTINUE
  END SUBROUTINE DINTYD
