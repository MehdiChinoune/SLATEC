!DECK D1UPDT
SUBROUTINE D1UPDT(M,N,S,Ls,U,V,W,Sing)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  D1UPDT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DNSQ and DNSQE
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (R1UPDT-S, D1UPDT-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Given an M by N lower trapezoidal matrix S, an M-vector U,
  !     and an N-vector V, the problem is to determine an
  !     orthogonal matrix Q such that
  !
  !                   t
  !           (S + U*V )*Q
  !
  !     is again lower trapezoidal.
  !
  !     This subroutine determines Q as the product of 2*(N - 1)
  !     transformations
  !
  !           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
  !
  !     where GV(I), GW(I) are Givens rotations in the (I,N) plane
  !     which eliminate elements in the I-th and N-th planes,
  !     respectively. Q itself is not accumulated, rather the
  !     information to recover the GV, GW rotations is returned.
  !
  !     The SUBROUTINE statement is
  !
  !       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
  !
  !     where
  !
  !       M is a positive integer input variable set to the number
  !         of rows of S.
  !
  !       N is a positive integer input variable set to the number
  !         of columns of S. N must not exceed M.
  !
  !       S is an array of length LS. On input S must contain the lower
  !         trapezoidal matrix S stored by columns. On output S contains
  !         the lower trapezoidal matrix produced as described above.
  !
  !       LS is a positive integer input variable not less than
  !         (N*(2*M-N+1))/2.
  !
  !       U is an input array of length M which must contain the
  !         vector U.
  !
  !       V is an array of length N. On input V must contain the vector
  !         V. On output V(I) contains the information necessary to
  !         recover the Givens rotation GV(I) described above.
  !
  !       W is an output array of length M. W(I) contains information
  !         necessary to recover the Givens rotation GW(I) described
  !         above.
  !
  !       SING is a LOGICAL output variable. SING is set TRUE if any
  !         of the diagonal elements of the output S are zero. Otherwise
  !         SING is set FALSE.
  !
  !***SEE ALSO  DNSQ, DNSQE
  !***ROUTINES CALLED  D1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  D1UPDT
  REAL(8) :: D1MACH
  INTEGER i, j, jj, l, Ls, M, N, nm1, nmj
  REAL(8) :: cos, cotan, giant, one, p25, p5, S(*), sin, tan, &
    tau, temp, U(*), V(*), W(*), zero
  LOGICAL Sing
  SAVE one, p5, p25, zero
  DATA one, p5, p25, zero/1.0D0, 5.0D-1, 2.5D-1, 0.0D0/
  !
  !     GIANT IS THE LARGEST MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  D1UPDT
  giant = D1MACH(2)
  !
  !     INITIALIZE THE DIAGONAL ELEMENT POINTER.
  !
  jj = (N*(2*M-N+1))/2 - (M-N)
  !
  !     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
  !
  l = jj
  DO i = N, M
    W(i) = S(l)
    l = l + 1
  ENDDO
  !
  !     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
  !     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
  !
  nm1 = N - 1
  IF ( nm1>=1 ) THEN
    DO nmj = 1, nm1
      j = N - nmj
      jj = jj - (M-j+1)
      W(j) = zero
      IF ( V(j)/=zero ) THEN
        !
        !        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
        !        J-TH ELEMENT OF V.
        !
        IF ( ABS(V(N))>=ABS(V(j)) ) THEN
          tan = V(j)/V(N)
          cos = p5/SQRT(p25+p25*tan**2)
          sin = cos*tan
          tau = sin
        ELSE
          cotan = V(N)/V(j)
          sin = p5/SQRT(p25+p25*cotan**2)
          cos = sin*cotan
          tau = one
          IF ( ABS(cos)*giant>one ) tau = one/cos
        ENDIF
        !
        !        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
        !        NECESSARY TO RECOVER THE GIVENS ROTATION.
        !
        V(N) = sin*V(j) + cos*V(N)
        V(j) = tau
        !
        !        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
        !
        l = jj
        DO i = j, M
          temp = cos*S(l) - sin*W(i)
          W(i) = sin*S(l) + cos*W(i)
          S(l) = temp
          l = l + 1
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !
  !     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
  !
  DO i = 1, M
    W(i) = W(i) + V(N)*U(i)
  ENDDO
  !
  !     ELIMINATE THE SPIKE.
  !
  Sing = .FALSE.
  IF ( nm1>=1 ) THEN
    DO j = 1, nm1
      IF ( W(j)/=zero ) THEN
        !
        !        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
        !        J-TH ELEMENT OF THE SPIKE.
        !
        IF ( ABS(S(jj))>=ABS(W(j)) ) THEN
          tan = W(j)/S(jj)
          cos = p5/SQRT(p25+p25*tan**2)
          sin = cos*tan
          tau = sin
        ELSE
          cotan = S(jj)/W(j)
          sin = p5/SQRT(p25+p25*cotan**2)
          cos = sin*cotan
          tau = one
          IF ( ABS(cos)*giant>one ) tau = one/cos
        ENDIF
        !
        !        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
        !
        l = jj
        DO i = j, M
          temp = cos*S(l) + sin*W(i)
          W(i) = -sin*S(l) + cos*W(i)
          S(l) = temp
          l = l + 1
        ENDDO
        !
        !        STORE THE INFORMATION NECESSARY TO RECOVER THE
        !        GIVENS ROTATION.
        !
        W(j) = tau
      ENDIF
      !
      !        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
      !
      IF ( S(jj)==zero ) Sing = .TRUE.
      jj = jj + (M-j+1)
    ENDDO
  ENDIF
  !
  !     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
  !
  l = jj
  DO i = N, M
    S(l) = W(i)
    l = l + 1
  ENDDO
  IF ( S(jj)==zero ) Sing = .TRUE.
  !
  !     LAST CARD OF SUBROUTINE D1UPDT.
  !
END SUBROUTINE D1UPDT
