!** DBSPDR
PURE SUBROUTINE DBSPDR(T,A,N,K,Nderiv,Ad)
  !> Use the B-representation to construct a divided difference
  !  table preparatory to a (right) derivative calculation.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E3, K6
  !***
  ! **Type:**      DOUBLE PRECISION (BSPDR-S, DBSPDR-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, DIFFERENTIATION OF SPLINES,
  !             INTERPOLATION
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract     **** a double precision routine ****
  !         DBSPDR is the BSPLDR routine of the reference.
  !
  !         DBSPDR uses the B-representation (T,A,N,K) to construct a
  !         divided difference table ADIF preparatory to a (right)
  !         derivative calculation in DBSPEV.  The lower triangular matrix
  !         ADIF is stored in vector AD by columns.  The arrays are
  !         related by
  !
  !           ADIF(I,J) = AD(I-J+1 + (2*N-J+2)*(J-1)/2)
  !
  !         I = J,N  ,   J=1,NDERIV.
  !
  !     Description of Arguments
  !
  !         Input      T,A are double precision
  !          T       - knot vector of length N+K
  !          A       - B-spline coefficient vector of length N
  !          N       - number of B-spline coefficients
  !                    N = sum of knot multiplicities-K
  !          K       - order of the spline, K >= 1
  !          NDERIV  - number of derivatives, 1 <= NDERIV <= K.
  !                    NDERIV=1 gives the zero-th derivative =
  !                    function value
  !
  !         Output     AD is double precision
  !          AD      - table of differences in a vector of length
  !                    (2*N-NDERIV+1)*NDERIV/2 for input to DBSPEV
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***
  ! **References:**  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER, INTENT(IN) :: K, N, Nderiv
  REAL(DP), INTENT(IN) :: A(N), T(N+K)
  REAL(DP), INTENT(OUT) :: Ad((2*N-Nderiv+1)*Nderiv/2)
  INTEGER :: i, id, ii, ipkmid, jj, jm, kmid
  REAL(DP) :: diff, fkmid
  !     DIMENSION T(N+K), AD((2*N-NDERIV+1)*NDERIV/2)
  !* FIRST EXECUTABLE STATEMENT  DBSPDR
  IF( K<1 ) THEN
    ERROR STOP 'DBSPDR : K DOES NOT SATISFY K>=1'
  ELSEIF( N<K ) THEN
    ERROR STOP 'DBSPDR : N DOES NOT SATISFY N>=K'
  ELSEIF( Nderiv<1 .OR. Nderiv>K ) THEN
    ERROR STOP 'DBSPDR : NDERIV DOES NOT SATISFY 1<=NDERIV<=K'
  END IF
  DO i = 1, N
    Ad(i) = A(i)
  END DO
  IF( Nderiv==1 ) RETURN
  kmid = K
  jj = N
  jm = 0
  DO id = 2, Nderiv
    kmid = kmid - 1
    fkmid = kmid
    ii = 1
    DO i = id, N
      ipkmid = i + kmid
      diff = T(ipkmid) - T(i)
      IF( diff/=0._DP ) Ad(ii+jj) = (Ad(ii+jm+1)-Ad(ii+jm))/diff*fkmid
      ii = ii + 1
    END DO
    jm = jj
    jj = jj + N - id + 1
  END DO

  RETURN
END SUBROUTINE DBSPDR