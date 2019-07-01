!** DPPVAL
REAL(DP) PURE  FUNCTION DPPVAL(Ldc,C,Xi,Lxi,K,Ideriv,X)
  !> Calculate the value of the IDERIV-th derivative of the
  !  B-spline from the PP-representation.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  E3, K6
  !***
  ! **Type:**      DOUBLE PRECISION (PPVAL-S, DPPVAL-D)
  !***
  ! **Keywords:**  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Written by Carl de Boor and modified by D. E. Amos
  !
  !     Abstract    **** a double precision routine ****
  !         DPPVAL is the PPVALU function of the reference.
  !
  !         DPPVAL calculates (at X) the value of the IDERIV-th
  !         derivative of the B-spline from the PP-representation
  !         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in
  !         the interval XI(J) <= X < XI(J+1) is evaluated, J=1,LXI.
  !         Right limiting values at X=XI(J) are obtained.  DPPVAL will
  !         extrapolate beyond XI(1) and XI(LXI+1).
  !
  !         To obtain left limiting values (left derivatives) at XI(J)
  !         replace LXI by J-1 and set X=XI(J),J=2,LXI+1.
  !
  !     Description of Arguments
  !
  !         Input      C,XI,X are double precision
  !          LDC     - leading dimension of C matrix, LDC >= K
  !          C       - matrix of dimension at least (K,LXI) containing
  !                    right derivatives at break points XI(*).
  !          XI      - break point vector of length LXI+1
  !          LXI     - number of polynomial pieces
  !          K       - order of B-spline, K >= 1
  !          IDERIV  - order of the derivative, 0 <= IDERIV <= K-1
  !                    IDERIV=0 gives the B-spline value
  !          X       - argument, XI(1) <= X <= XI(LXI+1)
  !          INPPV   - an initialization parameter which must be set
  !                    to 1 the first time DPPVAL is called.
  !
  !         Output     DPPVAL is double precision
  !          INPPV   - INPPV contains information for efficient process-
  !                    ing after the initial call and INPPV must not
  !                    be changed by the user.  Distinct splines require
  !                    distinct INPPV parameters.
  !          DPPVAL  - value of the IDERIV-th derivative at X
  !
  !     Error Conditions
  !         Improper input is a fatal error
  !
  !***
  ! **References:**  Carl de Boor, Package for calculating with B-splines,
  !                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
  !                 pp. 441-472.
  !***
  ! **Routines called:**  DINTRV, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER, INTENT(IN) :: Ideriv, K, Ldc, Lxi
  REAL(DP), INTENT(IN) :: C(Ldc,Lxi), X, Xi(Lxi+1)
  INTEGER :: i, j, ndummy, kk, inppv
  REAL(DP) :: dx
  !* FIRST EXECUTABLE STATEMENT  DPPVAL
  DPPVAL = 0._DP
  inppv = 1
  IF( K<1 ) THEN
    ERROR STOP 'DPPVAL : K DOES NOT SATISFY K>=1'
  ELSEIF( Ldc<K ) THEN
    ERROR STOP 'DPPVAL : LDC DOES NOT SATISFY LDC>=K'
  ELSEIF( Lxi<1 ) THEN
    ERROR STOP 'DPPVAL : LXI DOES NOT SATISFY LXI>=1'
  ELSEIF( Ideriv<0 .OR. Ideriv>=K ) THEN
    ERROR STOP 'DPPVAL : IDERIV DOES NOT SATISFY 0<=IDERIV<K'
  ELSE
    i = K - Ideriv
    kk = i
    CALL DINTRV(Xi,Lxi,X,inppv,i,ndummy)
    dx = X - Xi(i)
    j = K
    DO
      DPPVAL = (DPPVAL/kk)*dx + C(j,i)
      j = j - 1
      kk = kk - 1
      IF( kk<=0 ) EXIT
    END DO
  END IF

  RETURN
END FUNCTION DPPVAL