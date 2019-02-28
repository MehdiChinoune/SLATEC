!DECK SCOV
SUBROUTINE SCOV(FCN,Iopt,M,N,X,Fvec,R,Ldr,Info,Wa1,Wa2,Wa3,Wa4)
  IMPLICIT NONE
  REAL ENORM
  !***BEGIN PROLOGUE  SCOV
  !***PURPOSE  Calculate the covariance matrix for a nonlinear data
  !            fitting problem.  It is intended to be used after a
  !            successful return from either SNLS1 or SNLS1E.
  !***LIBRARY   SLATEC
  !***CATEGORY  K1B1
  !***TYPE      SINGLE PRECISION (SCOV-S, DCOV-D)
  !***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING,
  !             NONLINEAR LEAST SQUARES
  !***AUTHOR  Hiebert, K. L., (SNLA)
  !***DESCRIPTION
  !
  !  1. Purpose.
  !
  !     SCOV calculates the covariance matrix for a nonlinear data
  !     fitting problem.  It is intended to be used after a
  !     successful return from either SNLS1 or SNLS1E. SCOV
  !     and SNLS1 (and SNLS1E) have compatible parameters.  The
  !     required external subroutine, FCN, is the same
  !     for all three codes, SCOV, SNLS1, and SNLS1E.
  !
  !  2. Subroutine and Type Statements.
  !
  !     SUBROUTINE SCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO,
  !                     WA1,WA2,WA3,WA4)
  !     INTEGER IOPT,M,N,LDR,INFO
  !     REAL X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M)
  !     EXTERNAL FCN
  !
  !  3. Parameters.
  !
  !       FCN is the name of the user-supplied subroutine which calculates
  !         the functions.  If the user wants to supply the Jacobian
  !         (IOPT=2 or 3), then FCN must be written to calculate the
  !         Jacobian, as well as the functions.  See the explanation
  !         of the IOPT argument below.  FCN must be declared in an
  !         EXTERNAL statement in the calling program and should be
  !         written as follows.
  !
  !         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
  !         INTEGER IFLAG,LDFJAC,M,N
  !         REAL X(N),FVEC(M)
  !         ----------
  !         FJAC and LDFJAC may be ignored    , if IOPT=1.
  !         REAL FJAC(LDFJAC,N)               , if IOPT=2.
  !         REAL FJAC(N)                      , if IOPT=3.
  !         ----------
  !           IFLAG will never be zero when FCN is called by SCOV.
  !         RETURN
  !         ----------
  !           If IFLAG=1, calculate the functions at X and return
  !           this vector in FVEC.
  !         RETURN
  !         ----------
  !           If IFLAG=2, calculate the full Jacobian at X and return
  !           this matrix in FJAC.  Note that IFLAG will never be 2 unless
  !           IOPT=2.  FVEC contains the function values at X and must
  !           not be altered.  FJAC(I,J) must be set to the derivative
  !           of FVEC(I) with respect to X(J).
  !         RETURN
  !         ----------
  !           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
  !           and return this vector in FJAC.  Note that IFLAG will
  !           never be 3 unless IOPT=3.  FJAC(J) must be set to
  !           the derivative of FVEC(LDFJAC) with respect to X(J).
  !         RETURN
  !         ----------
  !         END
  !
  !
  !         The value of IFLAG should not be changed by FCN unless the
  !         user wants to terminate execution of SCOV.  In this case, set
  !         IFLAG to a negative integer.
  !
  !
  !    IOPT is an input variable which specifies how the Jacobian will
  !         be calculated.  If IOPT=2 or 3, then the user must supply the
  !         Jacobian, as well as the function values, through the
  !         subroutine FCN.  If IOPT=2, the user supplies the full
  !         Jacobian with one call to FCN.  If IOPT=3, the user supplies
  !         one row of the Jacobian with each call.  (In this manner,
  !         storage can be saved because the full Jacobian is not stored.)
  !         If IOPT=1, the code will approximate the Jacobian by forward
  !         differencing.
  !
  !       M is a positive integer input variable set to the number of
  !         functions.
  !
  !       N is a positive integer input variable set to the number of
  !         variables.  N must not exceed M.
  !
  !       X is an array of length N.  On input X must contain the value
  !         at which the covariance matrix is to be evaluated.  This is
  !         usually the value for X returned from a successful run of
  !         SNLS1 (or SNLS1E).  The value of X will not be changed.
  !
  !    FVEC is an output array of length M which contains the functions
  !         evaluated at X.
  !
  !       R is an output array.  For IOPT=1 and 2, R is an M by N array.
  !         For IOPT=3, R is an N by N array.  On output, if INFO=1,
  !         the upper N by N submatrix of R contains the covariance
  !         matrix evaluated at X.
  !
  !     LDR is a positive integer input variable which specifies
  !         the leading dimension of the array R.  For IOPT=1 and 2,
  !         LDR must not be less than M.  For IOPT=3, LDR must not
  !         be less than N.
  !
  !    INFO is an integer output variable.  If the user has terminated
  !         execution, INFO is set to the (negative) value of IFLAG.  See
  !         description of FCN. Otherwise, INFO is set as follows.
  !
  !         INFO = 0 Improper input parameters (M.LE.0 or N.LE.0).
  !
  !         INFO = 1 Successful return.  The covariance matrix has been
  !                  calculated and stored in the upper N by N
  !                  submatrix of R.
  !
  !         INFO = 2 The Jacobian matrix is singular for the input value
  !                  of X.  The covariance matrix cannot be calculated.
  !                  The upper N by N submatrix of R contains the QR
  !                  factorization of the Jacobian (probably not of
  !                  interest to the user).
  !
  !     WA1 is a work array of length N.
  !     WA2 is a work array of length N.
  !     WA3 is a work array of length N.
  !     WA4 is a work array of length M.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ENORM, FDJAC3, QRFAC, RWUPDT, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   810522  DATE WRITTEN
  !   890505  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Fixed an error message.  (RWC)
  !***END PROLOGUE  SCOV
  !
  !     REVISED 820707-1100
  !     REVISED YYMMDD HHMM
  !
  INTEGER i, idum, iflag, Info, Iopt, j, k, kp1, Ldr, M, N, nm1, &
    nrow
  REAL X(*), R(Ldr,*), Fvec(*), Wa1(*), Wa2(*), Wa3(*), Wa4(*)
  EXTERNAL FCN
  REAL one, sigma, temp, zero
  LOGICAL sing
  SAVE zero, one
  DATA zero/0.E0/, one/1.E0/
  !***FIRST EXECUTABLE STATEMENT  SCOV
  sing = .FALSE.
  iflag = 0
  IF ( M>0.AND.N>0 ) THEN
    !
    !     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N)
    iflag = 1
    CALL FCN(iflag,M,N,X,Fvec,R,Ldr)
    IF ( iflag>=0 ) THEN
      temp = ENORM(M,Fvec)
      sigma = one
      IF ( M/=N ) sigma = temp*temp/(M-N)
      !
      !     CALCULATE THE JACOBIAN
      IF ( Iopt==3 ) THEN
        !
        !     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE
        !     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R.
        !     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.)
        DO j = 1, N
          Wa2(j) = zero
          DO i = 1, N
            R(i,j) = zero
          ENDDO
        ENDDO
        iflag = 3
        DO i = 1, M
          nrow = i
          CALL FCN(iflag,M,N,X,Fvec,Wa1,nrow)
          IF ( iflag<0 ) GOTO 100
          temp = Fvec(i)
          CALL RWUPDT(N,R,Ldr,Wa1,Wa2,temp,Wa3,Wa4)
        ENDDO
      ELSE
        !
        !     STORE THE FULL JACOBIAN USING M*N STORAGE
        IF ( Iopt==1 ) THEN
          !
          !     CODE APPROXIMATES THE JACOBIAN
          CALL FDJAC3(FCN,M,N,X,Fvec,R,Ldr,iflag,zero,Wa4)
        ELSE
          !
          !     USER SUPPLIES THE JACOBIAN
          iflag = 2
          CALL FCN(iflag,M,N,X,Fvec,R,Ldr)
        ENDIF
        IF ( iflag<0 ) GOTO 100
        !
        !     COMPUTE THE QR DECOMPOSITION
        CALL QRFAC(M,N,R,Ldr,.FALSE.,[idum],1,Wa1,Wa1,Wa1)
        DO i = 1, N
          R(i,i) = Wa1(i)
        ENDDO
      ENDIF
      !
      !     CHECK IF R IS SINGULAR.
      DO i = 1, N
        IF ( R(i,i)==zero ) sing = .TRUE.
      ENDDO
      IF ( .NOT.(sing) ) THEN
        !
        !     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE
        !     IN THE UPPER TRIANGLE OF R.
        IF ( N/=1 ) THEN
          nm1 = N - 1
          DO k = 1, nm1
            !
            !     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF THE
            !     IDENTITY MATRIX.
            DO i = 1, N
              Wa1(i) = zero
            ENDDO
            Wa1(k) = one
            !
            R(k,k) = Wa1(k)/R(k,k)
            kp1 = k + 1
            DO i = kp1, N
              !
              !     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*).
              DO j = i, N
                Wa1(j) = Wa1(j) - R(k,i-1)*R(i-1,j)
              ENDDO
              R(k,i) = Wa1(i)/R(i,i)
            ENDDO
          ENDDO
        ENDIF
        R(N,N) = one/R(N,N)
        !
        !     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER
        !     TRIANGLE OF R.
        DO i = 1, N
          DO j = i, N
            temp = zero
            DO k = j, N
              temp = temp + R(i,k)*R(j,k)
            ENDDO
            R(i,j) = temp*sigma
          ENDDO
        ENDDO
        Info = 1
      ENDIF
    ENDIF
  ENDIF
  !
  100 CONTINUE
  IF ( M<=0.OR.N<=0 ) Info = 0
  IF ( iflag<0 ) Info = iflag
  IF ( sing ) Info = 2
  IF ( Info<0 ) CALL XERMSG('SLATEC','SCOV',&
    'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.'&
    ,1,1)
  IF ( Info==0 ) CALL XERMSG('SLATEC','SCOV','INVALID INPUT PARAMETER.',2,1)
  IF ( Info==2 ) CALL XERMSG('SLATEC','SCOV',&
    'SINGULAR JACOBIAN MATRIX, COVARIANCE MATRIX CANNOT BE '&
    //'CALCULATED.',1,1)
END SUBROUTINE SCOV
