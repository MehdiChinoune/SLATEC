!DECK CPZERO
SUBROUTINE CPZERO(In,A,R,T,Iflg,S)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CPZERO
  !***PURPOSE  Find the zeros of a polynomial with complex coefficients.
  !***LIBRARY   SLATEC
  !***CATEGORY  F1A1B
  !***TYPE      COMPLEX (RPZERO-S, CPZERO-C)
  !***KEYWORDS  POLYNOMIAL ROOTS, POLYNOMIAL ZEROS, REAL ROOTS
  !***AUTHOR  Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  !      Find the zeros of the complex polynomial
  !         P(Z)= A(1)*Z**N + A(2)*Z**(N-1) +...+ A(N+1)
  !
  !    Input...
  !       IN = degree of P(Z)
  !       A = complex vector containing coefficients of P(Z),
  !            A(I) = coefficient of Z**(N+1-i)
  !       R = N word complex vector containing initial estimates for zeros
  !            if these are known.
  !       T = 4(N+1) word array used for temporary storage
  !       IFLG = flag to indicate if initial estimates of
  !              zeros are input.
  !            If IFLG .EQ. 0, no estimates are input.
  !            If IFLG .NE. 0, the vector R contains estimates of
  !               the zeros
  !       ** WARNING ****** If estimates are input, they must
  !                         be separated, that is, distinct or
  !                         not repeated.
  !       S = an N word array
  !
  !    Output...
  !       R(I) = Ith zero,
  !       S(I) = bound for R(I) .
  !       IFLG = error diagnostic
  !    Error Diagnostics...
  !       If IFLG .EQ. 0 on return, all is well
  !       If IFLG .EQ. 1 on return, A(1)=0.0 or N=0 on input
  !       If IFLG .EQ. 2 on return, the program failed to converge
  !                after 25*N iterations.  Best current estimates of the
  !                zeros are in R(I).  Error bounds are not calculated.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CPEVL
  !***REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CPZERO
  INTEGER i, Iflg, imax, In, j, n, n1, nit, nmax, nr
  REAL u, v, x
  REAL S(*)
  COMPLEX R(*), T(*), A(*), pn(1), temp(1)
  !***FIRST EXECUTABLE STATEMENT  CPZERO
  IF ( In<=0.OR.ABS(A(1))==0.0 ) THEN
    Iflg = 1
    RETURN
  ELSE
    !
    !       CHECK FOR EASILY OBTAINED ZEROS
    !
    n = In
    n1 = n + 1
    IF ( Iflg==0 ) THEN
      DO
        n1 = n + 1
        IF ( n<=1 ) THEN
          R(1) = -A(2)/A(1)
          S(1) = 0.0
          RETURN
        ELSEIF ( ABS(A(n1))/=0.0 ) THEN
          !
          !          IF INITIAL ESTIMATES FOR ZEROS NOT GIVEN, FIND SOME
          !
          temp = -A(2)/(A(1)*n)
          CALL CPEVL(n,n,A,temp(1),T,T,.FALSE.)
          imax = n + 2
          T(n1) = ABS(T(n1))
          DO i = 2, n1
            T(n+i) = -ABS(T(n+2-i))
            IF ( REAL(T(n+i))<REAL(T(imax)) ) imax = n + i
          ENDDO
          x = (-REAL(T(imax))/REAL(T(n1)))**(1./(imax-n1))
          DO
            x = 2.*x
            CALL CPEVL(n,0,T(n1),CMPLX(x,0.0),pn,pn,.FALSE.)
            IF ( REAL(pn(1))>=0. ) THEN
              u = .5*x
              v = x
              DO
                x = .5*(u+v)
                CALL CPEVL(n,0,T(n1),CMPLX(x,0.0),pn,pn,.FALSE.)
                IF ( REAL(pn(1))>0. ) v = x
                IF ( REAL(pn(1))<=0. ) u = x
                IF ( (v-u)<=.001*(1.+v) ) THEN
                  DO i = 1, n
                    u = (3.14159265/n)*(2*i-1.5)
                    R(i) = MAX(x,.001*ABS(temp(1)))*CMPLX(COS(u),SIN(u)) + temp(1)
                  ENDDO
                  GOTO 50
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ELSE
          R(n) = 0.0
          S(n) = 0.0
          n = n - 1
        ENDIF
      ENDDO
    ENDIF
    !
    !          MAIN ITERATION LOOP STARTS HERE
    !
    50     nr = 0
    nmax = 25*n
    DO nit = 1, nmax
      DO i = 1, n
        IF ( nit==1.OR.ABS(T(i))/=0. ) THEN
          CALL CPEVL(n,0,A,R(i),pn,temp,.TRUE.)
          IF ( ABS(REAL(pn(1)))+ABS(AIMAG(pn(1)))>REAL(temp(1))+AIMAG(temp(1)) ) THEN
            temp = A(1)
            DO j = 1, n
              IF ( j/=i ) temp = temp*(R(i)-R(j))
            ENDDO
            T(i) = pn(1)/temp(1)
          ELSE
            T(i) = 0.0
            nr = nr + 1
          ENDIF
        ENDIF
      ENDDO
      DO i = 1, n
        R(i) = R(i) - T(i)
      ENDDO
      IF ( nr==n ) GOTO 100
    ENDDO
    GOTO 200
  ENDIF
  !
  !          CALCULATE ERROR BOUNDS FOR ZEROS
  !
  100 CONTINUE
  DO nr = 1, n
    CALL CPEVL(n,n,A,R(nr),T,T(n+2),.TRUE.)
    x = ABS(CMPLX(ABS(REAL(T(1))),ABS(AIMAG(T(1))))+T(n+2))
    S(nr) = 0.0
    DO i = 1, n
      x = x*REAL(n1-i)/i
      temp = CMPLX(MAX(ABS(REAL(T(i+1)))-REAL(T(n1+i)),0.0),&
        MAX(ABS(AIMAG(T(i+1)))-AIMAG(T(n1+i)),0.0))
      S(nr) = MAX(S(nr),(ABS(temp(1))/x)**(1./i))
    ENDDO
    S(nr) = 1./S(nr)
  ENDDO
  RETURN
  !        ERROR EXITS
  200 CONTINUE
  IFlg = 2
  RETURN
END SUBROUTINE CPZERO
