!*==TEVLC.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK TEVLC
SUBROUTINE TEVLC(N,D,E2,Ierr)
  IMPLICIT NONE
  !*--TEVLC5
  !*** Start of declarations inserted by SPAG
  REAL CNV , dhold
  INTEGER IK , K , NCMplx , nhalf , NM , NPP , ntop
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  TEVLC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (TEVLC-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This subroutine finds the eigenvalues of a symmetric tridiagonal
  !     matrix by the rational QL method.
  !
  !     On Input-
  !
  !        N is the order of the matrix,
  !
  !        D contains the diagonal elements of the input matrix,
  !
  !        E2 contains the subdiagonal elements of the input matrix
  !           in its last N-1 positions.  E2(1) is arbitrary.
  !
  !      On Output-
  !
  !        D contains the eigenvalues in ascending order.  If an
  !          error exit is made, the eigenvalues are correct and
  !          ordered for indices 1,2,...IERR-1, but may not be
  !          the smallest eigenvalues,
  !
  !        E2 has been destroyed,
  !
  !        IERR is set to
  !          ZERO       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !
  !***SEE ALSO  CBLKTR
  !***REFERENCES  C. H. Reinsch, Eigenvalues of a real, symmetric, tri-
  !                 diagonal matrix, Algorithm 464, Communications of the
  !                 ACM 16, 11 (November 1973), pp. 689.
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CCBLK
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920528  DESCRIPTION revised and REFERENCES section added.  (WRB)
  !***END PROLOGUE  TEVLC
  !
  INTEGER i , j , l , m , N , ii , l1 , mml , Ierr
  REAL D(*) , E2(*)
  REAL b , c , f , g , h , p , r , s , MAChep
  !
  COMMON /CCBLK / NPP , K , MAChep , CNV , NM , NCMplx , IK
  !***FIRST EXECUTABLE STATEMENT  TEVLC
  Ierr = 0
  IF ( N/=1 ) THEN
    !
    DO i = 2 , N
      E2(i-1) = E2(i)*E2(i)
    ENDDO
    !
    f = 0.0
    b = 0.0
    E2(N) = 0.0
    !
    DO l = 1 , N
      j = 0
      h = MAChep*(ABS(D(l))+SQRT(E2(l)))
      IF ( b<=h ) THEN
        b = h
        c = b*b
      ENDIF
      !
      !     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
      !
      DO m = l , N
        IF ( E2(m)<=c ) EXIT
        !
        !     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
        !                THROUGH THE BOTTOM OF THE LOOP **********
        !
      ENDDO
      !
      IF ( m/=l ) THEN
        DO WHILE ( j/=30 )
          j = j + 1
          !
          !     ********** FORM SHIFT **********
          !
          l1 = l + 1
          s = SQRT(E2(l))
          g = D(l)
          p = (D(l1)-g)/(2.0*s)
          r = SQRT(p*p+1.0)
          D(l) = s/(p+SIGN(r,p))
          h = g - D(l)
          !
          DO i = l1 , N
            D(i) = D(i) - h
          ENDDO
          !
          f = f + h
          !
          !     ********** RATIONAL QL TRANSFORMATION **********
          !
          g = D(m)
          IF ( g==0.0 ) g = b
          h = g
          s = 0.0
          mml = m - l
          !
          !     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
          !
          DO ii = 1 , mml
            i = m - ii
            p = g*h
            r = p + E2(i)
            E2(i+1) = s*r
            s = E2(i)/r
            D(i+1) = h + s*(h+D(i))
            g = D(i) - E2(i)/g
            IF ( g==0.0 ) g = b
            h = g*p/r
          ENDDO
          !
          E2(l) = s*g
          D(l) = h
          !
          !     ********** GUARD AGAINST UNDERFLOWED H **********
          !
          IF ( h==0.0 ) GOTO 20
          IF ( ABS(E2(l))<=ABS(c/h) ) GOTO 20
          E2(l) = h*E2(l)
          IF ( E2(l)==0.0 ) GOTO 20
        ENDDO
        GOTO 50
      ENDIF
      20       p = D(l) + f
      !
      !     ********** ORDER EIGENVALUES **********
      !
      IF ( l/=1 ) THEN
        !
        !     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
        !
        DO ii = 2 , l
          i = l + 2 - ii
          IF ( p>=D(i-1) ) GOTO 40
          D(i) = D(i-1)
        ENDDO
      ENDIF
      !
      i = 1
      40       D(i) = p
    ENDDO
    !
    IF ( ABS(D(N))<ABS(D(1)) ) THEN
      nhalf = N/2
      DO i = 1 , nhalf
        ntop = N - i
        dhold = D(i)
        D(i) = D(ntop+1)
        D(ntop+1) = dhold
      ENDDO
    ENDIF
    GOTO 99999
    !
    !     ********** SET ERROR -- NO CONVERGENCE TO AN
    !                EIGENVALUE AFTER 30 ITERATIONS **********
    !
    50     Ierr = l
  ENDIF
  !
  !     ********** LAST CARD OF TQLRAT **********
  !
  99999 END SUBROUTINE TEVLC
