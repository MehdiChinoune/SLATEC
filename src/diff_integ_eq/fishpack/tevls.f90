!** TEVLS
SUBROUTINE TEVLS(N,D,E2,Ierr)
  !>
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (TEVLS-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
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
  !***
  ! **See also:**  BLKTRI
  !***
  ! **References:**  C. H. Reinsch, Eigenvalues of a REAL(SP), symmetric, tri-
  !                 diagonal matrix, Algorithm 464, Communications of the
  !                 ACM 16, 11 (November 1973), pp. 689.
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CBLKT

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920528  DESCRIPTION revised and REFERENCES section added.  (WRB)
  USE CBLKT, ONLY : eps_com
  INTEGER :: N, Ierr
  REAL(SP) :: D(N), E2(N)
  INTEGER :: nhalf, ntop, i, j, l, m, ii, l1, mml
  REAL(SP) :: dhold, b, c, f, g, h, p, r, s
  !* FIRST EXECUTABLE STATEMENT  TEVLS
  Ierr = 0
  IF ( N/=1 ) THEN
    !
    DO i = 2, N
      E2(i-1) = E2(i)*E2(i)
    END DO
    !
    f = 0.0
    b = 0.0
    E2(N) = 0.0
    !
    DO l = 1, N
      j = 0
      h = eps_com*(ABS(D(l))+SQRT(E2(l)))
      IF ( b<=h ) THEN
        b = h
        c = b*b
      END IF
      !
      !     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
      !
      DO m = l, N
        IF ( E2(m)<=c ) EXIT
        !
        !     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
        !                THROUGH THE BOTTOM OF THE LOOP **********
        !
      END DO
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
          DO i = l1, N
            D(i) = D(i) - h
          END DO
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
          DO ii = 1, mml
            i = m - ii
            p = g*h
            r = p + E2(i)
            E2(i+1) = s*r
            s = E2(i)/r
            D(i+1) = h + s*(h+D(i))
            g = D(i) - E2(i)/g
            IF ( g==0.0 ) g = b
            h = g*p/r
          END DO
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
        END DO
        GOTO 50
      END IF
      20  p = D(l) + f
      !
      !     ********** ORDER EIGENVALUES **********
      !
      IF ( l/=1 ) THEN
        !
        !     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
        !
        DO ii = 2, l
          i = l + 2 - ii
          IF ( p>=D(i-1) ) GOTO 40
          D(i) = D(i-1)
        END DO
      END IF
      !
      i = 1
      40  D(i) = p
    END DO
    !
    IF ( ABS(D(N))<ABS(D(1)) ) THEN
      nhalf = N/2
      DO i = 1, nhalf
        ntop = N - i
        dhold = D(i)
        D(i) = D(ntop+1)
        D(ntop+1) = dhold
      END DO
    END IF
    RETURN
    !
    !     ********** SET ERROR -- NO CONVERGENCE TO AN
    !                EIGENVALUE AFTER 30 ITERATIONS **********
    !
    50  Ierr = l
  END IF
  !
  !     ********** LAST CARD OF TQLRAT **********
  !
  RETURN
END SUBROUTINE TEVLS
