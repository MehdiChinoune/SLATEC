!DECK QZIT
SUBROUTINE QZIT(Nm,N,A,B,Eps1,Matz,Z,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QZIT
  !***PURPOSE  The second step of the QZ algorithm for generalized
  !            eigenproblems.  Accepts an upper Hessenberg and an upper
  !            triangular matrix and reduces the former to
  !            quasi-triangular form while preserving the form of the
  !            latter.  Usually preceded by QZHES and followed by QZVAL
  !            and QZVEC.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B3
  !***TYPE      SINGLE PRECISION (QZIT-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is the second step of the QZ algorithm
  !     for solving generalized matrix eigenvalue problems,
  !     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART,
  !     as modified in technical note NASA TN D-7305(1973) by WARD.
  !
  !     This subroutine accepts a pair of REAL matrices, one of them
  !     in upper Hessenberg form and the other in upper triangular form.
  !     It reduces the Hessenberg matrix to quasi-triangular form using
  !     orthogonal transformations while maintaining the triangular form
  !     of the other matrix.  It is usually preceded by  QZHES  and
  !     followed by  QZVAL  and, possibly,  QZVEC.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A, B, and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrices A and B.  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        A contains a real upper Hessenberg matrix.  A is a two-
  !          dimensional REAL array, dimensioned A(NM,N).
  !
  !        B contains a real upper triangular matrix.  B is a two-
  !          dimensional REAL array, dimensioned B(NM,N).
  !
  !        EPS1 is a tolerance used to determine negligible elements.
  !          EPS1 = 0.0 (or negative) may be input, in which case an
  !          element will be neglected only if it is less than roundoff
  !          error times the norm of its matrix.  If the input EPS1 is
  !          positive, then an element will be considered negligible
  !          if it is less than EPS1 times the norm of its matrix.  A
  !          positive value of EPS1 may result in faster execution,
  !          but less accurate results.  EPS1 is a REAL variable.
  !
  !        MATZ should be set to .TRUE. if the right hand transformations
  !          are to be accumulated for later use in computing
  !          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
  !          variable.
  !
  !        Z contains, if MATZ has been set to .TRUE., the transformation
  !          matrix produced in the reduction by  QZHES, if performed, or
  !          else the identity matrix.  If MATZ has been set to .FALSE.,
  !          Z is not referenced.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,N).
  !
  !     On Output
  !
  !        A has been reduced to quasi-triangular form.  The elements
  !          below the first subdiagonal are still zero, and no two
  !          consecutive subdiagonal elements are nonzero.
  !
  !        B is still in upper triangular form, although its elements
  !          have been altered.  The location B(N,1) is used to store
  !          EPS1 times the norm of B for later use by  QZVAL  and  QZVEC.
  !
  !        Z contains the product of the right hand transformations
  !          (for both steps) if MATZ has been set to .TRUE.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if neither A(J,J-1) nor A(J-1,J-2) has become
  !                     zero after a total of 30*N iterations.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  QZIT
  !
  INTEGER i, j, k, l, N, en, k1, k2, ld, ll, l1, na, Nm, ish, &
    itn, its, km1, lm1
  INTEGER enm2, Ierr, lor1, enorn
  REAL A(Nm,*), B(Nm,*), Z(Nm,*)
  REAL r, s, t, a1, a2, a3, ep, sh, u1, u2, u3, v1, v2, v3, &
    ani
  REAL a11, a12, a21, a22, a33, a34, a43, a44, bni, b11
  REAL b12, b22, b33, b34, b44, epsa, epsb, Eps1, anorm, bnorm
  LOGICAL Matz, notlas
  !
  !***FIRST EXECUTABLE STATEMENT  QZIT
  Ierr = 0
  !     .......... COMPUTE EPSA,EPSB ..........
  anorm = 0.0E0
  bnorm = 0.0E0
  !
  DO i = 1, N
    ani = 0.0E0
    IF ( i/=1 ) ani = ABS(A(i,i-1))
    bni = 0.0E0
    !
    DO j = i, N
      ani = ani + ABS(A(i,j))
      bni = bni + ABS(B(i,j))
    ENDDO
    !
    IF ( ani>anorm ) anorm = ani
    IF ( bni>bnorm ) bnorm = bni
  ENDDO
  !
  IF ( anorm==0.0E0 ) anorm = 1.0E0
  IF ( bnorm==0.0E0 ) bnorm = 1.0E0
  ep = Eps1
  IF ( ep<=0.0E0 ) THEN
    !     .......... COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO ..........
    ep = 1.0E0
    DO
      ep = ep/2.0E0
      IF ( 1.0E0+ep<=1.0E0 ) EXIT
    ENDDO
  ENDIF
  epsa = ep*anorm
  epsb = ep*bnorm
  !     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
  !                KEEPING B TRIANGULAR ..........
  lor1 = 1
  enorn = N
  en = N
  itn = 30*N
  !     .......... BEGIN QZ STEP ..........
  100 CONTINUE
  IF ( en<=2 ) THEN
    !     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
    IF ( N>1 ) B(N,1) = epsb
    GOTO 99999
  ELSE
    IF ( .NOT.Matz ) enorn = en
    its = 0
    na = en - 1
    enm2 = na - 1
  ENDIF
  200  ish = 2
  !     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
  !                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
  DO ll = 1, en
    lm1 = en - ll
    l = lm1 + 1
    IF ( l==1 ) GOTO 400
    IF ( ABS(A(l,lm1))<=epsa ) EXIT
  ENDDO
  !
  300  A(l,lm1) = 0.0E0
  IF ( l>=na ) THEN
    !     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
    en = lm1
    GOTO 100
  ENDIF
  !     .......... CHECK FOR SMALL TOP OF B ..........
  400  ld = l
  DO
    l1 = l + 1
    b11 = B(l,l)
    IF ( ABS(b11)>epsb ) THEN
      a11 = A(l,l)/b11
      a21 = A(l1,l)/b11
      IF ( ish/=1 ) THEN
        !     .......... ITERATION STRATEGY ..........
        IF ( itn==0 ) THEN
          !     .......... SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT
          !                HAS BECOME NEGLIGIBLE AFTER 30*N ITERATIONS ..........
          Ierr = en
          IF ( N>1 ) B(N,1) = epsb
          GOTO 99999
        ELSEIF ( its==10 ) THEN
          !     .......... AD HOC SHIFT ..........
          a1 = 0.0E0
          a2 = 1.0E0
          a3 = 1.1605E0
          EXIT
        ELSE
          !     .......... DETERMINE TYPE OF SHIFT ..........
          b22 = B(l1,l1)
          IF ( ABS(b22)<epsb ) b22 = epsb
          b33 = B(na,na)
          IF ( ABS(b33)<epsb ) b33 = epsb
          b44 = B(en,en)
          IF ( ABS(b44)<epsb ) b44 = epsb
          a33 = A(na,na)/b33
          a34 = A(na,en)/b44
          a43 = A(en,na)/b33
          a44 = A(en,en)/b44
          b34 = B(na,en)/b44
          t = 0.5E0*(a43*b34-a33-a44)
          r = t*t + a34*a43 - a33*a44
          IF ( r<0.0E0 ) THEN
            !     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
            a12 = A(l,l1)/b22
            a22 = A(l1,l1)/b22
            b12 = B(l,l1)/b22
            a1 = ((a33-a11)*(a44-a11)-a34*a43+a43*b34*a11)/a21 + a12 - &
              a11*b12
            a2 = (a22-a11) - a21*b12 - (a33-a11) - (a44-a11) + a43*b34
            a3 = A(l1+1,l1)/b22
            EXIT
          ELSE
            !     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
            ish = 1
            r = SQRT(r)
            sh = -t + r
            s = -t - r
            IF ( ABS(s-a44)<ABS(sh-a44) ) sh = s
            !     .......... LOOK FOR TWO CONSECUTIVE SMALL
            !                SUB-DIAGONAL ELEMENTS OF A.
            !                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
            DO ll = ld, enm2
              l = enm2 + ld - ll
              IF ( l==ld ) EXIT
              lm1 = l - 1
              l1 = l + 1
              t = A(l,l)
              IF ( ABS(B(l,l))>epsb ) t = t - sh*B(l,l)
              IF ( ABS(A(l,lm1))<=ABS(t/A(l1,l))*epsa ) GOTO 500
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      !
      a1 = a11 - sh
      a2 = a21
      IF ( l/=ld ) A(l,lm1) = -A(l,lm1)
      EXIT
    ELSE
      B(l,l) = 0.0E0
      s = ABS(A(l,l)) + ABS(A(l1,l))
      u1 = A(l,l)/s
      u2 = A(l1,l)/s
      r = SIGN(SQRT(u1*u1+u2*u2),u1)
      v1 = -(u1+r)/r
      v2 = -u2/r
      u2 = v2/v1
      !
      DO j = l, enorn
        t = A(l,j) + u2*A(l1,j)
        A(l,j) = A(l,j) + t*v1
        A(l1,j) = A(l1,j) + t*v2
        t = B(l,j) + u2*B(l1,j)
        B(l,j) = B(l,j) + t*v1
        B(l1,j) = B(l1,j) + t*v2
      ENDDO
      !
      IF ( l/=1 ) A(l,lm1) = -A(l,lm1)
      lm1 = l
      l = l1
      GOTO 300
    ENDIF
    500  ENDDO
    its = its + 1
    itn = itn - 1
    IF ( .NOT.Matz ) lor1 = ld
    !     .......... MAIN LOOP ..........
    DO k = l, na
      notlas = k/=na .AND. ish==2
      k1 = k + 1
      k2 = k + 2
      km1 = MAX(k-1,l)
      ll = MIN(en,k1+ish)
      IF ( notlas ) THEN
        !     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
        IF ( k/=l ) THEN
          a1 = A(k,km1)
          a2 = A(k1,km1)
          a3 = A(k2,km1)
        ENDIF
        s = ABS(a1) + ABS(a2) + ABS(a3)
        IF ( s==0.0E0 ) CYCLE
        u1 = a1/s
        u2 = a2/s
        u3 = a3/s
        r = SIGN(SQRT(u1*u1+u2*u2+u3*u3),u1)
        v1 = -(u1+r)/r
        v2 = -u2/r
        v3 = -u3/r
        u2 = v2/v1
        u3 = v3/v1
        !
        DO j = km1, enorn
          t = A(k,j) + u2*A(k1,j) + u3*A(k2,j)
          A(k,j) = A(k,j) + t*v1
          A(k1,j) = A(k1,j) + t*v2
          A(k2,j) = A(k2,j) + t*v3
          t = B(k,j) + u2*B(k1,j) + u3*B(k2,j)
          B(k,j) = B(k,j) + t*v1
          B(k1,j) = B(k1,j) + t*v2
          B(k2,j) = B(k2,j) + t*v3
        ENDDO
        !
        IF ( k/=l ) THEN
          A(k1,km1) = 0.0E0
          A(k2,km1) = 0.0E0
        ENDIF
        !     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
        s = ABS(B(k2,k2)) + ABS(B(k2,k1)) + ABS(B(k2,k))
        IF ( s/=0.0E0 ) THEN
          u1 = B(k2,k2)/s
          u2 = B(k2,k1)/s
          u3 = B(k2,k)/s
          r = SIGN(SQRT(u1*u1+u2*u2+u3*u3),u1)
          v1 = -(u1+r)/r
          v2 = -u2/r
          v3 = -u3/r
          u2 = v2/v1
          u3 = v3/v1
          !
          DO i = lor1, ll
            t = A(i,k2) + u2*A(i,k1) + u3*A(i,k)
            A(i,k2) = A(i,k2) + t*v1
            A(i,k1) = A(i,k1) + t*v2
            A(i,k) = A(i,k) + t*v3
            t = B(i,k2) + u2*B(i,k1) + u3*B(i,k)
            B(i,k2) = B(i,k2) + t*v1
            B(i,k1) = B(i,k1) + t*v2
            B(i,k) = B(i,k) + t*v3
          ENDDO
          !
          B(k2,k) = 0.0E0
          B(k2,k1) = 0.0E0
          IF ( Matz ) THEN
            !
            DO i = 1, N
              t = Z(i,k2) + u2*Z(i,k1) + u3*Z(i,k)
              Z(i,k2) = Z(i,k2) + t*v1
              Z(i,k1) = Z(i,k1) + t*v2
              Z(i,k) = Z(i,k) + t*v3
            ENDDO
          ENDIF
        ENDIF
      ELSE
        !     .......... ZERO A(K+1,K-1) ..........
        IF ( k/=l ) THEN
          a1 = A(k,km1)
          a2 = A(k1,km1)
        ENDIF
        s = ABS(a1) + ABS(a2)
        IF ( s==0.0E0 ) EXIT
        u1 = a1/s
        u2 = a2/s
        r = SIGN(SQRT(u1*u1+u2*u2),u1)
        v1 = -(u1+r)/r
        v2 = -u2/r
        u2 = v2/v1
        !
        DO j = km1, enorn
          t = A(k,j) + u2*A(k1,j)
          A(k,j) = A(k,j) + t*v1
          A(k1,j) = A(k1,j) + t*v2
          t = B(k,j) + u2*B(k1,j)
          B(k,j) = B(k,j) + t*v1
          B(k1,j) = B(k1,j) + t*v2
        ENDDO
        !
        IF ( k/=l ) A(k1,km1) = 0.0E0
      ENDIF
      !     .......... ZERO B(K+1,K) ..........
      s = ABS(B(k1,k1)) + ABS(B(k1,k))
      IF ( s/=0.0E0 ) THEN
        u1 = B(k1,k1)/s
        u2 = B(k1,k)/s
        r = SIGN(SQRT(u1*u1+u2*u2),u1)
        v1 = -(u1+r)/r
        v2 = -u2/r
        u2 = v2/v1
        !
        DO i = lor1, ll
          t = A(i,k1) + u2*A(i,k)
          A(i,k1) = A(i,k1) + t*v1
          A(i,k) = A(i,k) + t*v2
          t = B(i,k1) + u2*B(i,k)
          B(i,k1) = B(i,k1) + t*v1
          B(i,k) = B(i,k) + t*v2
        ENDDO
        !
        B(k1,k) = 0.0E0
        IF ( Matz ) THEN
          !
          DO i = 1, N
            t = Z(i,k1) + u2*Z(i,k)
            Z(i,k1) = Z(i,k1) + t*v1
            Z(i,k) = Z(i,k) + t*v2
          ENDDO
        ENDIF
      ENDIF
      !
    ENDDO
    !     .......... END QZ STEP ..........
    GOTO 200
      99999 CONTINUE
  END SUBROUTINE QZIT
