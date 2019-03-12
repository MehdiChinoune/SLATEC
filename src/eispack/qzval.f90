!DECK QZVAL
SUBROUTINE QZVAL(Nm,N,A,B,Alfr,Alfi,Beta,Matz,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QZVAL
  !***PURPOSE  The third step of the QZ algorithm for generalized
  !            eigenproblems.  Accepts a pair of real matrices, one in
  !            quasi-triangular form and the other in upper triangular
  !            form and computes the eigenvalues of the associated
  !            eigenproblem.  Usually preceded by QZHES, QZIT, and
  !            followed by QZVEC.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C2C
  !***TYPE      SINGLE PRECISION (QZVAL-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is the third step of the QZ algorithm
  !     for solving generalized matrix eigenvalue problems,
  !     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
  !
  !     This subroutine accepts a pair of REAL matrices, one of them
  !     in quasi-triangular form and the other in upper triangular form.
  !     It reduces the quasi-triangular matrix further, so that any
  !     remaining 2-by-2 blocks correspond to pairs of complex
  !     eigenvalues, and returns quantities whose ratios give the
  !     generalized eigenvalues.  It is usually preceded by  QZHES
  !     and  QZIT  and may be followed by  QZVEC.
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
  !        A contains a real upper quasi-triangular matrix.  A is a two-
  !          dimensional REAL array, dimensioned A(NM,N).
  !
  !        B contains a real upper triangular matrix.  In addition,
  !          location B(N,1) contains the tolerance quantity (EPSB)
  !          computed and saved in  QZIT.  B is a two-dimensional REAL
  !          array, dimensioned B(NM,N).
  !
  !        MATZ should be set to .TRUE. if the right hand transformations
  !          are to be accumulated for later use in computing
  !          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
  !          variable.
  !
  !        Z contains, if MATZ has been set to .TRUE., the transformation
  !          matrix produced in the reductions by  QZHES  and  QZIT,  if
  !          performed, or else the identity matrix.  If MATZ has been set
  !          to .FALSE., Z is not referenced.  Z is a two-dimensional REAL
  !          array, dimensioned Z(NM,N).
  !
  !     On Output
  !
  !        A has been reduced further to a quasi-triangular matrix in
  !          which all nonzero subdiagonal elements correspond to pairs
  !          of complex eigenvalues.
  !
  !        B is still in upper triangular form, although its elements
  !          have been altered.  B(N,1) is unaltered.
  !
  !        ALFR and ALFI contain the real and imaginary parts of the
  !          diagonal elements of the triangular matrix that would be
  !          obtained if A were reduced completely to triangular form
  !          by unitary transformations.  Non-zero values of ALFI occur
  !          in pairs, the first member positive and the second negative.
  !          ALFR and ALFI are one-dimensional REAL arrays, dimensioned
  !          ALFR(N) and ALFI(N).
  !
  !        BETA contains the diagonal elements of the corresponding B,
  !          normalized to be real and non-negative.  The generalized
  !          eigenvalues are then the ratios ((ALFR+I*ALFI)/BETA).
  !          BETA is a one-dimensional REAL array, dimensioned BETA(N).
  !
  !        Z contains the product of the right hand transformations
  !          (for all three steps) if MATZ has been set to .TRUE.
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  QZVAL
  !
  INTEGER i, j, N, en, na, Nm, nn, isw
  REAL A(Nm,*), B(Nm,*), Alfr(*), Alfi(*), Beta(*), Z(Nm,*)
  REAL c, d, e, r, s, t, an, a1, a2, bn, cq, cz, di, dr, ei, &
    ti, tr
  REAL u1, u2, v1, v2, a1i, a11, a12, a2i, a21, a22, b11, b12, &
    b22
  REAL sqi, sqr, ssi, ssr, szi, szr, a11i, a11r, a12i, a12r
  REAL a22i, a22r, epsb
  LOGICAL Matz
  !
  !***FIRST EXECUTABLE STATEMENT  QZVAL
  epsb = B(N,1)
  isw = 1
  !     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
  !                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
  DO nn = 1, N
    en = N + 1 - nn
    na = en - 1
    IF ( isw==2 ) GOTO 250
    IF ( en/=1 ) THEN
      IF ( A(en,na)/=0.0E0 ) THEN
        !     .......... 2-BY-2 BLOCK ..........
        IF ( ABS(B(na,na))<=epsb ) GOTO 100
        IF ( ABS(B(en,en))>epsb ) THEN
          an = ABS(A(na,na)) + ABS(A(na,en)) + ABS(A(en,na)) + ABS(A(en,en))
          bn = ABS(B(na,na)) + ABS(B(na,en)) + ABS(B(en,en))
          a11 = A(na,na)/an
          a12 = A(na,en)/an
          a21 = A(en,na)/an
          a22 = A(en,en)/an
          b11 = B(na,na)/bn
          b12 = B(na,en)/bn
          b22 = B(en,en)/bn
          e = a11/b11
          ei = a22/b22
          s = a21/(b11*b22)
          t = (a22-e*b22)/b22
          IF ( ABS(e)>ABS(ei) ) THEN
            e = ei
            t = (a11-e*b11)/b11
          ENDIF
          c = 0.5E0*(t-s*b12)
          d = c*c + s*(a12-e*b12)
          IF ( d<0.0E0 ) THEN
            !     .......... TWO COMPLEX ROOTS ..........
            e = e + c
            ei = SQRT(-d)
            a11r = a11 - e*b11
            a11i = ei*b11
            a12r = a12 - e*b12
            a12i = ei*b12
            a22r = a22 - e*b22
            a22i = ei*b22
            IF ( ABS(a11r)+ABS(a11i)+ABS(a12r)+ABS(a12i)<ABS(a21)+ABS(a22r)&
                +ABS(a22i) ) THEN
              a1 = a22r
              a1i = a22i
              a2 = -a21
              a2i = 0.0E0
            ELSE
              a1 = a12r
              a1i = a12i
              a2 = -a11r
              a2i = -a11i
            ENDIF
            !     .......... CHOOSE COMPLEX Z ..........
            cz = SQRT(a1*a1+a1i*a1i)
            IF ( cz==0.0E0 ) THEN
              szr = 1.0E0
              szi = 0.0E0
            ELSE
              szr = (a1*a2+a1i*a2i)/cz
              szi = (a1*a2i-a1i*a2)/cz
              r = SQRT(cz*cz+szr*szr+szi*szi)
              cz = cz/r
              szr = szr/r
              szi = szi/r
            ENDIF
            IF ( an<(ABS(e)+ei)*bn ) THEN
              a1 = cz*a11 + szr*a12
              a1i = szi*a12
              a2 = cz*a21 + szr*a22
              a2i = szi*a22
            ELSE
              a1 = cz*b11 + szr*b12
              a1i = szi*b12
              a2 = szr*b22
              a2i = szi*b22
            ENDIF
            !     .......... CHOOSE COMPLEX Q ..........
            cq = SQRT(a1*a1+a1i*a1i)
            IF ( cq==0.0E0 ) THEN
              sqr = 1.0E0
              sqi = 0.0E0
            ELSE
              sqr = (a1*a2+a1i*a2i)/cq
              sqi = (a1*a2i-a1i*a2)/cq
              r = SQRT(cq*cq+sqr*sqr+sqi*sqi)
              cq = cq/r
              sqr = sqr/r
              sqi = sqi/r
            ENDIF
            !     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
            !                IF TRANSFORMATIONS WERE APPLIED ..........
            ssr = sqr*szr + sqi*szi
            ssi = sqr*szi - sqi*szr
            i = 1
            tr = cq*cz*a11 + cq*szr*a12 + sqr*cz*a21 + ssr*a22
            ti = cq*szi*a12 - sqi*cz*a21 + ssi*a22
            dr = cq*cz*b11 + cq*szr*b12 + ssr*b22
            di = cq*szi*b12 + ssi*b22
            DO
              t = ti*dr - tr*di
              j = na
              IF ( t<0.0E0 ) j = en
              r = SQRT(dr*dr+di*di)
              Beta(j) = bn*r
              Alfr(j) = an*(tr*dr+ti*di)/r
              Alfi(j) = an*t/r
              IF ( i/=1 ) GOTO 250
              i = 2
              tr = ssr*a11 - sqr*cz*a12 - cq*szr*a21 + cq*cz*a22
              ti = -ssi*a11 - sqi*cz*a12 + cq*szi*a21
              dr = ssr*b11 - sqr*cz*b12 + cq*cz*b22
              di = -ssi*b11 - sqi*cz*b12
            ENDDO
          ELSE
            !     .......... TWO REAL ROOTS.
            !                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
            e = e + (c+SIGN(SQRT(d),c))
            a11 = a11 - e*b11
            a12 = a12 - e*b12
            a22 = a22 - e*b22
            IF ( ABS(a11)+ABS(a12)<ABS(a21)+ABS(a22) ) THEN
              a1 = a22
              a2 = a21
            ELSE
              a1 = a12
              a2 = a11
            ENDIF
            GOTO 50
          ENDIF
        ELSE
          a1 = A(en,en)
          a2 = A(en,na)
          bn = 0.0E0
          GOTO 50
        ENDIF
      ENDIF
    ENDIF
    !     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
    Alfr(en) = A(en,en)
    IF ( B(en,en)<0.0E0 ) Alfr(en) = -Alfr(en)
    Beta(en) = ABS(B(en,en))
    Alfi(en) = 0.0E0
    CYCLE
    !     .......... CHOOSE AND APPLY REAL Z ..........
    50     s = ABS(a1) + ABS(a2)
    u1 = a1/s
    u2 = a2/s
    r = SIGN(SQRT(u1*u1+u2*u2),u1)
    v1 = -(u1+r)/r
    v2 = -u2/r
    u2 = v2/v1
    !
    DO i = 1, en
      t = A(i,en) + u2*A(i,na)
      A(i,en) = A(i,en) + t*v1
      A(i,na) = A(i,na) + t*v2
      t = B(i,en) + u2*B(i,na)
      B(i,en) = B(i,en) + t*v1
      B(i,na) = B(i,na) + t*v2
    ENDDO
    !
    IF ( Matz ) THEN
      !
      DO i = 1, N
        t = Z(i,en) + u2*Z(i,na)
        Z(i,en) = Z(i,en) + t*v1
        Z(i,na) = Z(i,na) + t*v2
      ENDDO
    ENDIF
    !
    IF ( bn==0.0E0 ) GOTO 200
    IF ( an>=ABS(e)*bn ) THEN
      a1 = B(na,na)
      a2 = B(en,na)
      GOTO 150
    ENDIF
    100    a1 = A(na,na)
    a2 = A(en,na)
    !     .......... CHOOSE AND APPLY REAL Q ..........
    150    s = ABS(a1) + ABS(a2)
    IF ( s/=0.0E0 ) THEN
      u1 = a1/s
      u2 = a2/s
      r = SIGN(SQRT(u1*u1+u2*u2),u1)
      v1 = -(u1+r)/r
      v2 = -u2/r
      u2 = v2/v1
      !
      DO j = na, N
        t = A(na,j) + u2*A(en,j)
        A(na,j) = A(na,j) + t*v1
        A(en,j) = A(en,j) + t*v2
        t = B(na,j) + u2*B(en,j)
        B(na,j) = B(na,j) + t*v1
        B(en,j) = B(en,j) + t*v2
      ENDDO
    ENDIF
    !
    200    A(en,na) = 0.0E0
    B(en,na) = 0.0E0
    Alfr(na) = A(na,na)
    Alfr(en) = A(en,en)
    IF ( B(na,na)<0.0E0 ) Alfr(na) = -Alfr(na)
    IF ( B(en,en)<0.0E0 ) Alfr(en) = -Alfr(en)
    Beta(na) = ABS(B(na,na))
    Beta(en) = ABS(B(en,en))
    Alfi(en) = 0.0E0
    Alfi(na) = 0.0E0
    250    isw = 3 - isw
  ENDDO
  !
END SUBROUTINE QZVAL
