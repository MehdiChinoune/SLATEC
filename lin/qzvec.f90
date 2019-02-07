!*==QZVEC.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK QZVEC
SUBROUTINE QZVEC(Nm,N,A,B,Alfr,Alfi,Beta,Z)
  IMPLICIT NONE
  !*--QZVEC5
  !***BEGIN PROLOGUE  QZVEC
  !***PURPOSE  The optional fourth step of the QZ algorithm for
  !            generalized eigenproblems.  Accepts a matrix in
  !            quasi-triangular form and another in upper triangular
  !            and computes the eigenvectors of the triangular problem
  !            and transforms them back to the original coordinates
  !            Usually preceded by QZHES, QZIT, and QZVAL.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C3
  !***TYPE      SINGLE PRECISION (QZVEC-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is the optional fourth step of the QZ algorithm
  !     for solving generalized matrix eigenvalue problems,
  !     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
  !
  !     This subroutine accepts a pair of REAL matrices, one of them in
  !     quasi-triangular form (in which each 2-by-2 block corresponds to
  !     a pair of complex eigenvalues) and the other in upper triangular
  !     form.  It computes the eigenvectors of the triangular problem and
  !     transforms the results back to the original coordinate system.
  !     It is usually preceded by  QZHES,  QZIT, and  QZVAL.
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
  !        ALFR, ALFI, and BETA are one-dimensional REAL arrays with
  !          components whose ratios ((ALFR+I*ALFI)/BETA) are the
  !          generalized eigenvalues.  They are usually obtained from
  !          QZVAL.  They are dimensioned ALFR(N), ALFI(N), and BETA(N).
  !
  !        Z contains the transformation matrix produced in the reductions
  !          by  QZHES,  QZIT, and  QZVAL,  if performed.  If the
  !          eigenvectors of the triangular problem are desired, Z must
  !          contain the identity matrix.  Z is a two-dimensional REAL
  !          array, dimensioned Z(NM,N).
  !
  !     On Output
  !
  !        A is unaltered.  Its subdiagonal elements provide information
  !           about the storage of the complex eigenvectors.
  !
  !        B has been destroyed.
  !
  !        ALFR, ALFI, and BETA are unaltered.
  !
  !        Z contains the real and imaginary parts of the eigenvectors.
  !          If ALFI(J) .EQ. 0.0, the J-th eigenvalue is real and
  !            the J-th column of Z contains its eigenvector.
  !          If ALFI(J) .NE. 0.0, the J-th eigenvalue is complex.
  !            If ALFI(J) .GT. 0.0, the eigenvalue is the first of
  !              a complex pair and the J-th and (J+1)-th columns
  !              of Z contain its eigenvector.
  !            If ALFI(J) .LT. 0.0, the eigenvalue is the second of
  !              a complex pair and the (J-1)-th and J-th columns
  !              of Z contain the conjugate of its eigenvector.
  !          Each eigenvector is normalized so that the modulus
  !          of its largest component is 1.0 .
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
  !***END PROLOGUE  QZVEC
  !
  INTEGER i, j, k, m, N, en, ii, jj, na, Nm, nn, isw, enm2
  REAL A(Nm,*), B(Nm,*), Alfr(*), Alfi(*), Beta(*), Z(Nm,*)
  REAL d, q, r, s, t, w, x, y, di, dr, ra, rr, sa, ti, tr, &
    t1, t2
  REAL w1, x1, zz, z1, alfm, almi, almr, betm, epsb
  !
  !***FIRST EXECUTABLE STATEMENT  QZVEC
  epsb = B(N,1)
  isw = 1
  !     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
  DO nn = 1, N
    en = N + 1 - nn
    na = en - 1
    IF ( isw/=2 ) THEN
      IF ( Alfi(en)/=0.0E0 ) THEN
        !     .......... COMPLEX VECTOR ..........
        m = na
        almr = Alfr(m)
        almi = Alfi(m)
        betm = Beta(m)
        !     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
        !                EIGENVECTOR MATRIX IS TRIANGULAR ..........
        y = betm*A(en,na)
        B(na,na) = -almi*B(en,en)/y
        B(na,en) = (almr*B(en,en)-betm*A(en,en))/y
        B(en,na) = 0.0E0
        B(en,en) = 1.0E0
        enm2 = na - 1
        IF ( enm2/=0 ) THEN
          !     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
          DO ii = 1, enm2
            i = na - ii
            w = betm*A(i,i) - almr*B(i,i)
            w1 = -almi*B(i,i)
            ra = 0.0E0
            sa = 0.0E0
            !
            DO j = m, en
              x = betm*A(i,j) - almr*B(i,j)
              x1 = -almi*B(i,j)
              ra = ra + x*B(j,na) - x1*B(j,en)
              sa = sa + x*B(j,en) + x1*B(j,na)
            ENDDO
            !
            IF ( i/=1.AND.isw/=2 ) THEN
              IF ( betm*A(i,i-1)/=0.0E0 ) THEN
                zz = w
                z1 = w1
                r = ra
                s = sa
                isw = 2
                CYCLE
              ENDIF
            ENDIF
            m = i
            IF ( isw==2 ) GOTO 6
            !     .......... COMPLEX 1-BY-1 BLOCK ..........
            tr = -ra
            ti = -sa
            2              dr = w
            di = w1
            !     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
            4              IF ( ABS(di)<=ABS(dr) ) THEN
            rr = di/dr
            d = dr + di*rr
            t1 = (tr+ti*rr)/d
            t2 = (ti-tr*rr)/d
            IF ( isw==1 ) GOTO 10
            IF ( isw==2 ) GOTO 8
        ENDIF
        rr = dr/di
        d = dr*rr + di
        t1 = (tr*rr+ti)/d
        t2 = (ti*rr-tr)/d
        IF ( isw==1 ) GOTO 10
        IF ( isw==2 ) GOTO 8
        !     .......... COMPLEX 2-BY-2 BLOCK ..........
        6              x = betm*A(i,i+1) - almr*B(i,i+1)
        x1 = -almi*B(i,i+1)
        y = betm*A(i+1,i)
        tr = y*ra - w*r + w1*s
        ti = y*sa - w*s - w1*r
        dr = w*zz - w1*z1 - x*y
        di = w*z1 + w1*zz - x1*y
        IF ( dr==0.0E0.AND.di==0.0E0 ) dr = epsb
        GOTO 4
        8              B(i+1,na) = t1
        B(i+1,en) = t2
        isw = 1
        IF ( ABS(y)>ABS(w)+ABS(w1) ) THEN
          t1 = (-r-zz*B(i+1,na)+z1*B(i+1,en))/y
          t2 = (-s-zz*B(i+1,en)-z1*B(i+1,na))/y
        ELSE
          tr = -ra - x*B(i+1,na) + x1*B(i+1,en)
          ti = -sa - x*B(i+1,en) - x1*B(i+1,na)
          GOTO 2
        ENDIF
        10             B(i,na) = t1
        B(i,en) = t2
  ENDDO
ENDIF
ELSE
!     .......... REAL VECTOR ..........
m = en
B(en,en) = 1.0E0
IF ( na/=0 ) THEN
  alfm = Alfr(m)
  betm = Beta(m)
  !     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
  DO ii = 1, na
    i = en - ii
    w = betm*A(i,i) - alfm*B(i,i)
    r = 0.0E0
    !
    DO j = m, en
      r = r + (betm*A(i,j)-alfm*B(i,j))*B(j,en)
    ENDDO
    !
    IF ( i/=1.AND.isw/=2 ) THEN
      IF ( betm*A(i,i-1)/=0.0E0 ) THEN
        zz = w
        s = r
        GOTO 12
      ENDIF
    ENDIF
    m = i
    IF ( isw==2 ) THEN
      !     .......... REAL 2-BY-2 BLOCK ..........
      x = betm*A(i,i+1) - alfm*B(i,i+1)
      y = betm*A(i+1,i)
      q = w*zz - x*y
      t = (x*s-zz*r)/q
      B(i,en) = t
      IF ( ABS(x)<=ABS(zz) ) THEN
        B(i+1,en) = (-s-y*t)/zz
      ELSE
        B(i+1,en) = (-r-w*t)/x
      ENDIF
    ELSE
      !     .......... REAL 1-BY-1 BLOCK ..........
      t = w
      IF ( w==0.0E0 ) t = epsb
      B(i,en) = -r/t
      CYCLE
    ENDIF
    12             isw = 3 - isw
    !     .......... END REAL VECTOR ..........
  ENDDO
ENDIF
CYCLE
ENDIF
ENDIF
!     .......... END COMPLEX VECTOR ..........
isw = 3 - isw
ENDDO
!     .......... END BACK SUBSTITUTION.
!                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
!                FOR J=N STEP -1 UNTIL 1 DO -- ..........
DO jj = 1, N
j = N + 1 - jj
!
DO i = 1, N
zz = 0.0E0
!
DO k = 1, j
zz = zz + Z(i,k)*B(k,j)
ENDDO
!
Z(i,j) = zz
ENDDO
ENDDO
!     .......... NORMALIZE SO THAT MODULUS OF LARGEST
!                COMPONENT OF EACH VECTOR IS 1.
!                (ISW IS 1 INITIALLY FROM BEFORE) ..........
DO j = 1, N
d = 0.0E0
IF ( isw==2 ) THEN
!
DO i = 1, N
r = ABS(Z(i,j-1)) + ABS(Z(i,j))
IF ( r/=0.0E0 ) r = r*SQRT((Z(i,j-1)/r)**2+(Z(i,j)/r)**2)
IF ( r>d ) d = r
ENDDO
!
DO i = 1, N
Z(i,j-1) = Z(i,j-1)/d
Z(i,j) = Z(i,j)/d
ENDDO
ELSEIF ( Alfi(j)==0.0E0 ) THEN
!
DO i = 1, N
IF ( ABS(Z(i,j))>d ) d = ABS(Z(i,j))
ENDDO
!
DO i = 1, N
Z(i,j) = Z(i,j)/d
ENDDO
!
CYCLE
ENDIF
!
isw = 3 - isw
ENDDO
!
END SUBROUTINE QZVEC
