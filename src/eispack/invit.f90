!** INVIT
SUBROUTINE INVIT(Nm,N,A,Wr,Wi,Select,Mm,M,Z,Ierr,Rm1,Rv1,Rv2)
  IMPLICIT NONE
  !>
  !***
  !  Compute the eigenvectors of a real upper Hessenberg
  !            matrix associated with specified eigenvalues by inverse
  !            iteration.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      SINGLE PRECISION (INVIT-S, CINVIT-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure INVIT
  !     by Peters and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
  !
  !     This subroutine finds those eigenvectors of a REAL UPPER
  !     Hessenberg matrix corresponding to specified eigenvalues,
  !     using inverse iteration.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        A contains the upper Hessenberg matrix.  A is a two-dimensional
  !          REAL array, dimensioned A(NM,N).
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues of the Hessenberg matrix.  The eigenvalues
  !          must be stored in a manner identical to that output by
  !          subroutine  HQR,  which recognizes possible splitting of the
  !          matrix.  WR and WI are one-dimensional REAL arrays,
  !          dimensioned WR(N) and WI(N).
  !
  !        SELECT specifies the eigenvectors to be found. The
  !          eigenvector corresponding to the J-th eigenvalue is
  !          specified by setting SELECT(J) to .TRUE.  SELECT is a
  !          one-dimensional LOGICAL array, dimensioned SELECT(N).
  !
  !        MM should be set to an upper bound for the number of
  !          columns required to store the eigenvectors to be found.
  !          NOTE that two columns are required to store the
  !          eigenvector corresponding to a complex eigenvalue.  One
  !          column is required to store the eigenvector corresponding
  !          to a real eigenvalue.  MM is an INTEGER variable.
  !
  !     On OUTPUT
  !
  !        A and WI are unaltered.
  !
  !        WR may have been altered since close eigenvalues are perturbed
  !          slightly in searching for independent eigenvectors.
  !
  !        SELECT may have been altered.  If the elements corresponding
  !          to a pair of conjugate complex eigenvalues were each
  !          initially set to .TRUE., the program resets the second of
  !          the two elements to .FALSE.
  !
  !        M is the number of columns actually used to store the
  !          eigenvectors.  M is an INTEGER variable.
  !
  !        Z contains the real and imaginary parts of the eigenvectors.
  !          The eigenvectors are packed into the columns of Z starting
  !          at the first column.  If the next selected eigenvalue is
  !          real, the next column of Z contains its eigenvector.  If the
  !          eigenvalue is complex, the next two columns of Z contain the
  !          real and imaginary parts of its eigenvector, with the real
  !          part first.  The eigenvectors are normalized so that the
  !          component of largest magnitude is 1. Any vector which fails
  !          the acceptance test is set to zero.  Z is a two-dimensional
  !          REAL array, dimensioned Z(NM,MM).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          -(2*N+1)   if more than MM columns of Z are necessary
  !                     to store the eigenvectors corresponding to
  !                     the specified eigenvalues (in this case, M is
  !                     equal to the number of columns of Z containing
  !                     eigenvectors already computed),
  !          -K         if the iteration corresponding to the K-th
  !                     value fails (if this occurs more than once, K
  !                     is the index of the last occurrence); the
  !                     corresponding columns of Z are set to zero
  !                     vectors,
  !          -(N+K)     if both error situations occur.
  !
  !        RM1 is a two-dimensional REAL array used for temporary storage.
  !          This array holds the triangularized form of the upper
  !          Hessenberg matrix used in the inverse iteration process.
  !          RM1 is dimensioned RM1(N,N).
  !
  !        RV1 and RV2 are one-dimensional REAL arrays used for temporary
  !          storage.  They hold the approximate eigenvectors during the
  !          inverse iteration process.  RV1 and RV2 are dimensioned
  !          RV1(N) and RV2(N).
  !
  !     The ALGOL procedure GUESSVEC appears in INVIT in-line.
  !
  !     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
  !     Calls CDIV for complex division.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***
  ! **References:**  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***
  ! **Routines called:**  CDIV, PYTHAG

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, M, N, s, ii, ip, Mm, mp, Nm, ns, n1, &
    uk, ip1, its, km1, Ierr
  REAL A(Nm,*), Wr(*), Wi(*), Z(Nm,*)
  REAL Rm1(N,*), Rv1(*), Rv2(*)
  REAL t, w, x, y, eps3
  REAL norm, normv, growto, ilambd, rlambd, ukroot
  REAL PYTHAG
  LOGICAL Select(N)
  !
  !* FIRST EXECUTABLE STATEMENT  INVIT
  Ierr = 0
  uk = 0
  s = 1
  !     .......... IP = 0, REAL EIGENVALUE
  !                     1, FIRST OF CONJUGATE COMPLEX PAIR
  !                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
  ip = 0
  n1 = N - 1
  !
  DO k = 1, N
    IF ( Wi(k)/=0.0E0.AND.ip>=0 ) THEN
      ip = 1
      IF ( Select(k).AND.Select(k+1) ) Select(k+1) = .FALSE.
    END IF
    IF ( .NOT.Select(k) ) GOTO 400
    IF ( Wi(k)/=0.0E0 ) s = s + 1
    IF ( s>Mm ) GOTO 500
    IF ( uk<k ) THEN
      !     .......... CHECK FOR POSSIBLE SPLITTING ..........
      DO uk = k, N
        IF ( uk==N ) EXIT
        IF ( A(uk+1,uk)==0.0E0 ) EXIT
      END DO
      !     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
      !                (HESSENBERG) MATRIX ..........
      norm = 0.0E0
      mp = 1
      !
      DO i = 1, uk
        x = 0.0E0
        !
        DO j = mp, uk
          x = x + ABS(A(i,j))
        END DO
        !
        IF ( x>norm ) norm = x
        mp = i
      END DO
      !     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
      !                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
      IF ( norm==0.0E0 ) norm = 1.0E0
      eps3 = norm
      DO
        eps3 = 0.5E0*eps3
        IF ( norm+eps3<=norm ) THEN
          eps3 = 2.0E0*eps3
          !     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
          ukroot = SQRT(REAL(uk))
          growto = 0.1E0/ukroot
          EXIT
        END IF
      END DO
    END IF
    rlambd = Wr(k)
    ilambd = Wi(k)
    IF ( k==1 ) GOTO 150
    km1 = k - 1
    GOTO 100
    !     .......... PERTURB EIGENVALUE IF IT IS CLOSE
    !                TO ANY PREVIOUS EIGENVALUE ..........
    50  rlambd = rlambd + eps3
    !     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
    100 CONTINUE
    DO ii = 1, km1
      i = k - ii
      IF ( Select(i).AND.ABS(Wr(i)-rlambd)<eps3.AND.ABS(Wi(i)-ilambd)<eps3 )&
        GOTO 50
    END DO
    !
    Wr(k) = rlambd
    !     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
    ip1 = k + ip
    Wr(ip1) = rlambd
    !     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
    !                AND INITIAL REAL VECTOR ..........
    150  mp = 1
    !
    DO i = 1, uk
      !
      DO j = mp, uk
        Rm1(j,i) = A(i,j)
      END DO
      !
      Rm1(i,i) = Rm1(i,i) - rlambd
      mp = i
      Rv1(i) = eps3
    END DO
    !
    its = 0
    IF ( ilambd/=0.0E0 ) THEN
      !     .......... COMPLEX EIGENVALUE.
      !                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
      !                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
      !                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
      ns = N - s
      Z(1,s-1) = -ilambd
      Z(1,s) = 0.0E0
      IF ( N/=2 ) THEN
        Rm1(1,3) = -ilambd
        Z(1,s-1) = 0.0E0
        IF ( N/=3 ) THEN
          !
          DO i = 4, N
            Rm1(1,i) = 0.0E0
          END DO
        END IF
      END IF
      !
      DO i = 2, uk
        mp = i - 1
        w = Rm1(mp,i)
        IF ( i<N ) t = Rm1(mp,i+1)
        IF ( i==N ) t = Z(mp,s-1)
        x = Rm1(mp,mp)*Rm1(mp,mp) + t*t
        IF ( w*w<=x ) THEN
          IF ( x==0.0E0 ) THEN
            Rm1(mp,mp) = eps3
            IF ( i<N ) Rm1(mp,i+1) = 0.0E0
            IF ( i==N ) Z(mp,s-1) = 0.0E0
            t = 0.0E0
            x = eps3*eps3
          END IF
          w = w/x
          x = Rm1(mp,mp)*w
          y = -t*w
          !
          DO j = i, uk
            IF ( j<n1 ) THEN
              t = Rm1(mp,j+2)
              Rm1(i,j+2) = -x*t - y*Rm1(j,mp)
            ELSE
              l = j - ns
              t = Z(mp,l)
              Z(i,l) = -x*t - y*Rm1(j,mp)
            END IF
            Rm1(j,i) = Rm1(j,i) - x*Rm1(j,mp) + y*t
          END DO
          !
          IF ( i<n1 ) THEN
            Rm1(i,i+2) = Rm1(i,i+2) - ilambd
          ELSE
            l = i - ns
            Z(i,l) = Z(i,l) - ilambd
          END IF
        ELSE
          x = Rm1(mp,mp)/w
          y = t/w
          Rm1(mp,mp) = w
          IF ( i<N ) Rm1(mp,i+1) = 0.0E0
          IF ( i==N ) Z(mp,s-1) = 0.0E0
          !
          DO j = i, uk
            w = Rm1(j,i)
            Rm1(j,i) = Rm1(j,mp) - x*w
            Rm1(j,mp) = w
            IF ( j<n1 ) THEN
              Rm1(i,j+2) = Rm1(mp,j+2) - y*w
              Rm1(mp,j+2) = 0.0E0
            ELSE
              l = j - ns
              Z(i,l) = Z(mp,l) - y*w
              Z(mp,l) = 0.0E0
            END IF
          END DO
          !
          Rm1(i,i) = Rm1(i,i) - y*ilambd
          IF ( i<n1 ) THEN
            Rm1(mp,i+2) = -ilambd
            Rm1(i,i+2) = Rm1(i,i+2) + x*ilambd
          ELSE
            l = i - ns
            Z(mp,l) = -ilambd
            Z(i,l) = Z(i,l) + x*ilambd
          END IF
        END IF
      END DO
      !
      IF ( uk<n1 ) THEN
        t = Rm1(uk,uk+2)
      ELSE
        l = uk - ns
        t = Z(uk,l)
      END IF
      IF ( Rm1(uk,uk)==0.0E0.AND.t==0.0E0 ) Rm1(uk,uk) = eps3
      GOTO 250
    ELSE
      !     .......... REAL EIGENVALUE.
      !                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
      !                REPLACING ZERO PIVOTS BY EPS3 ..........
      IF ( uk/=1 ) THEN
        !
        DO i = 2, uk
          mp = i - 1
          IF ( ABS(Rm1(mp,i))>ABS(Rm1(mp,mp)) ) THEN
            !
            DO j = mp, uk
              y = Rm1(j,i)
              Rm1(j,i) = Rm1(j,mp)
              Rm1(j,mp) = y
            END DO
          END IF
          !
          IF ( Rm1(mp,mp)==0.0E0 ) Rm1(mp,mp) = eps3
          x = Rm1(mp,i)/Rm1(mp,mp)
          IF ( x/=0.0E0 ) THEN
            !
            DO j = i, uk
              Rm1(j,i) = Rm1(j,i) - x*Rm1(j,mp)
            END DO
          END IF
          !
        END DO
      END IF
      !
      IF ( Rm1(uk,uk)==0.0E0 ) Rm1(uk,uk) = eps3
    END IF
    !     .......... BACK SUBSTITUTION FOR REAL VECTOR
    !                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
    200 CONTINUE
    DO ii = 1, uk
      i = uk + 1 - ii
      y = Rv1(i)
      IF ( i/=uk ) THEN
        ip1 = i + 1
        !
        DO j = ip1, uk
          y = y - Rm1(j,i)*Rv1(j)
        END DO
      END IF
      !
      Rv1(i) = y/Rm1(i,i)
    END DO
    !
    GOTO 300
    !     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
    !                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
    250 CONTINUE
    DO ii = 1, uk
      i = uk + 1 - ii
      x = Rv1(i)
      y = 0.0E0
      IF ( i/=uk ) THEN
        ip1 = i + 1
        !
        DO j = ip1, uk
          IF ( j<n1 ) THEN
            t = Rm1(i,j+2)
          ELSE
            l = j - ns
            t = Z(i,l)
          END IF
          x = x - Rm1(j,i)*Rv1(j) + t*Rv2(j)
          y = y - Rm1(j,i)*Rv2(j) - t*Rv1(j)
        END DO
      END IF
      !
      IF ( i<n1 ) THEN
        t = Rm1(i,i+2)
      ELSE
        l = i - ns
        t = Z(i,l)
      END IF
      CALL CDIV(x,y,Rm1(i,i),t,Rv1(i),Rv2(i))
    END DO
    !     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
    !                EIGENVECTOR AND NORMALIZATION ..........
    300  its = its + 1
    norm = 0.0E0
    normv = 0.0E0
    !
    DO i = 1, uk
      IF ( ilambd==0.0E0 ) x = ABS(Rv1(i))
      IF ( ilambd/=0.0E0 ) x = PYTHAG(Rv1(i),Rv2(i))
      IF ( normv<x ) THEN
        normv = x
        j = i
      END IF
      norm = norm + x
    END DO
    !
    IF ( norm>=growto ) THEN
      !     .......... ACCEPT VECTOR ..........
      x = Rv1(j)
      IF ( ilambd==0.0E0 ) x = 1.0E0/x
      IF ( ilambd/=0.0E0 ) y = Rv2(j)
      !
      DO i = 1, uk
        IF ( ilambd/=0.0E0 ) THEN
          CALL CDIV(Rv1(i),Rv2(i),x,y,Z(i,s-1),Z(i,s))
        ELSE
          Z(i,s) = Rv1(i)*x
        END IF
      END DO
      !
      IF ( uk==N ) GOTO 350
      j = uk + 1
      !     .......... IN-LINE PROCEDURE FOR CHOOSING
      !                A NEW STARTING VECTOR ..........
    ELSEIF ( its>=uk ) THEN
      !     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
      j = 1
      Ierr = -k
    ELSE
      x = ukroot
      y = eps3/(x+1.0E0)
      Rv1(1) = eps3
      !
      DO i = 2, uk
        Rv1(i) = y
      END DO
      !
      j = uk - its + 1
      Rv1(j) = Rv1(j) - eps3*x
      IF ( ilambd/=0.0E0 ) GOTO 250
      GOTO 200
    END IF
    !     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
    DO i = j, N
      Z(i,s) = 0.0E0
      IF ( ilambd/=0.0E0 ) Z(i,s-1) = 0.0E0
    END DO
    !
    350  s = s + 1
    400  IF ( ip==(-1) ) ip = 0
    IF ( ip==1 ) ip = -1
  END DO
  !
  GOTO 600
  !     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
  !                SPACE REQUIRED ..........
  500 CONTINUE
  IF ( Ierr/=0 ) Ierr = Ierr - N
  IF ( Ierr==0 ) Ierr = -(2*N+1)
  600  M = s - 1 - ABS(ip)
END SUBROUTINE INVIT
