!** CINVIT
SUBROUTINE CINVIT(Nm,N,Ar,Ai,Wr,Wi,Select,Mm,M,Zr,Zi,Ierr,Rm1,Rm2,Rv1,Rv2)
  IMPLICIT NONE
  !>
  !***
  !  Compute the eigenvectors of a complex upper Hessenberg
  !            associated with specified eigenvalues using inverse
  !            iteration.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C2B
  !***
  ! **Type:**      COMPLEX (INVIT-S, CINVIT-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure CXINVIT
  !     by Peters and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
  !
  !     This subroutine finds those eigenvectors of A COMPLEX UPPER
  !     Hessenberg matrix corresponding to specified eigenvalues,
  !     using inverse iteration.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR, AI, ZR and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        N is the order of the matrix A=(AR,AI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the complex upper Hessenberg matrix.  AR and AI are
  !          two-dimensional REAL arrays, dimensioned AR(NM,N)
  !          and AI(NM,N).
  !
  !        WR and WI contain the real and imaginary parts, respectively,
  !          of the eigenvalues of the matrix.  The eigenvalues must be
  !          stored in a manner identical to that of subroutine  COMLR,
  !          which recognizes possible splitting of the matrix.  WR and
  !          WI are one-dimensional REAL arrays, dimensioned WR(N) and
  !          WI(N).
  !
  !        SELECT specifies the eigenvectors to be found.  The
  !          eigenvector corresponding to the J-th eigenvalue is
  !          specified by setting SELECT(J) to .TRUE.  SELECT is a
  !          one-dimensional LOGICAL array, dimensioned SELECT(N).
  !
  !        MM should be set to an upper bound for the number of
  !          eigenvectors to be found.  MM is an INTEGER variable.
  !
  !     On OUTPUT
  !
  !        AR, AI, WI, and SELECT are unaltered.
  !
  !        WR may have been altered since close eigenvalues are perturbed
  !          slightly in searching for independent eigenvectors.
  !
  !        M is the number of eigenvectors actually found.  M is an
  !          INTEGER variable.
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors corresponding to the flagged eigenvalues.
  !          The eigenvectors are normalized so that the component of
  !          largest magnitude is 1.  Any vector which fails the
  !          acceptance test is set to zero.  ZR and ZI are
  !          two-dimensional REAL arrays, dimensioned ZR(NM,MM) and
  !          ZI(NM,MM).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          -(2*N+1)   if more than MM eigenvectors have been requested
  !                     (the MM eigenvectors calculated to this point are
  !                     in ZR and ZI),
  !          -K         if the iteration corresponding to the K-th
  !                     value fails (if this occurs more than once, K
  !                     is the index of the last occurrence); the
  !                     corresponding columns of ZR and ZI are set to
  !                     zero vectors,
  !          -(N+K)     if both error situations occur.
  !
  !        RV1 and RV2 are one-dimensional REAL arrays used for
  !          temporary storage, dimensioned RV1(N) and RV2(N).
  !          They hold the approximate eigenvectors during the inverse
  !          iteration process.
  !
  !        RM1 and RM2 are two-dimensional REAL arrays used for
  !          temporary storage, dimensioned RM1(N,N) and RM2(N,N).
  !          These arrays hold the triangularized form of the upper
  !          Hessenberg matrix used in the inverse iteration process.
  !
  !     The ALGOL procedure GUESSVEC appears in CINVIT in-line.
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
  INTEGER i, j, k, M, N, s, ii, Mm, mp, Nm, uk, ip1, its, km1, Ierr
  REAL Ar(Nm,*), Ai(Nm,*), Wr(*), Wi(*), Zr(Nm,*), Zi(Nm,*)
  REAL Rm1(N,*), Rm2(N,*), Rv1(*), Rv2(*)
  REAL x, y, eps3, norm, normv, growto, ilambd, rlambd, ukroot
  REAL PYTHAG
  LOGICAL Select(N)
  !
  !* FIRST EXECUTABLE STATEMENT  CINVIT
  Ierr = 0
  uk = 0
  s = 1
  !
  DO k = 1, N
    IF ( .NOT.Select(k) ) CYCLE
    IF ( s>Mm ) GOTO 300
    IF ( uk<k ) THEN
      !     .......... CHECK FOR POSSIBLE SPLITTING ..........
      DO uk = k, N
        IF ( uk==N ) EXIT
        IF ( Ar(uk+1,uk)==0.0E0.AND.Ai(uk+1,uk)==0.0E0 ) EXIT
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
          x = x + PYTHAG(Ar(i,j),Ai(i,j))
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
          !     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
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
    !     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
    !                AND INITIAL COMPLEX VECTOR ..........
    150  mp = 1
    !
    DO i = 1, uk
      !
      DO j = mp, uk
        Rm1(i,j) = Ar(i,j)
        Rm2(i,j) = Ai(i,j)
      END DO
      !
      Rm1(i,i) = Rm1(i,i) - rlambd
      Rm2(i,i) = Rm2(i,i) - ilambd
      mp = i
      Rv1(i) = eps3
    END DO
    !     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
    !                REPLACING ZERO PIVOTS BY EPS3 ..........
    IF ( uk/=1 ) THEN
      !
      DO i = 2, uk
        mp = i - 1
        IF ( PYTHAG(Rm1(i,mp),Rm2(i,mp))>PYTHAG(Rm1(mp,mp),Rm2(mp,mp)) ) THEN
          !
          DO j = mp, uk
            y = Rm1(i,j)
            Rm1(i,j) = Rm1(mp,j)
            Rm1(mp,j) = y
            y = Rm2(i,j)
            Rm2(i,j) = Rm2(mp,j)
            Rm2(mp,j) = y
          END DO
        END IF
        !
        IF ( Rm1(mp,mp)==0.0E0.AND.Rm2(mp,mp)==0.0E0 ) Rm1(mp,mp) = eps3
        CALL CDIV(Rm1(i,mp),Rm2(i,mp),Rm1(mp,mp),Rm2(mp,mp),x,y)
        IF ( x/=0.0E0.OR.y/=0.0E0 ) THEN
          !
          DO j = i, uk
            Rm1(i,j) = Rm1(i,j) - x*Rm1(mp,j) + y*Rm2(mp,j)
            Rm2(i,j) = Rm2(i,j) - x*Rm2(mp,j) - y*Rm1(mp,j)
          END DO
        END IF
        !
      END DO
    END IF
    !
    IF ( Rm1(uk,uk)==0.0E0.AND.Rm2(uk,uk)==0.0E0 ) Rm1(uk,uk) = eps3
    its = 0
    DO
      !     .......... BACK SUBSTITUTION
      !                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
      DO ii = 1, uk
        i = uk + 1 - ii
        x = Rv1(i)
        y = 0.0E0
        IF ( i/=uk ) THEN
          ip1 = i + 1
          !
          DO j = ip1, uk
            x = x - Rm1(i,j)*Rv1(j) + Rm2(i,j)*Rv2(j)
            y = y - Rm1(i,j)*Rv2(j) - Rm2(i,j)*Rv1(j)
          END DO
        END IF
        !
        CALL CDIV(x,y,Rm1(i,i),Rm2(i,i),Rv1(i),Rv2(i))
      END DO
      !     .......... ACCEPTANCE TEST FOR EIGENVECTOR
      !                AND NORMALIZATION ..........
      its = its + 1
      norm = 0.0E0
      normv = 0.0E0
      !
      DO i = 1, uk
        x = PYTHAG(Rv1(i),Rv2(i))
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
        y = Rv2(j)
        !
        DO i = 1, uk
          CALL CDIV(Rv1(i),Rv2(i),x,y,Zr(i,s),Zi(i,s))
        END DO
        !
        IF ( uk==N ) GOTO 200
        j = uk + 1
        EXIT
        !     .......... IN-LINE PROCEDURE FOR CHOOSING
        !                A NEW STARTING VECTOR ..........
      ELSEIF ( its>=uk ) THEN
        !     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
        j = 1
        Ierr = -k
        EXIT
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
      END IF
    END DO
    !     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
    DO i = j, N
      Zr(i,s) = 0.0E0
      Zi(i,s) = 0.0E0
    END DO
    !
    200  s = s + 1
  END DO
  !
  GOTO 400
  !     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
  !                SPACE REQUIRED ..........
  300 CONTINUE
  IF ( Ierr/=0 ) Ierr = Ierr - N
  IF ( Ierr==0 ) Ierr = -(2*N+1)
  400  M = s - 1
END SUBROUTINE CINVIT
