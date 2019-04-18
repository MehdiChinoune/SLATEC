!** DHFTI
SUBROUTINE DHFTI(A,Mda,M,N,B,Mdb,Nb,Tau,Krank,Rnorm,H,G,Ip)
  !>
  !***
  !  Solve a least squares problem for banded matrices using
  !            sequential accumulation of rows of the data matrix.
  !            Exactly one right-hand side vector is permitted.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D9
  !***
  ! **Type:**      DOUBLE PRECISION (HFTI-S, DHFTI-D)
  !***
  ! **Keywords:**  CURVE FITTING, LEAST SQUARES
  !***
  ! **Author:**  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
  !
  !     This subroutine solves a linear least squares problem or a set of
  !     linear least squares problems having the same matrix but different
  !     right-side vectors.  The problem data consists of an M by N matrix
  !     A, an M by NB matrix B, and an absolute tolerance parameter TAU
  !     whose usage is described below.  The NB column vectors of B
  !     represent right-side vectors for NB distinct linear least squares
  !     problems.
  !
  !     This set of problems can also be written as the matrix least
  !     squares problem
  !
  !                       AX = B,
  !
  !     where X is the N by NB solution matrix.
  !
  !     Note that if B is the M by M identity matrix, then X will be the
  !     pseudo-inverse of A.
  !
  !     This subroutine first transforms the augmented matrix (A B) to a
  !     matrix (R C) using premultiplying Householder transformations with
  !     column interchanges.  All subdiagonal elements in the matrix R are
  !     zero and its diagonal elements satisfy
  !
  !                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)),
  !
  !                       I = 1,...,L-1, where
  !
  !                       L = MIN(M,N).
  !
  !     The subroutine will compute an integer, KRANK, equal to the number
  !     of diagonal terms of R that exceed TAU in magnitude. Then a
  !     solution of minimum Euclidean length is computed using the first
  !     KRANK rows of (R C).
  !
  !     To be specific we suggest that the user consider an easily
  !     computable matrix norm, such as, the maximum of all column sums of
  !     magnitudes.
  !
  !     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
  !     norm of B), it is suggested that TAU be set approximately equal to
  !     EPS*(norm of A).
  !
  !     The user must dimension all arrays appearing in the call list..
  !     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
  !     permits the solution of a range of problems in the same array
  !     space.
  !
  !     The entire set of parameters for DHFTI are
  !
  !     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
  !
  !     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
  !                       matrix A of the least squares problem AX = B.
  !                       The first dimensioning parameter of the array
  !                       A(*,*) is MDA, which must satisfy MDA.GE.M
  !                       Either M.GE.N or M.LT.N is permitted.  There
  !                       is no restriction on the rank of A.  The
  !                       condition MDA.LT.M is considered an error.
  !
  !     B(*),MDB,NB       If NB = 0 the subroutine will perform the
  !                       orthogonal decomposition but will make no
  !                       references to the array B(*).  If NB.GT.0
  !                       the array B(*) must initially contain the M by
  !                       NB matrix B of the least squares problem AX =
  !                       B.  If NB.GE.2 the array B(*) must be doubly
  !                       subscripted with first dimensioning parameter
  !                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
  !                       be either doubly or singly subscripted.  In
  !                       the latter case the value of MDB is arbitrary
  !                       but it should be set to some valid integer
  !                       value such as MDB = M.
  !
  !                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
  !                       is considered an error.
  !
  !     TAU               Absolute tolerance parameter provided by user
  !                       for pseudorank determination.
  !
  !     H(*),G(*),IP(*)   Arrays of working space used by DHFTI.
  !
  !     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
  !
  !     A(*,*)            The contents of the array A(*,*) will be
  !                       modified by the subroutine. These contents
  !                       are not generally required by the user.
  !
  !     B(*)              On return the array B(*) will contain the N by
  !                       NB solution matrix X.
  !
  !     KRANK             Set by the subroutine to indicate the
  !                       pseudorank of A.
  !
  !     RNORM(*)          On return, RNORM(J) will contain the Euclidean
  !                       norm of the residual vector for the problem
  !                       defined by the J-th column vector of the array
  !                       B(*,*) for J = 1,...,NB.
  !
  !     H(*),G(*)         On return these arrays respectively contain
  !                       elements of the pre- and post-multiplying
  !                       Householder transformations used to compute
  !                       the minimum Euclidean length solution.
  !
  !     IP(*)             Array in which the subroutine records indices
  !                       describing the permutation of column vectors.
  !                       The contents of arrays H(*),G(*) and IP(*)
  !                       are not generally required by the user.
  !
  !***
  ! **References:**  C. L. Lawson and R. J. Hanson, Solving Least Squares
  !                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
  !***
  ! **Routines called:**  D1MACH, DH12, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : D1MACH, XERMSG
  INTEGER i, ii, iopt, Ip(*), ip1, j, jb, jj, k, kp1, Krank, l, &
    ldiag, lmax, M, Mda, Mdb, N, Nb, nerr
  REAL(8) :: A(Mda,*), B(Mdb,*), dzero, factor, G(*), H(*), hmax, &
    Rnorm(*), sm, sm1, szero, Tau, tmp
  REAL(8) :: releps = 0.D0
  !     BEGIN BLOCK PERMITTING ...EXITS TO 360
  !* FIRST EXECUTABLE STATEMENT  DHFTI
  IF ( releps==0.D0 ) releps = D1MACH(4)
  szero = 0.0D0
  dzero = 0.0D0
  factor = 0.001D0
  !
  k = 0
  ldiag = MIN(M,N)
  IF ( ldiag>0 ) THEN
    !           BEGIN BLOCK PERMITTING ...EXITS TO 130
    !              BEGIN BLOCK PERMITTING ...EXITS TO 120
    IF ( Mda>=M ) THEN
      !
      IF ( Nb<=1.OR.MAX(M,N)<=Mdb ) THEN
        !
        DO j = 1, ldiag
          !                    BEGIN BLOCK PERMITTING ...EXITS TO 70
          IF ( j/=1 ) THEN
            !
            !                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
            !                          ..
            lmax = j
            DO l = j, N
              H(l) = H(l) - A(j-1,l)**2
              IF ( H(l)>H(lmax) ) lmax = l
            END DO
            !                    ......EXIT
            IF ( factor*H(lmax)>hmax*releps ) GOTO 5
          END IF
          !
          !                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
          !                       ..
          lmax = j
          DO l = j, N
            H(l) = 0.0D0
            DO i = j, M
              H(l) = H(l) + A(i,l)**2
            END DO
            IF ( H(l)>H(lmax) ) lmax = l
          END DO
          hmax = H(lmax)
          !                    ..
          !                     LMAX HAS BEEN DETERMINED
          !
          !                     DO COLUMN INTERCHANGES IF NEEDED.
          !                    ..
          5  Ip(j) = lmax
          IF ( Ip(j)/=j ) THEN
            DO i = 1, M
              tmp = A(i,j)
              A(i,j) = A(i,lmax)
              A(i,lmax) = tmp
            END DO
            H(lmax) = H(j)
          END IF
          !
          !                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
          !                     AND B.
          !                    ..
          CALL DH12(1,j,j+1,M,A(1,j),1,H(j),A(1,j+1),1,Mda,N-j)
          CALL DH12(2,j,j+1,M,A(1,j),1,H(j),B,1,Mdb,Nb)
        END DO
        !
        !                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
        !                  TAU.
        !                 ..
        DO j = 1, ldiag
          !              ......EXIT
          IF ( ABS(A(j,j))<=Tau ) GOTO 20
        END DO
        k = ldiag
        !           ......EXIT
        GOTO 50
      ELSE
        nerr = 2
        iopt = 2
        CALL XERMSG('SLATEC','DHFTI',&
          'MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.',nerr,iopt)
        !     ...............EXIT
        RETURN
      END IF
      20  k = j - 1
    ELSE
      nerr = 1
      iopt = 2
      CALL XERMSG('SLATEC','DHFTI','MDA.LT.M, PROBABLE ERROR.',nerr,iopt)
      !     ...............EXIT
      RETURN
    END IF
    50  kp1 = k + 1
    !
    !           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
    !
    IF ( Nb>=1 ) THEN
      DO jb = 1, Nb
        tmp = szero
        IF ( M>=kp1 ) THEN
          DO i = kp1, M
            tmp = tmp + B(i,jb)**2
          END DO
        END IF
        Rnorm(jb) = SQRT(tmp)
      END DO
    END IF
    !           SPECIAL FOR PSEUDORANK = 0
    IF ( k>0 ) THEN
      !
      !               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
      !               DECOMPOSITION OF FIRST K ROWS.
      !              ..
      IF ( k/=N ) THEN
        DO ii = 1, k
          i = kp1 - ii
          CALL DH12(1,i,kp1,N,A(i,1),Mda,G(i),A,Mda,1,i-1)
        END DO
      END IF
      !
      !
      IF ( Nb>=1 ) THEN
        DO jb = 1, Nb
          !
          !                  SOLVE THE K BY K TRIANGULAR SYSTEM.
          !                 ..
          DO l = 1, k
            sm = dzero
            i = kp1 - l
            ip1 = i + 1
            IF ( k>=ip1 ) THEN
              DO j = ip1, k
                sm = sm + A(i,j)*B(j,jb)
              END DO
            END IF
            sm1 = sm
            B(i,jb) = (B(i,jb)-sm1)/A(i,i)
          END DO
          !
          !                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
          !                 ..
          IF ( k/=N ) THEN
            DO j = kp1, N
              B(j,jb) = szero
            END DO
            DO i = 1, k
              CALL DH12(2,i,kp1,N,A(i,1),Mda,G(i),B(1,jb),1,Mdb,1)
            END DO
          END IF
          !
          !                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
          !                   COLUMN INTERCHANGES.
          !                 ..
          DO jj = 1, ldiag
            j = ldiag + 1 - jj
            IF ( Ip(j)/=j ) THEN
              l = Ip(j)
              tmp = B(l,jb)
              B(l,jb) = B(j,jb)
              B(j,jb) = tmp
            END IF
          END DO
        END DO
      END IF
    ELSEIF ( Nb>=1 ) THEN
      DO jb = 1, Nb
        DO i = 1, N
          B(i,jb) = szero
        END DO
      END DO
    END IF
  END IF
  !        ..
  !         THE SOLUTION VECTORS, X, ARE NOW
  !         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
  !
  Krank = k
  RETURN
END SUBROUTINE DHFTI
