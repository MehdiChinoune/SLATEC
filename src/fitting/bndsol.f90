!DECK BNDSOL
SUBROUTINE BNDSOL(Mode,G,Mdg,Nb,Ip,Ir,X,N,Rnorm)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BNDSOL
  !***PURPOSE  Solve the least squares problem for a banded matrix using
  !            sequential accumulation of rows of the data matrix.
  !            Exactly one right-hand side vector is permitted.
  !***LIBRARY   SLATEC
  !***CATEGORY  D9
  !***TYPE      SINGLE PRECISION (BNDSOL-S, DBNDSL-D)
  !***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !     These subroutines solve the least squares problem Ax = b for
  !     banded matrices A using sequential accumulation of rows of the
  !     data matrix.  Exactly one right-hand side vector is permitted.
  !
  !     These subroutines are intended for the type of least squares
  !     systems that arise in applications such as curve or surface
  !     fitting of data.  The least squares equations are accumulated and
  !     processed using only part of the data.  This requires a certain
  !     user interaction during the solution of Ax = b.
  !
  !     Specifically, suppose the data matrix (A B) is row partitioned
  !     into Q submatrices.  Let (E F) be the T-th one of these
  !     submatrices where E = (0 C 0).  Here the dimension of E is MT by N
  !     and the dimension of C is MT by NB.  The value of NB is the
  !     bandwidth of A.  The dimensions of the leading block of zeros in E
  !     are MT by JT-1.
  !
  !     The user of the subroutine BNDACC provides MT,JT,C and F for
  !     T=1,...,Q.  Not all of this data must be supplied at once.
  !
  !     Following the processing of the various blocks (E F), the matrix
  !     (A B) has been transformed to the form (R D) where R is upper
  !     triangular and banded with bandwidth NB.  The least squares
  !     system Rx = d is then easily solved using back substitution by
  !     executing the statement CALL BNDSOL(1,...). The sequence of
  !     values for JT must be nondecreasing.  This may require some
  !     preliminary interchanges of rows and columns of the matrix A.
  !
  !     The primary reason for these subroutines is that the total
  !     processing can take place in a working array of dimension MU by
  !     NB+1.  An acceptable value for MU is
  !
  !                       MU = MAX(MT + N + 1),
  !
  !     where N is the number of unknowns.
  !
  !     Here the maximum is taken over all values of MT for T=1,...,Q.
  !     Notice that MT can be taken to be a small as one, showing that
  !     MU can be as small as N+2.  The subprogram BNDACC processes the
  !     rows more efficiently if MU is large enough so that each new
  !     block (C F) has a distinct value of JT.
  !
  !     The four principle parts of these algorithms are obtained by the
  !     following call statements
  !
  !     CALL BNDACC(...)  Introduce new blocks of data.
  !
  !     CALL BNDSOL(1,...)Compute solution vector and length of
  !                       residual vector.
  !
  !     CALL BNDSOL(2,...)Given any row vector H solve YR = H for the
  !                       row vector Y.
  !
  !     CALL BNDSOL(3,...)Given any column vector W solve RZ = W for
  !                       the column vector Z.
  !
  !     The dots in the above call statements indicate additional
  !     arguments that will be specified in the following paragraphs.
  !
  !     The user must dimension the array appearing in the call list..
  !     G(MDG,NB+1)
  !
  !     Description of calling sequence for BNDACC..
  !
  !     The entire set of parameters for BNDACC are
  !
  !     Input..
  !
  !     G(*,*)            The working array into which the user will
  !                       place the MT by NB+1 block (C F) in rows IR
  !                       through IR+MT-1, columns 1 through NB+1.
  !                       See descriptions of IR and MT below.
  !
  !     MDG               The number of rows in the working array
  !                       G(*,*).  The value of MDG should be .GE. MU.
  !                       The value of MU is defined in the abstract
  !                       of these subprograms.
  !
  !     NB                The bandwidth of the data matrix A.
  !
  !     IP                Set by the user to the value 1 before the
  !                       first call to BNDACC.  Its subsequent value
  !                       is controlled by BNDACC to set up for the
  !                       next call to BNDACC.
  !
  !     IR                Index of the row of G(*,*) where the user is
  !                       the user to the value 1 before the first call
  !                       to BNDACC.  Its subsequent value is controlled
  !                       by BNDACC. A value of IR .GT. MDG is considered
  !                       an error.
  !
  !     MT,JT             Set by the user to indicate respectively the
  !                       number of new rows of data in the block and
  !                       the index of the first nonzero column in that
  !                       set of rows (E F) = (0 C 0 F) being processed.
  !     Output..
  !
  !     G(*,*)            The working array which will contain the
  !                       processed rows of that part of the data
  !                       matrix which has been passed to BNDACC.
  !
  !     IP,IR             The values of these arguments are advanced by
  !                       BNDACC to be ready for storing and processing
  !                       a new block of data in G(*,*).
  !
  !     Description of calling sequence for BNDSOL..
  !
  !     The user must dimension the arrays appearing in the call list..
  !
  !     G(MDG,NB+1), X(N)
  !
  !     The entire set of parameters for BNDSOL are
  !
  !     Input..
  !
  !     MODE              Set by the user to one of the values 1, 2, or
  !                       3.  These values respectively indicate that
  !                       the solution of AX = B, YR = H or RZ = W is
  !                       required.
  !
  !     G(*,*),MDG,       These arguments all have the same meaning and
  !      NB,IP,IR         contents as following the last call to BNDACC.
  !
  !     X(*)              With mode=2 or 3 this array contains,
  !                       respectively, the right-side vectors H or W of
  !                       the systems YR = H or RZ = W.
  !
  !     N                 The number of variables in the solution
  !                       vector.  If any of the N diagonal terms are
  !                       zero the subroutine BNDSOL prints an
  !                       appropriate message.  This condition is
  !                       considered an error.
  !
  !     Output..
  !
  !     X(*)              This array contains the solution vectors X,
  !                       Y or Z of the systems AX = B, YR = H or
  !                       RZ = W depending on the value of MODE=1,
  !                       2 or 3.
  !
  !     RNORM             If MODE=1 RNORM is the Euclidean length of the
  !                       residual vector AX-B.  When MODE=2 or 3 RNORM
  !                       is set to zero.
  !
  !     Remarks..
  !
  !     To obtain the upper triangular matrix and transformed right-hand
  !     side vector D so that the super diagonals of R form the columns
  !     of G(*,*), execute the following Fortran statements.
  !
  !     NBP1=NB+1
  !
  !     DO 10 J=1, NBP1
  !
  !  10 G(IR,J) = 0.E0
  !
  !     MT=1
  !
  !     JT=N+1
  !
  !     CALL BNDACC(G,MDG,NB,IP,IR,MT,JT)
  !
  !***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
  !                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BNDSOL
  REAL G, Rnorm, rsq, s, X, zero
  INTEGER i, i1, i2, ie, ii, iopt, Ip, Ir, irm1, ix, j, jg, l, &
    Mdg, Mode, N, Nb, nerr, np1
  DIMENSION G(Mdg,*), X(*)
  !***FIRST EXECUTABLE STATEMENT  BNDSOL
  zero = 0.
  !
  Rnorm = zero
  SELECT CASE (Mode)
    CASE (2)
      !                                   ********************* MODE = 2
      DO j = 1, N
        s = zero
        IF ( j/=1 ) THEN
          i1 = MAX(1,j-Nb+1)
          i2 = j - 1
          DO i = i1, i2
            l = j - i + 1 + MAX(0,i-Ip)
            s = s + X(i)*G(i,l)
          ENDDO
        ENDIF
        l = MAX(0,j-Ip)
        IF ( G(j,l+1)==0 ) GOTO 100
        X(j) = (X(j)-s)/G(j,l+1)
      ENDDO
      RETURN
    CASE (3)
    CASE DEFAULT
      !                                   ********************* MODE = 1
      !                                   ALG. STEP 26
      DO j = 1, N
        X(j) = G(j,Nb+1)
      ENDDO
      rsq = zero
      np1 = N + 1
      irm1 = Ir - 1
      IF ( np1<=irm1 ) THEN
        DO j = np1, irm1
          rsq = rsq + G(j,Nb+1)**2
        ENDDO
        Rnorm = SQRT(rsq)
      ENDIF
  END SELECT
  !                                   ********************* MODE = 3
  !                                   ALG. STEP 27
  DO ii = 1, N
    i = N + 1 - ii
    !                                   ALG. STEP 28
    s = zero
    l = MAX(0,i-Ip)
    !                                   ALG. STEP 29
    IF ( i/=N ) THEN
      !                                   ALG. STEP 30
      ie = MIN(N+1-i,Nb)
      DO j = 2, ie
        jg = j + l
        ix = i - 1 + j
        s = s + G(i,jg)*X(ix)
      ENDDO
    ENDIF
    !                                   ALG. STEP 31
    IF ( G(i,l+1)==0 ) GOTO 100
    X(i) = (X(i)-s)/G(i,l+1)
  ENDDO
  !                                   ALG. STEP 32
  RETURN
  !
  100  nerr = 1
  iopt = 2
  CALL XERMSG('SLATEC','BNDSOL',&
    'A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR MATRIX.',nerr,iopt)
END SUBROUTINE BNDSOL
