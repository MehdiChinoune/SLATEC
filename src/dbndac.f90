!*==DBNDAC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBNDAC
SUBROUTINE DBNDAC(G,Mdg,Nb,Ip,Ir,Mt,Jt)
  !***BEGIN PROLOGUE  DBNDAC
  !***PURPOSE  Compute the LU factorization of a  banded matrices using
  !            sequential accumulation of rows of the data matrix.
  !            Exactly one right-hand side vector is permitted.
  !***LIBRARY   SLATEC
  !***CATEGORY  D9
  !***TYPE      DOUBLE PRECISION (BNDACC-S, DBNDAC-D)
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
  !     The user of the subroutine DBNDAC provides MT,JT,C and F for
  !     T=1,...,Q.  Not all of this data must be supplied at once.
  !
  !     Following the processing of the various blocks (E F), the matrix
  !     (A B) has been transformed to the form (R D) where R is upper
  !     triangular and banded with bandwidth NB.  The least squares
  !     system Rx = d is then easily solved using back substitution by
  !     executing the statement CALL DBNDSL(1,...). The sequence of
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
  !     MU can be as small as N+2.  The subprogram DBNDAC processes the
  !     rows more efficiently if MU is large enough so that each new
  !     block (C F) has a distinct value of JT.
  !
  !     The four principle parts of these algorithms are obtained by the
  !     following call statements
  !
  !     CALL DBNDAC(...)  Introduce new blocks of data.
  !
  !     CALL DBNDSL(1,...)Compute solution vector and length of
  !                       residual vector.
  !
  !     CALL DBNDSL(2,...)Given any row vector H solve YR = H for the
  !                       row vector Y.
  !
  !     CALL DBNDSL(3,...)Given any column vector W solve RZ = W for
  !                       the column vector Z.
  !
  !     The dots in the above call statements indicate additional
  !     arguments that will be specified in the following paragraphs.
  !
  !     The user must dimension the array appearing in the call list..
  !     G(MDG,NB+1)
  !
  !     Description of calling sequence for DBNDAC..
  !
  !     The entire set of parameters for DBNDAC are
  !
  !     Input.. All Type REAL variables are DOUBLE PRECISION
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
  !                       first call to DBNDAC.  Its subsequent value
  !                       is controlled by DBNDAC to set up for the
  !                       next call to DBNDAC.
  !
  !     IR                Index of the row of G(*,*) where the user is
  !                       to place the new block of data (C F).  Set by
  !                       the user to the value 1 before the first call
  !                       to DBNDAC.  Its subsequent value is controlled
  !                       by DBNDAC. A value of IR .GT. MDG is considered
  !                       an error.
  !
  !     MT,JT             Set by the user to indicate respectively the
  !                       number of new rows of data in the block and
  !                       the index of the first nonzero column in that
  !                       set of rows (E F) = (0 C 0 F) being processed.
  !
  !     Output.. All Type REAL variables are DOUBLE PRECISION
  !
  !     G(*,*)            The working array which will contain the
  !                       processed rows of that part of the data
  !                       matrix which has been passed to DBNDAC.
  !
  !     IP,IR             The values of these arguments are advanced by
  !                       DBNDAC to be ready for storing and processing
  !                       a new block of data in G(*,*).
  !
  !     Description of calling sequence for DBNDSL..
  !
  !     The user must dimension the arrays appearing in the call list..
  !
  !     G(MDG,NB+1), X(N)
  !
  !     The entire set of parameters for DBNDSL are
  !
  !     Input.. All Type REAL variables are DOUBLE PRECISION
  !
  !     MODE              Set by the user to one of the values 1, 2, or
  !                       3.  These values respectively indicate that
  !                       the solution of AX = B, YR = H or RZ = W is
  !                       required.
  !
  !     G(*,*),MDG,       These arguments all have the same meaning and
  !      NB,IP,IR         contents as following the last call to DBNDAC.
  !
  !     X(*)              With mode=2 or 3 this array contains,
  !                       respectively, the right-side vectors H or W of
  !                       the systems YR = H or RZ = W.
  !
  !     N                 The number of variables in the solution
  !                       vector.  If any of the N diagonal terms are
  !                       zero the subroutine DBNDSL prints an
  !                       appropriate message.  This condition is
  !                       considered an error.
  !
  !     Output.. All Type REAL variables are DOUBLE PRECISION
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
  !     CALL DBNDAC(G,MDG,NB,IP,IR,MT,JT)
  !
  !***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
  !                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
  !***ROUTINES CALLED  DH12, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DBNDAC
  IMPLICIT NONE
  !*--DBNDAC193
  !*** Start of declarations inserted by SPAG
  REAL(8) :: G, rho, zero
  INTEGER i, ie, ig, ig1, ig2, iopt, Ip, Ir, j, jg, Jt, k, kh, &
    l, lp1, Mdg, mh, Mt, mu, Nb
  INTEGER nbp1, nerr
  !*** End of declarations inserted by SPAG
  DIMENSION G(Mdg,*)
  !***FIRST EXECUTABLE STATEMENT  DBNDAC
  zero = 0.D0
  !
  !              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
  !
  nbp1 = Nb + 1
  IF ( Mt<=0.OR.Nb<=0 ) RETURN
  !
  IF ( .NOT.Mdg<Ir ) THEN
    !
    !                                             ALG. STEP 5
    IF ( Jt/=Ip ) THEN
      !                                             ALG. STEPS 6-7
      IF ( Jt>Ir ) THEN
        !                                             ALG. STEPS 8-9
        DO i = 1, Mt
          ig1 = Jt + Mt - i
          ig2 = Ir + Mt - i
          DO j = 1, nbp1
            G(ig1,j) = G(ig2,j)
          ENDDO
        ENDDO
        !                                             ALG. STEP 10
        ie = Jt - Ir
        DO i = 1, ie
          ig = Ir + i - 1
          DO j = 1, nbp1
            G(ig,j) = zero
          ENDDO
        ENDDO
        !                                             ALG. STEP 11
        Ir = Jt
      ENDIF
      !                                             ALG. STEP 12
      mu = MIN(Nb-1,Ir-Ip-1)
      IF ( mu/=0 ) THEN
        !                                             ALG. STEP 13
        DO l = 1, mu
          !                                             ALG. STEP 14
          k = MIN(l,Jt-Ip)
          !                                             ALG. STEP 15
          lp1 = l + 1
          ig = Ip + l
          DO i = lp1, Nb
            jg = i - k
            G(ig,jg) = G(ig,i)
          ENDDO
          !                                             ALG. STEP 16
          DO i = 1, k
            jg = nbp1 - i
            G(ig,jg) = zero
          ENDDO
        ENDDO
      ENDIF
      !                                             ALG. STEP 17
      Ip = Jt
    ENDIF
    !                                             ALG. STEPS 18-19
    mh = Ir + Mt - Ip
    kh = MIN(nbp1,mh)
    !                                             ALG. STEP 20
    DO i = 1, kh
      CALL DH12(1,i,MAX(i+1,Ir-Ip+1),mh,G(Ip,i),1,rho,G(Ip,i+1),1,Mdg,&
        nbp1-i)
    ENDDO
    !                                             ALG. STEP 21
    Ir = Ip + kh
    !                                             ALG. STEP 22
    IF ( kh>=nbp1 ) THEN
      !                                             ALG. STEP 23
      DO i = 1, Nb
        G(Ir-1,i) = zero
      ENDDO
    ENDIF
  ELSE
    nerr = 1
    iopt = 2
    CALL XERMSG('SLATEC','DBNDAC','MDG.LT.IR, PROBABLE ERROR.',nerr,iopt)
    RETURN
  ENDIF
  !                                             ALG. STEP 24
  !                                             ALG. STEP 25
END SUBROUTINE DBNDAC
