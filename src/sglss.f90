!DECK SGLSS
SUBROUTINE SGLSS(A,Mda,M,N,B,Mdb,Nb,Rnorm,Work,Lw,Iwork,Liw,Info)
  IMPLICIT NONE
  REAL A, ae, B, re, Rnorm, Work
  INTEGER Info, key, krank, ksure, Liw, Lw, M, Mda, Mdb, mode, N, &
    Nb, np
  !***BEGIN PROLOGUE  SGLSS
  !***PURPOSE  Solve a linear least squares problems by performing a QR
  !            factorization of the matrix using Householder
  !            transformations.  Emphasis is put on detecting possible
  !            rank deficiency.
  !***LIBRARY   SLATEC
  !***CATEGORY  D9, D5
  !***TYPE      SINGLE PRECISION (SGLSS-S, DGLSS-D)
  !***KEYWORDS  LINEAR LEAST SQUARES, LQ FACTORIZATION, QR FACTORIZATION,
  !             UNDERDETERMINED LINEAR SYSTEMS
  !***AUTHOR  Manteuffel, T. A., (LANL)
  !***DESCRIPTION
  !
  !     SGLSS solves both underdetermined and overdetermined
  !     LINEAR systems AX = B, where A is an M by N matrix
  !     and B is an M by NB matrix of right hand sides. If
  !     M.GE.N, the least squares solution is computed by
  !     decomposing the matrix A into the product of an
  !     orthogonal matrix Q and an upper triangular matrix
  !     R (QR factorization). If M.LT.N, the minimal
  !     length solution is computed by factoring the
  !     matrix A into the product of a lower triangular
  !     matrix L and an orthogonal matrix Q (LQ factor-
  !     ization). If the matrix A is determined to be rank
  !     deficient, that is the rank of A is less than
  !     MIN(M,N), then the minimal length least squares
  !     solution is computed.
  !
  !     SGLSS assumes full machine precision in the data.
  !     If more control over the uncertainty in the data
  !     is desired, the codes LLSIA and ULSIA are
  !     recommended.
  !
  !     SGLSS requires MDA*N + (MDB + 1)*NB + 5*MIN(M,N) dimensioned
  !     real space and M+N dimensioned integer space.
  !
  !
  !   ******************************************************************
  !   *                                                                *
  !   *         WARNING - All input arrays are changed on exit.        *
  !   *                                                                *
  !   ******************************************************************
  !     SUBROUTINE SGLSS(A,MDA,M,N,B,MDB,NB,RNORM,WORK,LW,IWORK,LIW,INFO)
  !
  !     Input..
  !
  !     A(,)          Linear coefficient matrix of AX=B, with MDA the
  !      MDA,M,N      actual first dimension of A in the calling program.
  !                   M is the row dimension (no. of EQUATIONS of the
  !                   problem) and N the col dimension (no. of UNKNOWNS).
  !
  !     B(,)          Right hand side(s), with MDB the actual first
  !      MDB,NB       dimension of B in the calling program. NB is the
  !                   number of M by 1 right hand sides. Must have
  !                   MDB.GE.MAX(M,N). If NB = 0, B is never accessed.
  !
  !
  !     RNORM()       Vector of length at least NB.  On input the contents
  !                   of RNORM are unused.
  !
  !     WORK()        A real work array dimensioned 5*MIN(M,N).
  !
  !     LW            Actual dimension of WORK.
  !
  !     IWORK()       Integer work array dimensioned at least N+M.
  !
  !     LIW           Actual dimension of IWORK.
  !
  !
  !     INFO          A flag which provides for the efficient
  !                   solution of subsequent problems involving the
  !                   same A but different B.
  !                   If INFO = 0 original call
  !                      INFO = 1 subsequent calls
  !                   On subsequent calls, the user must supply A, INFO,
  !                   LW, IWORK, LIW, and the first 2*MIN(M,N) locations
  !                   of WORK as output by the original call to SGLSS.
  !
  !
  !     Output..
  !
  !     A(,)          Contains the triangular part of the reduced matrix
  !                   and the transformation information. It together with
  !                   the first 2*MIN(M,N) elements of WORK (see below)
  !                   completely specify the factorization of A.
  !
  !     B(,)          Contains the N by NB solution matrix X.
  !
  !
  !     RNORM()       Contains the Euclidean length of the NB residual
  !                   vectors  B(I)-AX(I), I=1,NB.
  !
  !     WORK()        The first 2*MIN(M,N) locations of WORK contain value
  !                   necessary to reproduce the factorization of A.
  !
  !     IWORK()       The first M+N locations contain the order in
  !                   which the rows and columns of A were used.
  !                   If M.GE.N columns then rows. If M.LT.N rows
  !                   then columns.
  !
  !     INFO          Flag to indicate status of computation on completion
  !                  -1   Parameter error(s)
  !                   0 - Full rank
  !                   N.GT.0 - Reduced rank  rank=MIN(M,N)-INFO
  !
  !***REFERENCES  T. Manteuffel, An interval analysis approach to rank
  !                 determination in linear least squares problems,
  !                 Report SAND80-0655, Sandia Laboratories, June 1980.
  !***ROUTINES CALLED  LLSIA, ULSIA
  !***REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGLSS
  DIMENSION A(Mda,*), B(Mdb,*), Rnorm(*), Work(*)
  INTEGER Iwork(*)
  !
  !***FIRST EXECUTABLE STATEMENT  SGLSS
  re = 0.
  ae = 0.
  key = 0
  mode = 2
  np = 0
  !
  !     IF M.GE.N CALL LLSIA
  !     IF M.LT.N CALL ULSIA
  !
  IF ( M<N ) THEN
    CALL ULSIA(A,Mda,M,N,B,Mdb,Nb,[re],[ae],key,mode,np,krank,ksure,Rnorm,Work,&
      Lw,Iwork,Liw,Info)
    IF ( Info==-1 ) RETURN
    Info = M - krank
    GOTO 99999
  ENDIF
  CALL LLSIA(A,Mda,M,N,B,Mdb,Nb,[re],[ae],key,mode,np,krank,ksure,Rnorm,Work,Lw,&
    Iwork,Liw,Info)
  IF ( Info==-1 ) RETURN
  Info = N - krank
  RETURN
  99999 CONTINUE
  END SUBROUTINE SGLSS
