!*==DGEFS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DGEFS
SUBROUTINE DGEFS(A,Lda,N,V,Itask,Ind,Work,Iwork)
  IMPLICIT NONE
  !*--DGEFS5
  !***BEGIN PROLOGUE  DGEFS
  !***PURPOSE  Solve a general system of linear equations.
  !***LIBRARY   SLATEC
  !***CATEGORY  D2A1
  !***TYPE      DOUBLE PRECISION (SGEFS-S, DGEFS-D, CGEFS-C)
  !***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
  !             GENERAL SYSTEM OF LINEAR EQUATIONS
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !    Subroutine DGEFS solves a general NxN system of double
  !    precision linear equations using LINPACK subroutines DGECO
  !    and DGESL.  That is, if A is an NxN double precision matrix
  !    and if X and B are double precision N-vectors, then DGEFS
  !    solves the equation
  !
  !                          A*X=B.
  !
  !    The matrix A is first factored into upper and lower tri-
  !    angular matrices U and L using partial pivoting.  These
  !    factors and the pivoting information are used to find the
  !    solution vector X.  An approximate condition number is
  !    calculated to provide a rough estimate of the number of
  !    digits of accuracy in the computed solution.
  !
  !    If the equation A*X=B is to be solved for more than one vector
  !    B, the factoring of A does not need to be performed again and
  !    the option to only solve (ITASK.GT.1) will be faster for
  !    the succeeding solutions.  In this case, the contents of A,
  !    LDA, N and IWORK must not have been altered by the user follow-
  !    ing factorization (ITASK=1).  IND will not be changed by DGEFS
  !    in this case.
  !
  !  Argument Description ***
  !
  !    A      DOUBLE PRECISION(LDA,N)
  !             on entry, the doubly subscripted array with dimension
  !               (LDA,N) which contains the coefficient matrix.
  !             on return, an upper triangular matrix U and the
  !               multipliers necessary to construct a matrix L
  !               so that A=L*U.
  !    LDA    INTEGER
  !             the leading dimension of the array A.  LDA must be great-
  !             er than or equal to N.  (terminal error message IND=-1)
  !    N      INTEGER
  !             the order of the matrix A.  The first N elements of
  !             the array A are the elements of the first column of
  !             the matrix A.  N must be greater than or equal to 1.
  !             (terminal error message IND=-2)
  !    V      DOUBLE PRECISION(N)
  !             on entry, the singly subscripted array(vector) of di-
  !               mension N which contains the right hand side B of a
  !               system of simultaneous linear equations A*X=B.
  !             on return, V contains the solution vector, X .
  !    ITASK  INTEGER
  !             If ITASK=1, the matrix A is factored and then the
  !               linear equation is solved.
  !             If ITASK .GT. 1, the equation is solved using the existing
  !               factored matrix A and IWORK.
  !             If ITASK .LT. 1, then terminal error message IND=-3 is
  !               printed.
  !    IND    INTEGER
  !             GT. 0  IND is a rough estimate of the number of digits
  !                     of accuracy in the solution, X.
  !             LT. 0  see error message corresponding to IND below.
  !    WORK   DOUBLE PRECISION(N)
  !             a singly subscripted array of dimension at least N.
  !    IWORK  INTEGER(N)
  !             a singly subscripted array of dimension at least N.
  !
  !  Error Messages Printed ***
  !
  !    IND=-1  terminal   N is greater than LDA.
  !    IND=-2  terminal   N is less than 1.
  !    IND=-3  terminal   ITASK is less than 1.
  !    IND=-4  terminal   The matrix A is computationally singular.
  !                         A solution has not been computed.
  !    IND=-10 warning    The solution has no apparent significance.
  !                         The solution may be inaccurate or the matrix
  !                         A may be poorly scaled.
  !
  !               Note-  The above terminal(*fatal*) error messages are
  !                      designed to be handled by XERMSG in which
  !                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
  !                      for warning error messages from XERMSG.  Unless
  !                      the user provides otherwise, an error message
  !                      will be printed followed by an abort.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  D1MACH, DGECO, DGESL, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800326  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGEFS
  !
  INTEGER Lda , N , Itask , Ind , Iwork(*)
  DOUBLE PRECISION A(Lda,*) , V(*) , Work(*) , D1MACH
  DOUBLE PRECISION rcond
  CHARACTER(8) :: xern1 , xern2
  !***FIRST EXECUTABLE STATEMENT  DGEFS
  IF ( Lda<N ) THEN
    Ind = -1
    WRITE (xern1,'(I8)') Lda
    WRITE (xern2,'(I8)') N
    CALL XERMSG('SLATEC','DGEFS','LDA = '//xern1//' IS LESS THAN N = '//&
      xern2,-1,1)
    RETURN
  ENDIF
  !
  IF ( N<=0 ) THEN
    Ind = -2
    WRITE (xern1,'(I8)') N
    CALL XERMSG('SLATEC','DGEFS','N = '//xern1//' IS LESS THAN 1',-2,1)
    RETURN
  ENDIF
  !
  IF ( Itask<1 ) THEN
    Ind = -3
    WRITE (xern1,'(I8)') Itask
    CALL XERMSG('SLATEC','DGEFS','ITASK = '//xern1//' IS LESS THAN 1',-3,1)
    RETURN
  ENDIF
  !
  IF ( Itask==1 ) THEN
    !
    !        FACTOR MATRIX A INTO LU
    !
    CALL DGECO(A,Lda,N,Iwork,rcond,Work)
    !
    !        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
    !
    IF ( rcond==0.0D0 ) THEN
      Ind = -4
      CALL XERMSG('SLATEC','DGEFS','SINGULAR MATRIX A - NO SOLUTION',-4,1)
      RETURN
    ENDIF
    !
    !        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
    !        AND CHECK FOR IND GREATER THAN ZERO
    !
    Ind = -LOG10(D1MACH(4)/rcond)
    IF ( Ind<=0 ) THEN
      Ind = -10
      CALL XERMSG('SLATEC','DGEFS','SOLUTION MAY HAVE NO SIGNIFICANCE',-10,&
        0)
    ENDIF
  ENDIF
  !
  !     SOLVE AFTER FACTORING
  !
  CALL DGESL(A,Lda,N,Iwork,V,0)
END SUBROUTINE DGEFS
