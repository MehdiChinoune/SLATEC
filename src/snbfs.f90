!*==SNBFS.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SNBFS
SUBROUTINE SNBFS(Abe,Lda,N,Ml,Mu,V,Itask,Ind,Work,Iwork)
  IMPLICIT NONE
  !*--SNBFS5
  !***BEGIN PROLOGUE  SNBFS
  !***PURPOSE  Solve a general nonsymmetric banded system of linear
  !            equations.
  !***LIBRARY   SLATEC
  !***CATEGORY  D2A2
  !***TYPE      SINGLE PRECISION (SNBFS-S, DNBFS-D, CNBFS-C)
  !***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC
  !***AUTHOR  Voorhees, E. A., (LANL)
  !***DESCRIPTION
  !
  !    Subroutine SNBFS solves a general nonsymmetric banded NxN
  !    system of single precision real linear equations using
  !    SLATEC subroutines SNBCO and SNBSL.  These are adaptations
  !    of the LINPACK subroutines SBGCO and SGBSL, which require
  !    a different format for storing the matrix elements.  If
  !    A  is an NxN real matrix and if  X  and  B  are real
  !    N-vectors, then SNBFS solves the equation
  !
  !                          A*X=B.
  !
  !    A band matrix is a matrix whose nonzero elements are all
  !    fairly near the main diagonal, specifically  A(I,J) = 0
  !    if  I-J is greater than  ML  or  J-I  is greater than
  !    MU .  The integers ML and MU are called the lower and upper
  !    band widths and  M = ML+MU+1  is the total band width.
  !    SNBFS uses less time and storage than the corresponding
  !    program for general matrices (SGEFS) if 2*ML+MU .LT. N  .
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
  !    the option to only solve (ITASK .GT. 1) will be faster for
  !    the succeeding solutions.  In this case, the contents of A,
  !    LDA, N and IWORK must not have been altered by the user follow-
  !    ing factorization (ITASK=1).  IND will not be changed by SNBFS
  !    in this case.
  !
  !
  !    Band Storage
  !
  !          If  A  is a band matrix, the following program segment
  !          will set up the input.
  !
  !                  ML = (band width below the diagonal)
  !                  MU = (band width above the diagonal)
  !                  DO 20 I = 1, N
  !                     J1 = MAX(1, I-ML)
  !                     J2 = MIN(N, I+MU)
  !                     DO 10 J = J1, J2
  !                        K = J - I + ML + 1
  !                        ABE(I,K) = A(I,J)
  !               10    CONTINUE
  !               20 CONTINUE
  !
  !          This uses columns  1  through  ML+MU+1  of ABE .
  !          Furthermore,  ML  additional columns are needed in
  !          ABE  starting with column  ML+MU+2  for elements
  !          generated during the triangularization.  The total
  !          number of columns needed in  ABE  is  2*ML+MU+1 .
  !
  !    Example:  If the original matrix is
  !
  !          11 12 13  0  0  0
  !          21 22 23 24  0  0
  !           0 32 33 34 35  0
  !           0  0 43 44 45 46
  !           0  0  0 54 55 56
  !           0  0  0  0 65 66
  !
  !     then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain
  !
  !           * 11 12 13  +     , * = not used
  !          21 22 23 24  +     , + = used for pivoting
  !          32 33 34 35  +
  !          43 44 45 46  +
  !          54 55 56  *  +
  !          65 66  *  *  +
  !
  !
  !  Argument Description ***
  !
  !    ABE    REAL(LDA,NC)
  !             on entry, contains the matrix in band storage as
  !               described above.  NC  must not be less than
  !               2*ML+MU+1 .  The user is cautioned to specify  NC
  !               with care since it is not an argument and cannot
  !               be checked by SNBFS.  The rows of the original
  !               matrix are stored in the rows of  ABE  and the
  !               diagonals of the original matrix are stored in
  !               columns  1  through  ML+MU+1  of  ABE .
  !             on return, contains an upper triangular matrix U and
  !               the multipliers necessary to construct a matrix L
  !               so that A=L*U.
  !    LDA    INTEGER
  !             the leading dimension of array ABE.  LDA must be great-
  !             er than or equal to N.  (terminal error message IND=-1)
  !    N      INTEGER
  !             the order of the matrix A.  N must be greater
  !             than or equal to 1 .  (terminal error message IND=-2)
  !    ML     INTEGER
  !             the number of diagonals below the main diagonal.
  !             ML  must not be less than zero nor greater than or
  !             equal to  N .  (terminal error message IND=-5)
  !    MU     INTEGER
  !             the number of diagonals above the main diagonal.
  !             MU  must not be less than zero nor greater than or
  !             equal to  N .  (terminal error message IND=-6)
  !    V      REAL(N)
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
  !             LT. 0  See error message corresponding to IND below.
  !    WORK   REAL(N)
  !             a singly subscripted array of dimension at least N.
  !    IWORK  INTEGER(N)
  !             a singly subscripted array of dimension at least N.
  !
  !  Error Messages Printed ***
  !
  !    IND=-1  terminal   N is greater than LDA.
  !    IND=-2  terminal   N is less than 1.
  !    IND=-3  terminal   ITASK is less than 1.
  !    IND=-4  terminal   the matrix A is computationally singular.
  !                         A solution has not been computed.
  !    IND=-5  terminal   ML is less than zero or is greater than
  !                         or equal to N .
  !    IND=-6  terminal   MU is less than zero or is greater than
  !                         or equal to N .
  !    IND=-10 warning    the solution has no apparent significance.
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
  !***ROUTINES CALLED  R1MACH, SNBCO, SNBSL, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SNBFS
  !
  INTEGER Lda , N , Itask , Ind , Iwork(*) , Ml , Mu
  REAL Abe(Lda,*) , V(*) , Work(*) , R1MACH
  REAL rcond
  CHARACTER(8) :: xern1 , xern2
  !***FIRST EXECUTABLE STATEMENT  SNBFS
  IF ( Lda<N ) THEN
    Ind = -1
    WRITE (xern1,'(I8)') Lda
    WRITE (xern2,'(I8)') N
    CALL XERMSG('SLATEC','SNBFS','LDA = '//xern1//' IS LESS THAN N = '//&
      xern2,-1,1)
    RETURN
  ENDIF
  !
  IF ( N<=0 ) THEN
    Ind = -2
    WRITE (xern1,'(I8)') N
    CALL XERMSG('SLATEC','SNBFS','N = '//xern1//' IS LESS THAN 1',-2,1)
    RETURN
  ENDIF
  !
  IF ( Itask<1 ) THEN
    Ind = -3
    WRITE (xern1,'(I8)') Itask
    CALL XERMSG('SLATEC','SNBFS','ITASK = '//xern1//' IS LESS THAN 1',-3,1)
    RETURN
  ENDIF
  !
  IF ( Ml<0.OR.Ml>=N ) THEN
    Ind = -5
    WRITE (xern1,'(I8)') Ml
    CALL XERMSG('SLATEC','SNBFS','ML = '//xern1//' IS OUT OF RANGE',-5,1)
    RETURN
  ENDIF
  !
  IF ( Mu<0.OR.Mu>=N ) THEN
    Ind = -6
    WRITE (xern1,'(I8)') Mu
    CALL XERMSG('SLATEC','SNBFS','MU = '//xern1//' IS OUT OF RANGE',-6,1)
    RETURN
  ENDIF
  !
  IF ( Itask==1 ) THEN
    !
    !        FACTOR MATRIX A INTO LU
    !
    CALL SNBCO(Abe,Lda,N,Ml,Mu,Iwork,rcond,Work)
    !
    !        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
    !
    IF ( rcond==0.0 ) THEN
      Ind = -4
      CALL XERMSG('SLATEC','SNBFS','SINGULAR MATRIX A - NO SOLUTION',-4,1)
      RETURN
    ENDIF
    !
    !        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
    !        AND CHECK FOR IND GREATER THAN ZERO
    !
    Ind = -LOG10(R1MACH(4)/rcond)
    IF ( Ind<=0 ) THEN
      Ind = -10
      CALL XERMSG('SLATEC','SNBFS','SOLUTION MAY HAVE NO SIGNIFICANCE',-10,&
        0)
    ENDIF
  ENDIF
  !
  !     SOLVE AFTER FACTORING
  !
  CALL SNBSL(Abe,Lda,N,Ml,Mu,Iwork,V,0)
END SUBROUTINE SNBFS
