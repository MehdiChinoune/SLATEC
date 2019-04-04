!** CGEIR
SUBROUTINE CGEIR(A,Lda,N,V,Itask,Ind,Work,Iwork)
  IMPLICIT NONE
  !>
  !***
  !  Solve a general system of linear equations.  Iterative
  !            refinement is used to obtain an error estimate.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2C1
  !***
  ! **Type:**      COMPLEX (SGEIR-S, CGEIR-C)
  !***
  ! **Keywords:**  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
  !             GENERAL SYSTEM OF LINEAR EQUATIONS
  !***
  ! **Author:**  Voorhees, E. A., (LANL)
  !***
  ! **Description:**
  !
  !    Subroutine CGEIR solves a general NxN system of complex
  !    linear equations using LINPACK subroutines CGEFA and CGESL.
  !    One pass of iterative refinement is used only to obtain an
  !    estimate of the accuracy.  That is, if A is an NxN complex
  !    matrix and if X and B are complex N-vectors, then CGEIR solves
  !    the equation
  !
  !                          A*X=B.
  !
  !    The matrix A is first factored into upper and lower tri-
  !    angular matrices U and L using partial pivoting.  These
  !    factors and the pivoting information are used to calculate
  !    the solution, X.  Then the residual vector is found and
  !    used to calculate an estimate of the relative error, IND.
  !    IND estimates the accuracy of the solution only when the
  !    input matrix and the right hand side are represented
  !    exactly in the computer and does not take into
  !    account any errors in the input data.
  !
  !    If the equation A*X=B is to be solved for more than one vector
  !    B, the factoring of A does not need to be performed again and
  !    the option to only solve (ITASK .GT. 1) will be faster for
  !    the succeeding solutions.  In this case, the contents of A,
  !    LDA, N, WORK, and IWORK must not have been altered by the
  !    user following factorization (ITASK=1).  IND will not be
  !    changed by CGEIR in this case.
  !
  !  Argument Description ***
  !
  !    A      COMPLEX(LDA,N)
  !             the doubly subscripted array with dimension (LDA,N)
  !             which contains the coefficient matrix.  A is not
  !             altered by the routine.
  !    LDA    INTEGER
  !             the leading dimension of the array A.  LDA must be great-
  !             er than or equal to N.  (Terminal error message IND=-1)
  !    N      INTEGER
  !             the order of the matrix A.  The first N elements of
  !             the array A are the elements of the first column of
  !             matrix A.  N must be greater than or equal to 1.
  !             (Terminal error message IND=-2)
  !    V      COMPLEX(N)
  !             on entry, the singly subscripted array(vector) of di-
  !               mension N which contains the right hand side B of a
  !               system of simultaneous linear equations A*X=B.
  !             on return, V contains the solution vector, X .
  !    ITASK  INTEGER
  !             if ITASK=1, the matrix A is factored and then the
  !               linear equation is solved.
  !             if ITASK .GT. 1, the equation is solved using the existing
  !               factored matrix A (stored in work).
  !             if ITASK .LT. 1, then terminal error message IND=-3 is
  !               printed.
  !    IND    INTEGER
  !             GT.0  IND is a rough estimate of the number of digits
  !                     of accuracy in the solution, X.  IND=75 means
  !                     that the solution vector X is zero.
  !             LT.0  see error message corresponding to IND below.
  !    WORK   COMPLEX(N*(N+1))
  !             a singly subscripted array of dimension at least N*(N+1).
  !    IWORK  INTEGER(N)
  !             a singly subscripted array of dimension at least N.
  !
  !  Error Messages Printed ***
  !
  !    IND=-1  terminal   N is greater than LDA.
  !    IND=-2  terminal   N is less than one.
  !    IND=-3  terminal   ITASK is less than one.
  !    IND=-4  terminal   The matrix A is computationally singular.
  !                         A solution has not been computed.
  !    IND=-10 warning    The solution has no apparent significance.
  !                         The solution may be inaccurate or the matrix
  !                         A may be poorly scaled.
  !
  !               NOTE-  The above terminal(*fatal*) error messages are
  !                      designed to be handled by XERMSG in which
  !                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
  !                      for warning error messages from XERMSG.  Unless
  !                      the user provides otherwise, an error message
  !                      will be printed followed by an abort.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CCOPY, CDCDOT, CGEFA, CGESL, R1MACH, SCASUM, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800502  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
  !           IF-THEN-ELSE.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER Lda, N, Itask, Ind, Iwork(*), info, j
  COMPLEX A(Lda,*), V(*), Work(N,*), CDCDOT
  REAL SCASUM, xnorm, dnorm, R1MACH
  CHARACTER(8) :: xern1, xern2
  !* FIRST EXECUTABLE STATEMENT  CGEIR
  IF ( Lda<N ) THEN
    Ind = -1
    WRITE (xern1,'(I8)') Lda
    WRITE (xern2,'(I8)') N
    CALL XERMSG('SLATEC','CGEIR','LDA = '//xern1//' IS LESS THAN N = '//&
      xern2,-1,1)
    RETURN
  END IF
  !
  IF ( N<=0 ) THEN
    Ind = -2
    WRITE (xern1,'(I8)') N
    CALL XERMSG('SLATEC','CGEIR','N = '//xern1//' IS LESS THAN 1',-2,1)
    RETURN
  END IF
  !
  IF ( Itask<1 ) THEN
    Ind = -3
    WRITE (xern1,'(I8)') Itask
    CALL XERMSG('SLATEC','CGEIR','ITASK = '//xern1//' IS LESS THAN 1',-3,1)
    RETURN
  END IF
  !
  IF ( Itask==1 ) THEN
    !        MOVE MATRIX A TO WORK
    DO j = 1, N
      CALL CCOPY(N,A(1,j),1,Work(1,j),1)
    END DO
    !
    !        FACTOR MATRIX A INTO LU
    !
    CALL CGEFA(Work,N,N,Iwork,info)
    !
    !        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
    !
    IF ( info/=0 ) THEN
      Ind = -4
      CALL XERMSG('SLATEC','CGEIR','SINGULAR MATRIX A - NO SOLUTION',-4,1)
      RETURN
    END IF
  END IF
  !
  !     SOLVE WHEN FACTORING COMPLETE
  !     MOVE VECTOR B TO WORK
  !
  CALL CCOPY(N,V(1),1,Work(1,N+1),1)
  CALL CGESL(Work,N,N,Iwork,V,0)
  !
  !     FORM NORM OF X0
  !
  xnorm = SCASUM(N,V(1),1)
  IF ( xnorm==0.0 ) THEN
    Ind = 75
    RETURN
  END IF
  !
  !     COMPUTE  RESIDUAL
  !
  DO j = 1, N
    Work(j,N+1) = CDCDOT(N,-Work(j,N+1),A(j,1),Lda,V,1)
  END DO
  !
  !     SOLVE A*DELTA=R
  !
  CALL CGESL(Work,N,N,Iwork,Work(1,N+1),0)
  !
  !     FORM NORM OF DELTA
  !
  dnorm = SCASUM(N,Work(1,N+1),1)
  !
  !     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
  !     AND CHECK FOR IND GREATER THAN ZERO
  !
  Ind = INT( -LOG10(MAX(R1MACH(4),dnorm/xnorm)) )
  IF ( Ind<=0 ) THEN
    Ind = -10
    CALL XERMSG('SLATEC','CGEIR','SOLUTION MAY HAVE NO SIGNIFICANCE',-10,0)
  END IF
END SUBROUTINE CGEIR
