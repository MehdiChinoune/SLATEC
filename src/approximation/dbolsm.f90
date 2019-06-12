!** DBOLSM
SUBROUTINE DBOLSM(W,Mdw,Minput,Ncols,Bl,Bu,Ind,Iopt,X,Rnorm,Mode,Rw,Ww,&
    Scl,Ibasis,Ibb)
  !>
  !  Subsidiary to DBOCLS and DBOLS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SBOLSM-S, DBOLSM-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !            **** Double Precision Version of SBOLSM ****
  !   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION ****
  !
  !          Solve E*X = F (least squares sense) with bounds on
  !            selected X values.
  !     The user must have DIMENSION statements of the form:
  !
  !       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
  !      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS)
  !       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS)
  !
  !     (Here NX=number of extra locations required for options 1,...,7;
  !     NX=0 for no options; here NI=number of extra locations possibly
  !     required for options 1-7; NI=0 for no options; NI=14 if all the
  !     options are simultaneously in use.)
  !
  !    INPUT
  !    -----
  !
  !    --------------------
  !    W(MDW,*),MINPUT,NCOLS
  !    --------------------
  !     The array W(*,*) contains the matrix [E:F] on entry. The matrix
  !     [E:F] has MINPUT rows and NCOLS+1 columns. This data is placed in
  !     the array W(*,*) with E occupying the first NCOLS columns and the
  !     right side vector F in column NCOLS+1. The row dimension, MDW, of
  !     the array W(*,*) must satisfy the inequality MDW .ge. MINPUT.
  !     Other values of MDW are errors. The values of MINPUT and NCOLS
  !     must be positive. Other values are errors.
  !
  !    ------------------
  !    BL(*),BU(*),IND(*)
  !    ------------------
  !     These arrays contain the information about the bounds that the
  !     solution values are to satisfy. The value of IND(J) tells the
  !     type of bound and BL(J) and BU(J) give the explicit values for
  !     the respective upper and lower bounds.
  !
  !    1.    For IND(J)=1, require X(J) .ge. BL(J).
  !    2.    For IND(J)=2, require X(J) .le. BU(J).
  !    3.    For IND(J)=3, require X(J) .ge. BL(J) and
  !                                X(J) .le. BU(J).
  !    4.    For IND(J)=4, no bounds on X(J) are required.
  !     The values of BL(*),BL(*) are modified by the subprogram. Values
  !     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3
  !     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an
  !     error.
  !
  !    -------
  !    IOPT(*)
  !    -------
  !     This is the array where the user can specify nonstandard options
  !     for DBOLSM. Most of the time this feature can be ignored by
  !     setting the input value IOPT(1)=99. Occasionally users may have
  !     needs that require use of the following subprogram options. For
  !     details about how to use the options see below: IOPT(*) CONTENTS.
  !
  !     Option Number   Brief Statement of Purpose
  !     ----- ------   ----- --------- -- -------
  !           1         Move the IOPT(*) processing pointer.
  !           2         Change rank determination tolerance.
  !           3         Change blow-up factor that determines the
  !                     size of variables being dropped from active
  !                     status.
  !           4         Reset the maximum number of iterations to use
  !                     in solving the problem.
  !           5         The data matrix is triangularized before the
  !                     problem is solved whenever (NCOLS/MINPUT) .lt.
  !                     FAC. Change the value of FAC.
  !           6         Redefine the weighting matrix used for
  !                     linear independence checking.
  !           7         Debug output is desired.
  !          99         No more options to change.
  !
  !    ----
  !    X(*)
  !    ----
  !     This array is used to pass data associated with options 1,2,3 and
  !     5. Ignore this input parameter if none of these options are used.
  !     Otherwise see below: IOPT(*) CONTENTS.
  !
  !    ----------------
  !    IBASIS(*),IBB(*)
  !    ----------------
  !     These arrays must be initialized by the user. The values
  !         IBASIS(J)=J, J=1,...,NCOLS
  !         IBB(J)   =1, J=1,...,NCOLS
  !     are appropriate except when using nonstandard features.
  !
  !    ------
  !    SCL(*)
  !    ------
  !     This is the array of scaling factors to use on the columns of the
  !     matrix E. These values must be defined by the user. To suppress
  !     any column scaling set SCL(J)=1.0, J=1,...,NCOLS.
  !
  !    OUTPUT
  !    ------
  !
  !    ----------
  !    X(*),RNORM
  !    ----------
  !     The array X(*) contains a solution (if MODE .ge. 0 or .eq. -22)
  !     for the constrained least squares problem. The value RNORM is the
  !     minimum residual vector length.
  !
  !    ----
  !    MODE
  !    ----
  !     The sign of mode determines whether the subprogram has completed
  !     normally, or encountered an error condition or abnormal status.
  !     A value of MODE .ge. 0 signifies that the subprogram has completed
  !     normally. The value of MODE (.ge. 0) is the number of variables
  !     in an active status: not at a bound nor at the value ZERO, for
  !     the case of free variables. A negative value of MODE will be one
  !     of the 18 cases -38,-37,...,-22, or -1. Values .lt. -1 correspond
  !     to an abnormal completion of the subprogram. To understand the
  !     abnormal completion codes see below: ERROR MESSAGES for DBOLSM
  !     An approximate solution will be returned to the user only when
  !     maximum iterations is reached, MODE=-22.
  !
  !    -----------
  !    RW(*),WW(*)
  !    -----------
  !     These are working arrays each with NCOLS entries. The array RW(*)
  !     contains the working (scaled, nonactive) solution values. The
  !     array WW(*) contains the working (scaled, active) gradient vector
  !     values.
  !
  !    ----------------
  !    IBASIS(*),IBB(*)
  !    ----------------
  !     These arrays contain information about the status of the solution
  !     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the
  !     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are
  !     the active variables. The value (IBB(J)-1) is the number of times
  !     variable J was reflected from its upper bound. (Normally the user
  !     can ignore these parameters.)
  !
  !    IOPT(*) CONTENTS
  !    ------- --------
  !     The option array allows a user to modify internal variables in
  !     the subprogram without recompiling the source code. A central
  !     goal of the initial software design was to do a good job for most
  !     people. Thus the use of options will be restricted to a select
  !     group of users. The processing of the option array proceeds as
  !     follows: a pointer, here called LP, is initially set to the value
  !     1. The value is updated as the options are processed.  At the
  !     pointer position the option number is extracted and used for
  !     locating other information that allows for options to be changed.
  !     The portion of the array IOPT(*) that is used for each option is
  !     fixed; the user and the subprogram both know how many locations
  !     are needed for each option. A great deal of error checking is
  !     done by the subprogram on the contents of the option array.
  !     Nevertheless it is still possible to give the subprogram optional
  !     input that is meaningless. For example, some of the options use
  !     the location X(NCOLS+IOFF) for passing data. The user must manage
  !     the allocation of these locations when more than one piece of
  !     option data is being passed to the subprogram.
  !
  !   1
  !   -
  !     Move the processing pointer (either forward or backward) to the
  !     location IOPT(LP+1). The processing pointer is moved to location
  !     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over
  !     locations 3,...,NCOLS+2 of IOPT(*),
  !
  !       IOPT(1)=1
  !       IOPT(2)=NCOLS+3
  !       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
  !       IOPT(NCOLS+3)=99
  !       CALL DBOLSM
  !
  !     CAUTION: Misuse of this option can yield some very hard-to-find
  !     bugs.  Use it with care.
  !
  !   2
  !   -
  !     The algorithm that solves the bounded least squares problem
  !     iteratively drops columns from the active set. This has the
  !     effect of joining a new column vector to the QR factorization of
  !     the rectangular matrix consisting of the partially triangularized
  !     nonactive columns. After triangularizing this matrix a test is
  !     made on the size of the pivot element. The column vector is
  !     rejected as dependent if the magnitude of the pivot element is
  !     .le. TOL* magnitude of the column in components strictly above
  !     the pivot element. Nominally the value of this (rank) tolerance
  !     is TOL = SQRT(R1MACH(4)). To change only the value of TOL, for
  !     example,
  !
  !       X(NCOLS+1)=TOL
  !       IOPT(1)=2
  !       IOPT(2)=1
  !       IOPT(3)=99
  !       CALL DBOLSM
  !
  !     Generally, if LP is the processing pointer for IOPT(*),
  !
  !       X(NCOLS+IOFF)=TOL
  !       IOPT(LP)=2
  !       IOPT(LP+1)=IOFF
  !        .
  !       CALL DBOLSM
  !
  !     The required length of IOPT(*) is increased by 2 if option 2 is
  !     used; The required length of X(*) is increased by 1. A value of
  !     IOFF .le. 0 is an error. A value of TOL .le. R1MACH(4) gives a
  !     warning message; it is not considered an error.
  !
  !   3
  !   -
  !     A solution component is left active (not used) if, roughly
  !     speaking, it seems too large. Mathematically the new component is
  !     left active if the magnitude is .ge.((vector norm of F)/(matrix
  !     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(R1MACH(4)).
  !     To change only the value of BLOWUP, for example,
  !
  !       X(NCOLS+2)=BLOWUP
  !       IOPT(1)=3
  !       IOPT(2)=2
  !       IOPT(3)=99
  !       CALL DBOLSM
  !
  !     Generally, if LP is the processing pointer for IOPT(*),
  !
  !       X(NCOLS+IOFF)=BLOWUP
  !       IOPT(LP)=3
  !       IOPT(LP+1)=IOFF
  !        .
  !       CALL DBOLSM
  !
  !     The required length of IOPT(*) is increased by 2 if option 3 is
  !     used; the required length of X(*) is increased by 1. A value of
  !     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error.
  !
  !   4
  !   -
  !     Normally the algorithm for solving the bounded least squares
  !     problem requires between NCOLS/3 and NCOLS drop-add steps to
  !     converge. (this remark is based on examining a small number of
  !     test cases.) The amount of arithmetic for such problems is
  !     typically about twice that required for linear least squares if
  !     there are no bounds and if plane rotations are used in the
  !     solution method. Convergence of the algorithm, while
  !     mathematically certain, can be much slower than indicated. To
  !     avoid this potential but unlikely event ITMAX drop-add steps are
  !     permitted. Nominally ITMAX=5*(MAX(MINPUT,NCOLS)). To change the
  !     value of ITMAX, for example,
  !
  !       IOPT(1)=4
  !       IOPT(2)=ITMAX
  !       IOPT(3)=99
  !       CALL DBOLSM
  !
  !     Generally, if LP is the processing pointer for IOPT(*),
  !
  !       IOPT(LP)=4
  !       IOPT(LP+1)=ITMAX
  !        .
  !       CALL DBOLSM
  !
  !     The value of ITMAX must be .gt. 0. Other values are errors. Use
  !     of this option increases the required length of IOPT(*) by 2.
  !
  !   5
  !   -
  !     For purposes of increased efficiency the MINPUT by NCOLS+1 data
  !     matrix [E:F] is triangularized as a first step whenever MINPUT
  !     satisfies FAC*MINPUT .gt. NCOLS. Nominally FAC=0.75. To change the
  !     value of FAC,
  !
  !       X(NCOLS+3)=FAC
  !       IOPT(1)=5
  !       IOPT(2)=3
  !       IOPT(3)=99
  !       CALL DBOLSM
  !
  !     Generally, if LP is the processing pointer for IOPT(*),
  !
  !       X(NCOLS+IOFF)=FAC
  !       IOPT(LP)=5
  !       IOPT(LP+1)=IOFF
  !        .
  !       CALL DBOLSM
  !
  !     The value of FAC must be nonnegative. Other values are errors.
  !     Resetting FAC=0.0 suppresses the initial triangularization step.
  !     Use of this option increases the required length of IOPT(*) by 2;
  !     The required length of of X(*) is increased by 1.
  !
  !   6
  !   -
  !     The norm used in testing the magnitudes of the pivot element
  !     compared to the mass of the column above the pivot line can be
  !     changed. The type of change that this option allows is to weight
  !     the components with an index larger than MVAL by the parameter
  !     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and
  !     WT, where LP is the processing pointer for IOPT(*),
  !
  !       X(NCOLS+IOFF)=WT
  !       IOPT(LP)=6
  !       IOPT(LP+1)=IOFF
  !       IOPT(LP+2)=MVAL
  !
  !     Use of this option increases the required length of IOPT(*) by 3.
  !     The length of X(*) is increased by 1. Values of MVAL must be
  !     nonnegative and not greater than MINPUT. Other values are errors.
  !     The value of WT must be positive. Any other value is an error. If
  !     either error condition is present a message will be printed.
  !
  !   7
  !   -
  !     Debug output, showing the detailed add-drop steps for the
  !     constrained least squares problem, is desired. This option is
  !     intended to be used to locate suspected bugs.
  !
  !   99
  !   --
  !     There are no more options to change.
  !
  !     The values for options are 1,...,7,99, and are the only ones
  !     permitted. Other values are errors. Options -99,-1,...,-7 mean
  !     that the repective options 99,1,...,7 are left at their default
  !     values. An example is the option to modify the (rank) tolerance:
  !
  !       X(NCOLS+1)=TOL
  !       IOPT(1)=-2
  !       IOPT(2)=1
  !       IOPT(3)=99
  !
  !    Error Messages for DBOLSM
  !    ----- -------- --- ---------
  !    -22    MORE THAN ITMAX = ... ITERATIONS SOLVING BOUNDED LEAST
  !           SQUARES PROBLEM.
  !
  !    -23    THE OPTION NUMBER = ... IS NOT DEFINED.
  !
  !    -24    THE OFFSET = ... BEYOND POSTION NCOLS = ... MUST BE POSITIVE
  !           FOR OPTION NUMBER 2.
  !
  !    -25    THE TOLERANCE FOR RANK DETERMINATION = ... IS LESS THAN
  !           MACHINE PRECISION = ....
  !
  !    -26    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
  !           FOR OPTION NUMBER 3.
  !
  !    -27    THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING VARIABLES
  !           MUST BE POSITIVE. NOW = ....
  !
  !    -28    THE MAXIMUM NUMBER OF ITERATIONS = ... MUST BE POSITIVE.
  !
  !    -29    THE OFFSET = ... BEYOND POSITION NCOLS = ... MUST BE POSTIVE
  !           FOR OPTION NUMBER 5.
  !
  !    -30    THE FACTOR (NCOLS/MINPUT) WHERE PRETRIANGULARIZING IS
  !           PERFORMED MUST BE NONNEGATIVE. NOW = ....
  !
  !    -31    THE NUMBER OF ROWS = ... MUST BE POSITIVE.
  !
  !    -32    THE NUMBER OF COLUMNS = ... MUST BE POSTIVE.
  !
  !    -33    THE ROW DIMENSION OF W(,) = ... MUST BE .GE. THE NUMBER OF
  !           ROWS = ....
  !
  !    -34    FOR J = ... THE CONSTRAINT INDICATOR MUST BE 1-4.
  !
  !    -35    FOR J = ... THE LOWER BOUND = ... IS .GT. THE UPPER BOUND =
  !           ....
  !
  !    -36    THE INPUT ORDER OF COLUMNS = ... IS NOT BETWEEN 1 AND NCOLS
  !           = ....
  !
  !    -37    THE BOUND POLARITY FLAG IN COMPONENT J = ... MUST BE
  !           POSITIVE. NOW = ....
  !
  !    -38    THE ROW SEPARATOR TO APPLY WEIGHTING (...) MUST LIE BETWEEN
  !           0 AND MINPUT = .... WEIGHT = ... MUST BE POSITIVE.
  !
  !***
  ! **See also:**  DBOCLS, DBOLS
  !***
  ! **Routines called:**  D1MACH, DAXPY, DCOPY, DDOT, DMOUT, DNRM2, DROT,
  !                    DROTG, DSWAP, DVOUT, IVOUT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   821220  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920422  Fixed usage of MINPUT.  (WRB)
  !   901009  Editorial changes, code now reads from top to bottom.  (RWC)

  !
  !     PURPOSE
  !     -------
  !     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED
  !     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS:
  !
  !     SOLVE E*X =  F  (LEAST SQUARES SENSE)
  !     WITH BOUNDS ON SELECTED X VALUES.
  !
  !     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
  !     EDITING AT THE CARD 'C++'.
  !     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS
  !     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/,
  !     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/,
  !     /SROT/ TO /DROT/, /SROTG/ TO /DROTG/, /R1MACH/ TO /D1MACH/,
  !     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/,
  !     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/,
  !     /REAL            / TO /DOUBLE PRECISION/.
  !++
  USE service, ONLY : XERMSG, D1MACH
  USE blas, ONLY : DROT, DROTG, DSWAP, DAXPY
  USE optimization, ONLY : DVOUT, IVOUT
  INTEGER :: Mdw, Minput, Mode, Ncols
  INTEGER :: Ibasis(Ncols), Ibb(Ncols), Ind(Ncols), Iopt(*)
  REAL(DP) :: Rnorm, W(Mdw,Ncols+1), Bl(Ncols), Bu(Ncols), X(2*Ncols+7), Rw(Ncols), &
    Ww(Ncols), Scl(Ncols)
  INTEGER :: i, igopr, ioff, ip, iprint, itemp, iter, itmax, j, jbig, jcol, &
    jdrop, jdrop1, jdrop2, jlarge, jmag, jp, lds, lgopr, lp, mrows, mval, nsetb, &
    i2(1), jbig2(1)
  REAL(DP) :: alpha, beta, bou, colabv, colblo, cl1, cl2, cl3, big, fac, &
    sc, ss, t, tolind, wt, t1, t2, wbig, wlarge, wmag, xnew, tolsze
  LOGICAL :: found, constr
  CHARACTER(8) :: xern1, xern2
  CHARACTER(16) :: xern3, xern4
  !
  REAL(DP), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0
  !
  !* FIRST EXECUTABLE STATEMENT  DBOLSM
  !
  !     Verify that the problem dimensions are defined properly.
  !
  IF ( Minput<=0 ) THEN
    WRITE (xern1,'(I8)') Minput
    CALL XERMSG('DBOLSM','THE NUMBER OF ROWS = '//xern1//&
      ' MUST BE POSITIVE.',31,1)
    Mode = -31
    RETURN
  END IF
  !
  IF ( Ncols<=0 ) THEN
    WRITE (xern1,'(I8)') Ncols
    CALL XERMSG('DBOLSM','THE NUMBER OF COLUMNS = '//xern1//&
      ' MUST BE POSITIVE.',32,1)
    Mode = -32
    RETURN
  END IF
  !
  IF ( Mdw<Minput ) THEN
    WRITE (xern1,'(I8)') Mdw
    WRITE (xern2,'(I8)') Minput
    CALL XERMSG('DBOLSM','THE ROW DIMENSION OF W(,) = '//xern1//&
      ' MUST BE .GE. THE NUMBER OF ROWS = '//xern2,33,1)
    Mode = -33
    RETURN
  END IF
  !
  !     Verify that bound information is correct.
  !
  DO j = 1, Ncols
    IF ( Ind(j)<1.OR.Ind(j)>4 ) THEN
      WRITE (xern1,'(I8)') j
      WRITE (xern2,'(I8)') Ind(j)
      CALL XERMSG('DBOLSM','FOR J = '//xern1//&
        ' THE CONSTRAINT INDICATOR MUST BE 1-4',34,1)
      Mode = -34
      RETURN
    END IF
  END DO
  !
  DO j = 1, Ncols
    IF ( Ind(j)==3 ) THEN
      IF ( Bu(j)<Bl(j) ) THEN
        WRITE (xern1,'(I8)') j
        WRITE (xern3,'(1PD15.6)') Bl(j)
        WRITE (xern4,'(1PD15.6)') Bu(j)
        CALL XERMSG('DBOLSM','FOR J = '//xern1//&
          ' THE LOWER BOUND = '//xern3//&
          ' IS .GT. THE UPPER BOUND = '//xern4,35,1)
        Mode = -35
        RETURN
      END IF
    END IF
  END DO
  !
  !     Check that permutation and polarity arrays have been set.
  !
  DO j = 1, Ncols
    IF ( Ibasis(j)<1.OR.Ibasis(j)>Ncols ) THEN
      WRITE (xern1,'(I8)') Ibasis(j)
      WRITE (xern2,'(I8)') Ncols
      CALL XERMSG('DBOLSM','THE INPUT ORDER OF COLUMNS = '//xern1//&
        ' IS NOT BETWEEN 1 AND NCOLS = '//xern2,36,1)
      Mode = -36
      RETURN
    END IF
    !
    IF ( Ibb(j)<=0 ) THEN
      WRITE (xern1,'(I8)') j
      WRITE (xern2,'(I8)') Ibb(j)
      CALL XERMSG('DBOLSM',&
        'THE BOUND POLARITY FLAG IN COMPONENT J = '//xern1//&
        ' MUST BE POSITIVE.$$NOW = '//xern2,37,1)
      Mode = -37
      RETURN
    END IF
  END DO
  !
  !     Process the option array.
  !
  fac = 0.75D0
  tolind = SQRT(D1MACH(4))
  tolsze = SQRT(D1MACH(4))
  itmax = 5*MAX(Minput,Ncols)
  wt = ONE
  mval = 0
  iprint = 0
  !
  !     Changes to some parameters can occur through the option array,
  !     IOPT(*).  Process this array looking carefully for input data
  !     errors.
  !
  lp = 0
  lds = 0
  DO
    !
    !     Test for no more options.
    !
    lp = lp + lds
    ip = Iopt(lp+1)
    jp = ABS(ip)
    IF ( ip==99 ) THEN
      !
      !     Pretriangularize rectangular arrays of certain sizes for
      !     increased efficiency.
      !
      IF ( fac*Minput>Ncols ) THEN
        DO j = 1, Ncols + 1
          DO i = Minput, j + mval + 1, -1
            CALL DROTG(W(i-1,j),W(i,j),sc,ss)
            W(i,j) = ZERO
            CALL DROT(Ncols-j+1,W(i-1,j+1),Mdw,W(i,j+1),Mdw,sc,ss)
          END DO
        END DO
        mrows = Ncols + mval + 1
      ELSE
        mrows = Minput
      END IF
      !
      !     Set the X(*) array to zero so all components are defined.
      !
      X(1:Ncols) = ZERO
      !
      !     The arrays IBASIS(*) and IBB(*) are initialized by the calling
      !     program and the column scaling is defined in the calling program.
      !     'BIG' is plus infinity on this machine.
      !
      big = D1MACH(2)
      DO j = 1, Ncols
        IF ( Ind(j)==1 ) THEN
          Bu(j) = big
        ELSEIF ( Ind(j)==2 ) THEN
          Bl(j) = -big
        ELSEIF ( Ind(j)==4 ) THEN
          Bl(j) = -big
          Bu(j) = big
        END IF
      END DO
      !
      DO j = 1, Ncols
        IF ( (Bl(j)<=ZERO.AND.ZERO<=Bu(j).AND.ABS(Bu(j))<ABS(Bl(j))).OR.&
            Bu(j)<ZERO ) THEN
          t = Bu(j)
          Bu(j) = -Bl(j)
          Bl(j) = -t
          Scl(j) = -Scl(j)
          DO i = 1, mrows
            W(i,j) = -W(i,j)
          END DO
        END IF
        !
        !         Indices in set T(=TIGHT) are denoted by negative values
        !         of IBASIS(*).
        !
        IF ( Bl(j)>=ZERO ) THEN
          Ibasis(j) = -Ibasis(j)
          t = -Bl(j)
          Bu(j) = Bu(j) + t
          CALL DAXPY(mrows,t,W(1,j),1,W(1,Ncols+1),1)
        END IF
      END DO
      !
      nsetb = 0
      iter = 0
      !
      IF ( iprint>0 ) THEN
        CALL DMOUT(mrows,Ncols+1,Mdw,W,'('' PRETRI. INPUT MATRIX'')',-4)
        CALL DVOUT(Ncols,Bl,'('' LOWER BOUNDS'')',-4)
        CALL DVOUT(Ncols,Bu,'('' UPPER BOUNDS'')',-4)
      END IF
      EXIT
    ELSEIF ( jp==99 ) THEN
      lds = 1
    ELSEIF ( jp==1 ) THEN
      !
      !         Move the IOPT(*) processing pointer.
      !
      IF ( ip>0 ) THEN
        lp = Iopt(lp+2) - 1
        lds = 0
      ELSE
        lds = 2
      END IF
    ELSEIF ( jp==2 ) THEN
      !
      !         Change tolerance for rank determination.
      !
      IF ( ip>0 ) THEN
        ioff = Iopt(lp+2)
        IF ( ioff<=0 ) THEN
          WRITE (xern1,'(I8)') ioff
          WRITE (xern2,'(I8)') Ncols
          CALL XERMSG('DBOLSM','THE OFFSET = '//xern1//&
            ' BEYOND POSITION NCOLS = '//xern2//&
            ' MUST BE POSITIVE FOR OPTION NUMBER 2.',24,1)
          Mode = -24
          RETURN
        END IF
        !
        tolind = X(Ncols+ioff)
        IF ( tolind<D1MACH(4) ) THEN
          WRITE (xern3,'(1PD15.6)') tolind
          WRITE (xern4,'(1PD15.6)') D1MACH(4)
          CALL XERMSG('DBOLSM',&
            'THE TOLERANCE FOR RANK DETERMINATION = '//xern3//&
            ' IS LESS THAN MACHINE PRECISION = '//xern4,25,0)
          Mode = -25
        END IF
      END IF
      lds = 2
    ELSEIF ( jp==3 ) THEN
      !
      !         Change blowup factor for allowing variables to become
      !         inactive.
      !
      IF ( ip>0 ) THEN
        ioff = Iopt(lp+2)
        IF ( ioff<=0 ) THEN
          WRITE (xern1,'(I8)') ioff
          WRITE (xern2,'(I8)') Ncols
          CALL XERMSG('DBOLSM','THE OFFSET = '//xern1//&
            ' BEYOND POSITION NCOLS = '//xern2//&
            ' MUST BE POSITIVE FOR OPTION NUMBER 3.',26,1)
          Mode = -26
          RETURN
        END IF
        !
        tolsze = X(Ncols+ioff)
        IF ( tolsze<=ZERO ) THEN
          WRITE (xern3,'(1PD15.6)') tolsze
          CALL XERMSG('DBOLSM','THE RECIPROCAL OF THE BLOW-UP FACTOR&
            & FOR REJECTING VARIABLES MUST BE POSITIVE.$$NOW = '//xern3,27,1)
          Mode = -27
          RETURN
        END IF
      END IF
      lds = 2
    ELSEIF ( jp==4 ) THEN
      !
      !         Change the maximum number of iterations allowed.
      !
      IF ( ip>0 ) THEN
        itmax = Iopt(lp+2)
        IF ( itmax<=0 ) THEN
          WRITE (xern1,'(I8)') itmax
          CALL XERMSG('DBOLSM',&
            'THE MAXIMUM NUMBER OF ITERATIONS = '//xern1//&
            ' MUST BE POSITIVE.',28,1)
          Mode = -28
          RETURN
        END IF
      END IF
      lds = 2
    ELSEIF ( jp==5 ) THEN
      !
      !         Change the factor for pretriangularizing the data matrix.
      !
      IF ( ip>0 ) THEN
        ioff = Iopt(lp+2)
        IF ( ioff<=0 ) THEN
          WRITE (xern1,'(I8)') ioff
          WRITE (xern2,'(I8)') Ncols
          CALL XERMSG('DBOLSM','THE OFFSET = '//xern1//&
            ' BEYOND POSITION NCOLS = '//xern2//&
            ' MUST BE POSITIVE FOR OPTION NUMBER 5.',29,1)
          Mode = -29
          RETURN
        END IF
        !
        fac = X(Ncols+ioff)
        IF ( fac<ZERO ) THEN
          WRITE (xern3,'(1PD15.6)') fac
          CALL XERMSG('DBOLSM',&
            'THE FACTOR (NCOLS/MINPUT) WHERE PRE-TRIANGULARIZING IS PERFORMED&
            & MUST BE NON-NEGATIVE.$$NOW = '//xern3,30,0)
          Mode = -30
          RETURN
        END IF
      END IF
      lds = 2
    ELSEIF ( jp==6 ) THEN
      !
      !         Change the weighting factor (from 1.0) to apply to components
      !         numbered .gt. MVAL (initially set to 1.)  This trick is needed
      !         for applications of this subprogram to the heavily weighted
      !         least squares problem that come from equality constraints.
      !
      IF ( ip>0 ) THEN
        ioff = Iopt(lp+2)
        mval = Iopt(lp+3)
        wt = X(Ncols+ioff)
      END IF
      !
      IF ( mval<0.OR.mval>Minput.OR.wt<=ZERO ) THEN
        WRITE (xern1,'(I8)') mval
        WRITE (xern2,'(I8)') Minput
        WRITE (xern3,'(1PD15.6)') wt
        CALL XERMSG('DBOLSM',&
          'THE ROW SEPARATOR TO APPLY WEIGHTING ('//xern1//&
          ') MUST LIE BETWEEN 0 AND MINPUT = '//xern2//&
          '.$$WEIGHT = '//xern3//' MUST BE POSITIVE.',38,0)
        Mode = -38
        RETURN
      END IF
      lds = 3
    ELSEIF ( jp==7 ) THEN
      !
      !         Turn on debug output.
      !
      IF ( ip>0 ) iprint = 1
      lds = 2
    ELSE
      WRITE (xern1,'(I8)') ip
      CALL XERMSG('DBOLSM','THE OPTION NUMBER = '//xern1//&
        ' IS NOT DEFINED.',23,1)
      Mode = -23
      RETURN
    END IF
  END DO
  !
  100  iter = iter + 1
  IF ( iter>itmax ) THEN
    WRITE (xern1,'(I8)') itmax
    CALL XERMSG('DBOLSM','MORE THAN ITMAX = '//xern1//&
      ' ITERATIONS SOLVING BOUNDED LEAST SQUARES PROBLEM.',22,1)
    Mode = -22
    !
    !        Rescale and translate variables.
    !
    igopr = 1
    GOTO 700
  END IF
  !
  !     Find a variable to become non-active.
  !                                                 T
  !     Compute (negative) of gradient vector, W = E *(F-E*X).
  !
  Ww(1:Ncols) = ZERO
  DO j = nsetb + 1, Ncols
    jcol = ABS(Ibasis(j))
    Ww(j) = DOT_PRODUCT(W(INEXT(nsetb):INEXT(nsetb)+mrows-nsetb-1,j), &
      W(INEXT(nsetb):INEXT(nsetb)+mrows-nsetb-1,Ncols+1)) *ABS(Scl(jcol))
  END DO
  !
  IF ( iprint>0 ) THEN
    CALL DVOUT(Ncols,Ww,'('' GRADIENT VALUES'')',-4)
    CALL IVOUT(Ncols,Ibasis,'('' INTERNAL VARIABLE ORDER'')',-4)
    CALL IVOUT(Ncols,Ibb,'('' BOUND POLARITY'')',-4)
  END IF
  DO
    !
    !     If active set = number of total rows, quit.
    !
    IF ( nsetb==mrows ) THEN
      found = .FALSE.
      GOTO 600
    END IF
    !
    !     Choose an extremal component of gradient vector for a candidate
    !     to become non-active.
    !
    wlarge = -big
    wmag = -big
    DO j = nsetb + 1, Ncols
      t = Ww(j)
      IF ( t/=big ) THEN
        itemp = Ibasis(j)
        jcol = ABS(itemp)
        t1 = NORM2(W(INEXT(nsetb):INEXT(nsetb)+mval-nsetb-1,j))
        IF ( itemp<0 ) THEN
          IF ( MOD(Ibb(jcol),2)==0 ) t = -t
          IF ( t>=ZERO ) THEN
            IF ( mval>nsetb ) t = t1
            IF ( t>wlarge ) THEN
              wlarge = t
              jlarge = j
            END IF
          END IF
        ELSE
          IF ( mval>nsetb ) t = t1
          IF ( ABS(t)>wmag ) THEN
            wmag = ABS(t)
            jmag = j
          END IF
        END IF
      END IF
    END DO
    !
    !     Choose magnitude of largest component of gradient for candidate.
    !
    jbig = 0
    wbig = ZERO
    IF ( wlarge>ZERO ) THEN
      jbig = jlarge
      wbig = wlarge
    END IF
    !
    IF ( wmag>=wbig ) THEN
      jbig = jmag
      wbig = wmag
    END IF
    !
    IF ( jbig==0 ) THEN
      found = .FALSE.
      IF ( iprint>0 ) CALL IVOUT(0,i2,'('' FOUND NO VARIABLE TO ENTER'')',-4)
      GOTO 600
    END IF
    !
    !     See if the incoming column is sufficiently independent.  This
    !     test is made before an elimination is performed.
    !
    IF ( iprint>0 ) CALL IVOUT(1,jbig2,'('' TRY TO BRING IN THIS COL.'')',-4)
    !
    IF ( mval<=nsetb ) THEN
      cl1 = NORM2(W(1:mval,jbig))
      cl2 = ABS(wt)*NORM2(W(INEXT(mval):INEXT(mval)+nsetb-mval-1,jbig))
      cl3 = ABS(wt)*NORM2(W(INEXT(nsetb):INEXT(nsetb)+mrows-nsetb-1,jbig))
      CALL DROTG(cl1,cl2,sc,ss)
      colabv = ABS(cl1)
      colblo = cl3
    ELSE
      cl1 = NORM2(W(1:nsetb,jbig))
      cl2 = NORM2(W(INEXT(nsetb):INEXT(nsetb)+mval-nsetb-1,jbig))
      cl3 = ABS(wt)*NORM2(W(INEXT(mval):INEXT(mval)+mrows-mval-1,jbig))
      colabv = cl1
      CALL DROTG(cl2,cl3,sc,ss)
      colblo = ABS(cl2)
    END IF
    !
    IF ( colblo<=tolind*colabv ) THEN
      Ww(jbig) = big
      IF ( iprint>0 ) CALL IVOUT(0,i2,'('' VARIABLE IS DEPENDENT, NOT USED.'')',-4)
      CYCLE
    END IF
    !
    !     Swap matrix columns NSETB+1 and JBIG, plus pointer information,
    !     and gradient values.
    !
    nsetb = nsetb + 1
    IF ( nsetb/=jbig ) THEN
      CALL DSWAP(mrows,W(1,nsetb),1,W(1,jbig),1)
      CALL DSWAP(1,Ww(nsetb),1,Ww(jbig),1)
      itemp = Ibasis(nsetb)
      Ibasis(nsetb) = Ibasis(jbig)
      Ibasis(jbig) = itemp
    END IF
    !
    !     Eliminate entries below the pivot line in column NSETB.
    !
    IF ( mrows>nsetb ) THEN
      DO i = mrows, nsetb + 1, -1
        IF ( i/=mval+1 ) THEN
          CALL DROTG(W(i-1,nsetb),W(i,nsetb),sc,ss)
          W(i,nsetb) = ZERO
          CALL DROT(Ncols-nsetb+1,W(i-1,nsetb+1),Mdw,W(i,nsetb+1),Mdw,sc,ss)
        END IF
      END DO
      !
      IF ( mval>=nsetb.AND.mval<mrows ) THEN
        CALL DROTG(W(nsetb,nsetb),W(mval+1,nsetb),sc,ss)
        W(mval+1,nsetb) = ZERO
        CALL DROT(Ncols-nsetb+1,W(nsetb,nsetb+1),Mdw,W(mval+1,nsetb+1),Mdw,&
          sc,ss)
      END IF
    END IF
    !
    IF ( W(nsetb,nsetb)==ZERO ) THEN
      Ww(nsetb) = big
      nsetb = nsetb - 1
      IF ( iprint>0 ) CALL IVOUT(0,i2,'('' PIVOT IS ZERO, NOT USED.'')',-4)
      CYCLE
    END IF
    !
    !     Check that new variable is moving in the right direction.
    !
    itemp = Ibasis(nsetb)
    jcol = ABS(itemp)
    xnew = (W(nsetb,Ncols+1)/W(nsetb,nsetb))/ABS(Scl(jcol))
    IF ( itemp<0 ) THEN
      !
      !         IF(WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) exit(quit)
      !         IF(WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO) exit(quit)
      !
      IF ( (Ww(nsetb)>=ZERO.AND.xnew<=ZERO).OR.&
          (Ww(nsetb)<=ZERO.AND.xnew>=ZERO) ) THEN
        !
        Ww(nsetb) = big
        nsetb = nsetb - 1
        IF ( iprint>0 ) CALL IVOUT(0,i2,&
          '('' VARIABLE HAS BAD DIRECTION, NOT USED.'')',-4)
        CYCLE
      END IF
    END IF
    EXIT
  END DO
  found = .TRUE.
  GOTO 600
  !
  !     Solve the triangular system.
  !
  200  Rw(1:nsetb) = W(1:nsetb,Ncols+1)
  DO j = nsetb, 1, -1
    Rw(j) = Rw(j)/W(j,j)
    jcol = ABS(Ibasis(j))
    t = Rw(j)
    IF ( MOD(Ibb(jcol),2)==0 ) Rw(j) = -Rw(j)
    CALL DAXPY(j-1,-t,W(1,j),1,Rw,1)
    Rw(j) = Rw(j)/ABS(Scl(jcol))
  END DO
  !
  IF ( iprint>0 ) THEN
    CALL DVOUT(nsetb,Rw,'('' SOLN. VALUES'')',-4)
    CALL IVOUT(nsetb,Ibasis,'('' COLS. USED'')',-4)
  END IF
  !
  IF ( lgopr==2 ) THEN
    X(1:nsetb) = Rw(1:nsetb)
    DO j = 1, nsetb
      itemp = Ibasis(j)
      jcol = ABS(itemp)
      IF ( itemp<0 ) THEN
        bou = ZERO
      ELSE
        bou = Bl(jcol)
      END IF
      !
      IF ( (-bou)/=big ) bou = bou/ABS(Scl(jcol))
      IF ( X(j)<=bou ) THEN
        jdrop1 = j
        EXIT
      END IF
      !
      bou = Bu(jcol)
      IF ( bou/=big ) bou = bou/ABS(Scl(jcol))
      IF ( X(j)>=bou ) THEN
        jdrop2 = j
        EXIT
      END IF
    END DO
    GOTO 300
  END IF
  !
  !     See if the unconstrained solution (obtained by solving the
  !     triangular system) satisfies the problem bounds.
  !
  alpha = TWO
  beta = TWO
  X(nsetb) = ZERO
  DO j = 1, nsetb
    itemp = Ibasis(j)
    jcol = ABS(itemp)
    t1 = TWO
    t2 = TWO
    IF ( itemp<0 ) THEN
      bou = ZERO
    ELSE
      bou = Bl(jcol)
    END IF
    IF ( (-bou)/=big ) bou = bou/ABS(Scl(jcol))
    IF ( Rw(j)<=bou ) t1 = (X(j)-bou)/(X(j)-Rw(j))
    bou = Bu(jcol)
    IF ( bou/=big ) bou = bou/ABS(Scl(jcol))
    IF ( Rw(j)>=bou ) t2 = (bou-X(j))/(Rw(j)-X(j))
    !
    !     If not, then compute a step length so that the variables remain
    !     feasible.
    !
    IF ( t1<alpha ) THEN
      alpha = t1
      jdrop1 = j
    END IF
    !
    IF ( t2<beta ) THEN
      beta = t2
      jdrop2 = j
    END IF
  END DO
  !
  constr = alpha<TWO .OR. beta<TWO
  IF ( .NOT.constr ) THEN
    !
    !         Accept the candidate because it satisfies the stated bounds
    !         on the variables.
    !
    X(1:nsetb) = Rw(1:nsetb)
    GOTO 100
  END IF
  !
  !     Take a step that is as large as possible with all variables
  !     remaining feasible.
  !
  DO j = 1, nsetb
    X(j) = X(j) + MIN(alpha,beta)*(Rw(j)-X(j))
  END DO
  !
  IF ( alpha<=beta ) THEN
    jdrop2 = 0
  ELSE
    jdrop1 = 0
  END IF
  !
  300 CONTINUE
  IF ( jdrop1+jdrop2<=0.OR.nsetb<=0 ) GOTO 100
  jdrop = jdrop1 + jdrop2
  itemp = Ibasis(jdrop)
  jcol = ABS(itemp)
  IF ( jdrop2>0 ) THEN
    !
    !         Variable is at an upper bound.  Subtract multiple of this
    !         column from right hand side.
    !
    t = Bu(jcol)
    IF ( itemp>0 ) THEN
      Bu(jcol) = t - Bl(jcol)
      Bl(jcol) = -t
      itemp = -itemp
      Scl(jcol) = -Scl(jcol)
      DO i = 1, jdrop
        W(i,jdrop) = -W(i,jdrop)
      END DO
    ELSE
      Ibb(jcol) = Ibb(jcol) + 1
      IF ( MOD(Ibb(jcol),2)==0 ) t = -t
    END IF
    !
    !     Variable is at a lower bound.
    !
  ELSEIF ( itemp<ZERO ) THEN
    t = ZERO
  ELSE
    t = -Bl(jcol)
    Bu(jcol) = Bu(jcol) + t
    itemp = -itemp
  END IF
  !
  CALL DAXPY(jdrop,t,W(1,jdrop),1,W(1,Ncols+1),1)
  !
  !     Move certain columns left to achieve upper Hessenberg form.
  !
  Rw(1:jdrop) = W(1:jdrop,jdrop)
  DO j = jdrop + 1, nsetb
    Ibasis(j-1) = Ibasis(j)
    X(j-1) = X(j)
    W(1:j,j-1) = W(1:j,j)
  END DO
  !
  Ibasis(nsetb) = itemp
  W(1,nsetb) = ZERO
  W(jdrop+1:mrows,nsetb) = ZERO
  W(1:jdrop,nsetb) = Rw(1:jdrop)
  !
  !     Transform the matrix from upper Hessenberg form to upper
  !     triangular form.
  !
  nsetb = nsetb - 1
  DO i = jdrop, nsetb
    !
    !         Look for small pivots and avoid mixing weighted and
    !         nonweighted rows.
    !
    IF ( i==mval ) THEN
      t = ZERO
      DO j = i, nsetb
        jcol = ABS(Ibasis(j))
        t1 = ABS(W(i,j)*Scl(jcol))
        IF ( t1>t ) THEN
          jbig = j
          t = t1
        END IF
      END DO
      GOTO 400
    END IF
    CALL DROTG(W(i,i),W(i+1,i),sc,ss)
    W(i+1,i) = ZERO
    CALL DROT(Ncols-i+1,W(i,i+1),Mdw,W(i+1,i+1),Mdw,sc,ss)
  END DO
  GOTO 500
  !
  !     The triangularization is completed by giving up the Hessenberg
  !     form and triangularizing a rectangular matrix.
  !
  400  CALL DSWAP(mrows,W(1,i),1,W(1,jbig),1)
  CALL DSWAP(1,Ww(i),1,Ww(jbig),1)
  CALL DSWAP(1,X(i),1,X(jbig),1)
  itemp = Ibasis(i)
  Ibasis(i) = Ibasis(jbig)
  Ibasis(jbig) = itemp
  jbig = i
  DO j = jbig, nsetb
    DO i = j + 1, mrows
      CALL DROTG(W(j,j),W(i,j),sc,ss)
      W(i,j) = ZERO
      CALL DROT(Ncols-j+1,W(j,j+1),Mdw,W(i,j+1),Mdw,sc,ss)
    END DO
  END DO
  !
  !     See if the remaining coefficients are feasible.  They should be
  !     because of the way MIN(ALPHA,BETA) was chosen.  Any that are not
  !     feasible will be set to their bounds and appropriately translated.
  !
  500  jdrop1 = 0
  jdrop2 = 0
  lgopr = 2
  GOTO 200
  !
  !     Find a variable to become non-active.
  !
  600 CONTINUE
  IF ( found ) THEN
    lgopr = 1
    GOTO 200
  END IF
  !
  !     Rescale and translate variables.
  !
  igopr = 2
  700  Rw(1:nsetb) = X(1:nsetb)
  X(1:Ncols) = ZERO
  DO j = 1, nsetb
    jcol = ABS(Ibasis(j))
    X(jcol) = Rw(j)*ABS(Scl(jcol))
  END DO
  !
  DO j = 1, Ncols
    IF ( MOD(Ibb(j),2)==0 ) X(j) = Bu(j) - X(j)
  END DO
  !
  DO j = 1, Ncols
    jcol = Ibasis(j)
    IF ( jcol<0 ) X(-jcol) = Bl(-jcol) + X(-jcol)
  END DO
  !
  DO j = 1, Ncols
    IF ( Scl(j)<ZERO ) X(j) = -X(j)
  END DO
  !
  i = MAX(nsetb,mval)
  Rnorm = NORM2(W(INEXT(i):INEXT(i)+mrows-i-1,Ncols+1))
  !
  IF ( igopr==2 ) Mode = nsetb
CONTAINS
  INTEGER FUNCTION INEXT(idum)
    INTEGER, INTENT(IN) :: idum
    INEXT = MIN(idum+1,mrows)
  END FUNCTION INEXT
END SUBROUTINE DBOLSM
