!** EFC
SUBROUTINE EFC(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Mdein,Mdeout,Coeff,Lw,W)
  !> Fit a piecewise polynomial curve to discrete data.
  !            The piecewise polynomials are represented as B-splines.
  !            The fitting is done in a weighted least squares sense.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A1A1, K1A2A, L8A3
  !***
  ! **Type:**      SINGLE PRECISION (EFC-S, DEFC-D)
  !***
  ! **Keywords:**  B-SPLINE, CURVE FITTING, WEIGHTED LEAST SQUARES
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !      This subprogram fits a piecewise polynomial curve
  !      to discrete data.  The piecewise polynomials are
  !      represented as B-splines.
  !      The fitting is done in a weighted least squares sense.
  !
  !      The data can be processed in groups of modest size.
  !      The size of the group is chosen by the user.  This feature
  !      may be necessary for purposes of using constrained curve fitting
  !      with subprogram FC( ) on a very large data set.
  !
  !      For a description of the B-splines and usage instructions to
  !      evaluate them, see
  !
  !      C. W. de Boor, Package for Calculating with B-Splines.
  !                     SIAM J. Numer. Anal., p. 441, (June, 1977).
  !
  !      For further discussion of (constrained) curve fitting using
  !      B-splines, see
  !
  !      R. J. Hanson, Constrained Least Squares Curve Fitting
  !                   to Discrete Data Using B-Splines, a User's
  !                   Guide. Sandia Labs. Tech. Rept. SAND-78-1291,
  !                   December, (1978).
  !
  !  Input..
  !      NDATA,XDATA(*),
  !      YDATA(*),
  !      SDDATA(*)
  !                         The NDATA discrete (X,Y) pairs and the Y value
  !                         standard deviation or uncertainty, SD, are in
  !                         the respective arrays XDATA(*), YDATA(*), and
  !                         SDDATA(*).  No sorting of XDATA(*) is
  !                         required.  Any non-negative value of NDATA is
  !                         allowed.  A negative value of NDATA is an
  !                         error.  A zero value for any entry of
  !                         SDDATA(*) will weight that data point as 1.
  !                         Otherwise the weight of that data point is
  !                         the reciprocal of this entry.
  !
  !      NORD,NBKPT,
  !      BKPT(*)
  !                         The NBKPT knots of the B-spline of order NORD
  !                         are in the array BKPT(*).  Normally the
  !                         problem data interval will be included between
  !                         the limits BKPT(NORD) and BKPT(NBKPT-NORD+1).
  !                         The additional end knots BKPT(I),I=1,...,
  !                         NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
  !                         required to compute the functions used to fit
  !                         the data.  No sorting of BKPT(*) is required.
  !                         Internal to  EFC( ) the extreme end knots may
  !                         be reduced and increased respectively to
  !                         accommodate any data values that are exterior
  !                         to the given knot values.  The contents of
  !                         BKPT(*) is not changed.
  !
  !                         NORD must be in the range 1 <= NORD <= 20.
  !                         The value of NBKPT must satisfy the condition
  !                         NBKPT >= 2*NORD.
  !                         Other values are considered errors.
  !
  !                         (The order of the spline is one more than the
  !                         degree of the piecewise polynomial defined on
  !                         each interval.  This is consistent with the
  !                         B-spline package convention.  For example,
  !                         NORD=4 when we are using piecewise cubics.)
  !
  !        MDEIN
  !                         An integer flag, with one of two possible
  !                         values (1 or 2), that directs the subprogram
  !                         action with regard to new data points provided
  !                         by the user.
  !
  !                         =1  The first time that EFC( ) has been
  !                         entered.  There are NDATA points to process.
  !
  !                         =2  This is another entry to EFC( ).  The sub-
  !                         program EFC( ) has been entered with MDEIN=1
  !                         exactly once before for this problem.  There
  !                         are NDATA new additional points to merge and
  !                         process with any previous points.
  !                         (When using EFC( ) with MDEIN=2 it is import-
  !                         ant that the set of knots remain fixed at the
  !                         same values for all entries to EFC( ).)
  !       LW
  !                         The amount of working storage actually
  !                         allocated for the working array W(*).
  !                         This quantity is compared with the
  !                         actual amount of storage needed in EFC( ).
  !                         Insufficient storage allocated for W(*) is
  !                         an error.  This feature was included in EFC( )
  !                         because misreading the storage formula
  !                         for W(*) might very well lead to subtle
  !                         and hard-to-find programming bugs.
  !
  !                         The length of the array W(*) must satisfy
  !
  !                         LW >= (NBKPT-NORD+3)*(NORD+1)+
  !                                 (NBKPT+1)*(NORD+1)+
  !                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
  !
  !  Output..
  !      MDEOUT
  !                         An output flag that indicates the status
  !                         of the curve fit.
  !
  !                         =-1  A usage error of EFC( ) occurred.  The
  !                         offending condition is noted with the SLATEC
  !                         library error processor, XERMSG( ).  In case
  !                         the working array W(*) is not long enough, the
  !                         minimal acceptable length is printed.
  !
  !                         =1  The B-spline coefficients for the fitted
  !                         curve have been returned in array COEFF(*).
  !
  !                         =2  Not enough data has been processed to
  !                         determine the B-spline coefficients.
  !                         The user has one of two options.  Continue
  !                         to process more data until a unique set
  !                         of coefficients is obtained, or use the
  !                         subprogram FC( ) to obtain a specific
  !                         set of coefficients.  The user should read
  !                         the usage instructions for FC( ) for further
  !                         details if this second option is chosen.
  !      COEFF(*)
  !                         If the output value of MDEOUT=1, this array
  !                         contains the unknowns obtained from the least
  !                         squares fitting process.  These N=NBKPT-NORD
  !                         parameters are the B-spline coefficients.
  !                         For MDEOUT=2, not enough data was processed to
  !                         uniquely determine the B-spline coefficients.
  !                         In this case, and also when MDEOUT=-1, all
  !                         values of COEFF(*) are set to zero.
  !
  !                         If the user is not satisfied with the fitted
  !                         curve returned by EFC( ), the constrained
  !                         least squares curve fitting subprogram FC( )
  !                         may be required.  The work done within EFC( )
  !                         to accumulate the data can be utilized by
  !                         the user, if so desired.  This involves
  !                         saving the first (NBKPT-NORD+3)*(NORD+1)
  !                         entries of W(*) and providing this data
  !                         to FC( ) with the "old problem" designation.
  !                         The user should read the usage instructions
  !                         for subprogram FC( ) for further details.
  !
  !  Working Array..
  !      W(*)
  !                         This array is typed REAL.
  !                         Its length is  specified as an input parameter
  !                         in LW as noted above.  The contents of W(*)
  !                         must not be modified by the user between calls
  !                         to EFC( ) with values of MDEIN=1,2,2,... .
  !                         The first (NBKPT-NORD+3)*(NORD+1) entries of
  !                         W(*) are acceptable as direct input to FC( )
  !                         for an "old problem" only when MDEOUT=1 or 2.
  !
  !  Evaluating the
  !  Fitted Curve..
  !                         To evaluate derivative number IDER at XVAL,
  !                         use the function subprogram BVALU( ).
  !
  !                         F = BVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
  !                                      XVAL,INBV,WORKB)
  !
  !                         The output of this subprogram will not be
  !                         defined unless an output value of MDEOUT=1
  !                         was obtained from EFC( ), XVAL is in the data
  !                         interval, and IDER is nonnegative and <
  !                         NORD.
  !
  !                         The first time BVALU( ) is called, INBV=1
  !                         must be specified.  This value of INBV is the
  !                         overwritten by BVALU( ).  The array WORKB(*)
  !                         must be of length at least 3*NORD, and must
  !                         not be the same as the W(*) array used in the
  !                         call to EFC( ).
  !
  !                         BVALU( ) expects the breakpoint array BKPT(*)
  !                         to be sorted.
  !
  !***
  ! **References:**  R. J. Hanson, Constrained least squares curve fitting
  !                 to discrete data using B-splines, a users guide,
  !                 Report SAND78-1291, Sandia Laboratories, December
  !                 1978.
  !***
  ! **Routines called:**  EFCMN

  !* REVISION HISTORY  (YYMMDD)
  !   800801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Change Prologue comments to refer to XERMSG.  (RWC)
  !   900607  Editorial changes to Prologue to make Prologues for EFC,
  !           DEFC, FC, and DFC look as much the same as possible.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  !      SUBROUTINE           FUNCTION/REMARKS
  !
  !      BSPLVN( )          COMPUTE FUNCTION VALUES OF B-SPLINES.  FROM
  !                         THE B-SPLINE PACKAGE OF DE BOOR NOTED ABOVE.
  !
  !      BNDACC( ),         BANDED LEAST SQUARES MATRIX PROCESSORS.
  !      BNDSOL( )          FROM LAWSON-HANSON, SOLVING LEAST
  !                         SQUARES PROBLEMS.
  !
  !      SSORT( )           DATA SORTING SUBROUTINE, FROM THE
  !                         SANDIA MATH. LIBRARY, SAND77-1441.
  !
  !      XERMSG( )          ERROR HANDLING ROUTINE
  !                         FOR THE SLATEC MATH. LIBRARY.
  !                         SEE SAND78-1189, BY R. E. JONES.
  !
  !      SCOPY( ),SSCAL( )  SUBROUTINES FROM THE BLAS PACKAGE.
  !
  !                         WRITTEN BY R. HANSON, SANDIA NATL. LABS.,
  !                         ALB., N. M., AUGUST-SEPTEMBER, 1980.
  !
  INTEGER :: Lw, Mdein, Mdeout, Nbkpt, Ndata, Nord
  REAL(SP) :: Bkpt(Nbkpt), Coeff(Nbkpt-Nord), Sddata(Ndata), W(Lw), Xdata(Ndata), &
    Ydata(Ndata)
  !
  INTEGER :: lbf, lbkpt, lg, lptemp, lww, lxtemp, mdg, mdw
  !
  !* FIRST EXECUTABLE STATEMENT  EFC
  !     LWW=1               USAGE IN EFCMN( ) OF W(*)..
  !     LWW,...,LG-1        W(*,*)
  !
  !     LG,...,LXTEMP-1     G(*,*)
  !
  !     LXTEMP,...,LPTEMP-1 XTEMP(*)
  !
  !     LPTEMP,...,LBKPT-1  PTEMP(*)
  !
  !     LBKPT,...,LBF       BKPT(*) (LOCAL TO EFCMN( ))
  !
  !     LBF,...,LBF+NORD**2 BF(*,*)
  !
  mdg = Nbkpt + 1
  mdw = Nbkpt - Nord + 3
  lww = 1
  lg = lww + mdw*(Nord+1)
  lxtemp = lg + mdg*(Nord+1)
  lptemp = lxtemp + MAX(Ndata,Nbkpt)
  lbkpt = lptemp + MAX(Ndata,Nbkpt)
  lbf = lbkpt + Nbkpt
  CALL EFCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkpt,Mdein,Mdeout,Coeff,&
    W(lbf),W(lxtemp),W(lptemp),W(lbkpt),W(lg),mdg,W(lww),mdw,Lw)
END SUBROUTINE EFC
