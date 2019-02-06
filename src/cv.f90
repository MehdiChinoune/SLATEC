!*==CV.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CV
      REAL FUNCTION CV(Xval,Ndata,Nconst,Nord,Nbkpt,Bkpt,W)
      IMPLICIT NONE
!*--CV5
!*** Start of declarations inserted by SPAG
      REAL Bkpt , SDOT , v , W , Xval , zero
      INTEGER i , ileft , ip , is , last , mdg , mdw , n , Nbkpt , Nconst , 
     &        Ndata , Nord
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CV
!***PURPOSE  Evaluate the variance function of the curve obtained
!            by the constrained B-spline fitting subprogram FC.
!***LIBRARY   SLATEC
!***CATEGORY  L7A3
!***TYPE      SINGLE PRECISION (CV-S, DCV-D)
!***KEYWORDS  ANALYSIS OF COVARIANCE, B-SPLINE,
!             CONSTRAINED LEAST SQUARES, CURVE FITTING
!***AUTHOR  Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     CV( ) is a companion function subprogram for FC( ).  The
!     documentation for FC( ) has complete usage instructions.
!
!     CV( ) is used to evaluate the variance function of the curve
!     obtained by the constrained B-spline fitting subprogram, FC( ).
!     The variance function defines the square of the probable error
!     of the fitted curve at any point, XVAL.  One can use the square
!     root of this variance function to determine a probable error band
!     around the fitted curve.
!
!     CV( ) is used after a call to FC( ).  MODE, an input variable to
!     FC( ), is used to indicate if the variance function is desired.
!     In order to use CV( ), MODE must equal 2 or 4 on input to FC( ).
!     MODE is also used as an output flag from FC( ).  Check to make
!     sure that MODE = 0 after calling FC( ), indicating a successful
!     constrained curve fit.  The array SDDATA, as input to FC( ), must
!     also be defined with the standard deviation or uncertainty of the
!     Y values to use CV( ).
!
!     To evaluate the variance function after calling FC( ) as stated
!     above, use CV( ) as shown here
!
!          VAR=CV(XVAL,NDATA,NCONST,NORD,NBKPT,BKPT,W)
!
!     The variance function is given by
!
!          VAR=(transpose of B(XVAL))*C*B(XVAL)/MAX(NDATA-N,1)
!
!     where N = NBKPT - NORD.
!
!     The vector B(XVAL) is the B-spline basis function values at
!     X=XVAL.  The covariance matrix, C, of the solution coefficients
!     accounts only for the least squares equations and the explicitly
!     stated equality constraints.  This fact must be considered when
!     interpreting the variance function from a data fitting problem
!     that has inequality constraints on the fitted curve.
!
!     All the variables in the calling sequence for CV( ) are used in
!     FC( ) except the variable XVAL.  Do not change the values of these
!     variables between the call to FC( ) and the use of CV( ).
!
!     The following is a brief description of the variables
!
!     XVAL    The point where the variance is desired.
!
!     NDATA   The number of discrete (X,Y) pairs for which FC( )
!             calculated a piece-wise polynomial curve.
!
!     NCONST  The number of conditions that constrained the B-spline in
!             FC( ).
!
!     NORD    The order of the B-spline used in FC( ).
!             The value of NORD must satisfy 1 < NORD < 20 .
!
!             (The order of the spline is one more than the degree of
!             the piece-wise polynomial defined on each interval.  This
!             is consistent with the B-spline package convention.  For
!             example, NORD=4 when we are using piece-wise cubics.)
!
!     NBKPT   The number of knots in the array BKPT(*).
!             The value of NBKPT must satisfy NBKPT .GE. 2*NORD.
!
!     BKPT(*) The real array of knots.  Normally the problem data
!             interval will be included between the limits BKPT(NORD)
!             and BKPT(NBKPT-NORD+1).  The additional end knots
!             BKPT(I),I=1,...,NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
!             required by FC( ) to compute the functions used to fit
!             the data.
!
!     W(*)    Real work array as used in FC( ).  See FC( ) for the
!             required length of W(*).  The contents of W(*) must not
!             be modified by the user if the variance function is
!             desired.
!
!***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
!                 to discrete data using B-splines, a users guide,
!                 Report SAND78-1291, Sandia Laboratories, December
!                 1978.
!***ROUTINES CALLED  BSPLVN, SDOT
!***REVISION HISTORY  (YYMMDD)
!   780801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CV
      DIMENSION Bkpt(Nbkpt) , W(*) , v(40)
!***FIRST EXECUTABLE STATEMENT  CV
      zero = 0.
      mdg = Nbkpt - Nord + 3
      mdw = Nbkpt - Nord + 1 + Nconst
      is = mdg*(Nord+1) + 2*MAX(Ndata,Nbkpt) + Nbkpt + Nord**2
      last = Nbkpt - Nord + 1
      ileft = Nord
      DO WHILE ( Xval>=Bkpt(ileft+1).AND.ileft<last-1 )
        ileft = ileft + 1
      ENDDO
      CALL BSPLVN(Bkpt,Nord,1,Xval,ileft,v(Nord+1))
      ileft = ileft - Nord + 1
      ip = mdw*(ileft-1) + ileft + is
      n = Nbkpt - Nord
      DO i = 1 , Nord
        v(i) = SDOT(Nord,W(ip),1,v(Nord+1),1)
        ip = ip + mdw
      ENDDO
      CV = MAX(SDOT(Nord,v,1,v(Nord+1),1),zero)
!
!     SCALE THE VARIANCE SO IT IS AN UNBIASED ESTIMATE.
      CV = CV/MAX(Ndata-n,1)
      END FUNCTION CV
