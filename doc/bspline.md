 Documentation for BSPLINE, a package of subprograms for working with piecewise 
 polynomial functions in B-representation.
*****
 ****Category:****  E, E1A, K, Z
*****
 ****Type:****      ALL (BSPDOC-A)
*****
 ****Keywords:****  B-SPLINE, DOCUMENTATION, SPLINES
*****
 ****Author:****  Amos, D. E., (SNLA)
*****
 ****Description:****

     Abstract
         BSPDOC is a non-executable, B-spline documentary routine.
         The narrative describes a B-spline and the routines
         necessary to manipulate B-splines at a fairly high level.
         The basic package described herein is that of reference
         5 with names altered to prevent duplication and conflicts
         with routines from reference 3.  The call lists used here
         are also different.  Work vectors were added to ensure
         portability and proper execution in an overlay environ-
         ment.  These work arrays can be used for other purposes
         except as noted in BSPVN.  While most of the original
         routines in reference 5 were restricted to orders 20
         or less, this restriction was removed from all routines
         except the quadrature routine BSQAD.  (See the section
         below on differentiation and integration for details.)

         The subroutines referenced below are single precision
         routines.  Corresponding double precision versions are also
         part of the package, and these are referenced by prefixing
         a D in front of the single precision name.  For example,
         BVALU and DBVALU are the single and double precision
         versions for evaluating a B-spline or any of its deriva-
         tives in the B-representation.

                ****Description of B-Splines****

     A collection of polynomials of fixed degree K-1 defined on a
     subdivision (X(I),X(I+1)), I=1,...,M-1 of (A,B) with X(1)=A,
     X(M)=B is called a B-spline of order K.  If the spline has K-2
     continuous derivatives on (A,B), then the B-spline is simply
     called a spline of order K.  Each of the M-1 polynomial pieces
     has K coefficients, making a total of K(M-1) parameters.  This
     B-spline and its derivatives have M-2 jumps at the subdivision
     points X(I), I=2,...,M-1.  Continuity requirements at these
     subdivision points add constraints and reduce the number of free
     parameters.  If a B-spline is continuous at each of the M-2 sub-
     division points, there are K(M-1)-(M-2) free parameters; if in
     addition the B-spline has continuous first derivatives, there
     are K(M-1)-2(M-2) free parameters, etc., until we get to a
     spline where we have K(M-1)-(K-1)(M-2) = M+K-2 free parameters.
     Thus, the principle is that increasing the continuity of
     derivatives decreases the number of free parameters and
     conversely.

     The points at which the polynomials are tied together by the
     continuity conditions are called knots.  If two knots are
     allowed to come together at some X(I), then we say that we
     have a knot of multiplicity 2 there, and the knot values are
     the X(I) value.  If we reverse the procedure of the first
     paragraph, we find that adding a knot to increase multiplicity
     increases the number of free parameters and, according to the
     principle above, we thereby introduce a discontinuity in what
     was the highest continuous derivative at that knot.  Thus, the
     number of free parameters is N = NU+K-2 where NU is the sum
     of multiplicities at the X(I) values with X(1) and X(M) of
     multiplicity 1 (NU = M if all knots are simple, i.e., for a
     spline, all knots have multiplicity 1.)  Each knot can have a
     multiplicity of at most K.  A B-spline is commonly written in the
     B-representation

               Y(X) = sum( A(I)*B(I,X), I=1, N)

     to show the explicit dependence of the spline on the free
     parameters or coefficients A(I)=BCOEF(I) and basis functions
     B(I,X).  These basis functions are themselves special B-splines
     which are zero except on (at most) K adjoining intervals where
     each B(I,X) is positive and, in most cases, hat or bell-
     shaped.  In order for the nonzero part of B(1,X) to be a spline
     covering (X(1),X(2)), it is necessary to put K-1 knots to the
     left of A and similarly for B(N,X) to the right of B.  Thus, the
     total number of knots for this representation is NU+2K-2 = N+K.
     These knots are carried in an array T(*) dimensioned by at least
     N+K.  From the construction, A=T(K) and B=T(N+1) and the spline is
     defined on T(K)<=X<=T(N+1).  The nonzero part of each basis
     function lies in the  Interval (T(I),T(I+K)).  In many problems
     where extrapolation beyond A or B is not anticipated, it is common
     practice to set T(1)=T(2)=...=T(K)=A and T(N+1)=T(N+2)=...=
     T(N+K)=B.  In summary, since T(K) and T(N+1) as well as
     interior knots can have multiplicity K, the number of free
     parameters N = sum of multiplicities - K.  The fact that each
     B(I,X) function is nonzero over at most K intervals means that
     for a given X value, there are at most K nonzero terms of the
     sum.  This leads to banded matrices in linear algebra problems,
     and references 3 and 6 take advantage of this in con-
     structing higher level routines to achieve speed and avoid
     ill-conditioning.

                     ****Basic Routines****

     The basic routines which most casual users will need are those
     concerned with direct evaluation of splines or B-splines.
     Since the B-representation, denoted by (T,BCOEF,N,K), is
     preferred because of numerical stability, the knots T(*), the
     B-spline coefficients BCOEF(*), the number of coefficients N,
     and the order K of the polynomial pieces (of degree K-1) are
     usually given.  While the knot array runs from T(1) to T(N+K),
     the B-spline is normally defined on the interval T(K)<=X<=
     T(N+1).  To evaluate the B-spline or any of its derivatives
     on this interval, one can use

                  Y = BVALU(T,BCOEF,N,K,ID,X,INBV,WORK)

     where ID is an integer for the ID-th derivative, 0<=ID<=K-1.
     ID=0 gives the zero-th derivative or B-spline value at X.
     If X<T(K) or X>T(N+1), whether by mistake or the result
     of round off accumulation in incrementing X, BVALU gives a
     diagnostic.  INBV is an initialization parameter which is set
     to 1 on the first call.  Distinct splines require distinct
     INBV parameters.  WORK is a scratch vector of length at least
     3*K.

     When more conventional communication is needed for publication,
     physical interpretation, etc., the B-spline coefficients can
     be converted to piecewise polynomial (PP) coefficients.  Thus,
     the breakpoints (distinct knots) XI(*), the number of
     polynomial pieces LXI, and the (right) derivatives C(*,J) at
     each breakpoint XI(J) are needed to define the Taylor
     expansion to the right of XI(J) on each interval XI(J)<=
     X<XI(J+1), J=1,LXI where XI(1)=A and XI(LXI+1)=B.
     These are obtained from the (T,BCOEF,N,K) representation by

                CALL BSPPP(T,BCOEF,N,K,LDC,C,XI,LXI,WORK)

     where LDC>=K is the leading dimension of the matrix C and
     WORK is a scratch vector of length at least K*(N+3).
     Then the PP-representation (C,XI,LXI,K) of Y(X), denoted
     by Y(J,X) on each interval XI(J)<=X<XI(J+1), is

     Y(J,X) = sum( C(I,J)*((X-XI(J))****(I-1))/factorial(I-1), I=1,K)

     for J=1,...,LXI.  One must view this conversion from the B-
     to the PP-representation with some skepticism because the
     conversion may lose significant digits when the B-spline
     varies in an almost discontinuous fashion.  To evaluate
     the B-spline or any of its derivatives using the PP-
     representation, one uses

                Y = PPVAL(LDC,C,XI,LXI,K,ID,X,INPPV)

     where ID and INPPV have the same meaning and usage as ID and
     INBV in BVALU.

     To determine to what extent the conversion process loses
     digits, compute the relative error ABS((Y1-Y2)/Y2) over
     the X interval with Y1 from PPVAL and Y2 from BVALU.  A
     major reason for considering PPVAL is that evaluation is
     much faster than that from BVALU.

     Recall that when multiple knots are encountered, jump type
     discontinuities in the B-spline or its derivatives occur
     at these knots, and we need to know that BVALU and PPVAL
     return right limiting values at these knots except at
     X=B where left limiting values are returned.  These values
     are used for the Taylor expansions about left end points of
     breakpoint intervals.  That is, the derivatives C(*,J) are
     right derivatives.  Note also that a computed X value which,
     mathematically, would be a knot value may differ from the knot
     by a round off error.  When this happens in evaluating a dis-
     continuous B-spline or some discontinuous derivative, the
     value at the knot and the value at X can be radically
     different.  In this case, setting X to a T or XI value makes
     the computation precise.  For left limiting values at knots
     other than X=B, see the prologues to BVALU and other
     routines.

                     ****Interpolation****

     BINTK is used to generate B-spline parameters (T,BCOEF,N,K)
     which will interpolate the data by calls to BVALU.  A similar
     interpolation can also be done for cubic splines using BINT4
     or the code in reference 7.  If the PP-representation is given,
     one can evaluate this representation at an appropriate number of
     abscissas to create data then use BINTK or BINT4 to generate
     the B-representation.

               ****Differentiation and Integration****

     Derivatives of B-splines are obtained from BVALU or PPVAL.
     Integrals are obtained from BSQAD using the B-representation
     (T,BCOEF,N,K) and PPQAD using the PP-representation (C,XI,LXI,
     K).  More complicated integrals involving the product of a
     of a function F and some derivative of a B-spline can be
     evaluated with BFQAD or PFQAD using the B- or PP- represen-
     tations respectively.  All quadrature routines, except for PPQAD,
     are limited in accuracy to 18 digits or working precision,
     whichever is smaller.  PPQAD is limited to working precision
     only.  In addition, the order K for BSQAD is limited to 20 or
     less.  If orders greater than 20 are required, use BFQAD with
     F(X) = 1.

                      ****Extrapolation****

     Extrapolation outside the interval (A,B) can be accomplished
     easily by the PP-representation using PPVAL.  However,
     caution should be exercised, especially when several knots
     are located at A or B or when the extrapolation is carried
     significantly beyond A or B.  On the other hand, direct
     evaluation with BVALU outside A=T(K)<=X<=T(N+1)=B
     produces an error message, and some manipulation of the knots
     and coefficients are needed to extrapolate with BVALU.  This
     process is described in reference 6.

                ****Curve Fitting and Smoothing****

     Unless one has many accurate data points, direct inter-
     polation is not recommended for summarizing data.  The
     results are often not in accordance with intuition since the
     fitted curve tends to oscillate through the set of points.
     Monotone splines (reference 7) can help curb this undulating
     tendency but constrained least squares is more likely to give an
     acceptable fit with fewer parameters.  Subroutine FC, des-
     cribed in reference 6, is recommended for this purpose.  The
     output from this fitting process is the B-representation.

              **** Routines in the B-Spline Package ****

                      Single Precision Routines

         The subroutines referenced below are SINGLE PRECISION
         routines. Corresponding DOUBLE PRECISION versions are also
         part of the package and these are referenced by prefixing
         a D in front of the single precision name. For example,
         BVALU and DBVALU are the SINGLE and DOUBLE PRECISION
         versions for evaluating a B-spline or any of its deriva-
         tives in the B-representation.

     BINT4 - interpolates with splines of order 4
     BINTK - interpolates with splines of order k
     BSQAD - integrates the B-representation on subintervals
     PPQAD - integrates the PP-representation
     BFQAD - integrates the product of a function F and any spline
             derivative in the B-representation
     PFQAD - integrates the product of a function F and any spline
             derivative in the PP-representation
     BVALU - evaluates the B-representation or a derivative
     PPVAL - evaluates the PP-representation or a derivative
     INTRV - gets the largest index of the knot to the left of x
     BSPPP - converts from B- to PP-representation
     BSPVD - computes nonzero basis functions and derivatives at x
     BSPDR - sets up difference array for BSPEV
     BSPEV - evaluates the B-representation and derivatives
     BSPVN - called by BSPEV, BSPVD, BSPPP and BINTK for function and
             derivative evaluations
                        Auxiliary Routines

       BSGQ8,PPGQ8,BNSLV,BNFAC,XERMSG,DBSGQ8,DPPGQ8,DBNSLV,DBNFAC

                    Machine Dependent Routines

                      I1MACH, R1MACH, D1MACH

*****
 ****References:****   
    1. D. E. Amos, Computation with splines and B-splines, Report SAND78-1968, Sandia Laboratories, March 1979.   
    2. D. E. Amos, Quadrature subroutines for splines and B-splines, Report SAND79-1825, Sandia Laboratories, December 1979.   
    3. Carl de Boor, A Practical Guide to Splines, Applied Mathematics Series 27, Springer-Verlag, New York, 1978.   
    4. Carl de Boor, On calculating with B-Splines, Journal of Approximation Theory 6, (1972), pp. 50-62.   
    5. Carl de Boor, Package for calculating with B-splines, SIAM Journal on Numerical Analysis 14, 3 (June 1977), pp. 441-472.   
    6. R. J. Hanson, Constrained least squares curve fitting to discrete data using B-splines, a users guide,  Report SAND78-1291, Sandia Laboratories, December 1978.   
    7. F. N. Fritsch and R. E. Carlson, Monotone piecewise cubic interpolation, SIAM Journal on Numerical Analysis 17, 2 (April 1980), pp. 238-246.