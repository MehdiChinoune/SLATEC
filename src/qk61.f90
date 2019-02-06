!*==QK61.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QK61
      SUBROUTINE QK61(F,A,B,Result,Abserr,Resabs,Resasc)
      IMPLICIT NONE
!*--QK615
!***BEGIN PROLOGUE  QK61
!***PURPOSE  To compute I = Integral of F over (A,B) with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      SINGLE PRECISION (QK61-S, DQK61-D)
!***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Integration rule
!        Standard fortran subroutine
!        Real version
!
!
!        PARAMETERS
!         ON ENTRY
!           F      - Real
!                    Function subprogram defining the integrand
!                    function F(X). The actual name for F needs to be
!                    declared E X T E R N A L in the calling program.
!
!           A      - Real
!                    Lower limit of integration
!
!           B      - Real
!                    Upper limit of integration
!
!         ON RETURN
!           RESULT - Real
!                    Approximation to the integral I
!                    RESULT is computed by applying the 61-point
!                    Kronrod rule (RESK) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (RESG).
!
!           ABSERR - Real
!                    Estimate of the modulus of the absolute error,
!                    which should equal or exceed ABS(I-RESULT)
!
!           RESABS - Real
!                    Approximation to the integral J
!
!           RESASC - Real
!                    Approximation to the integral of ABS(F-I/(B-A))
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QK61
!
      REAL A , absc , Abserr , B , centr , dhlgth , epmach , F , fc , fsum , 
     &     fval1 , fval2 , fv1 , fv2 , hlgth , Resabs , Resasc , resg , resk , 
     &     reskh , Result , R1MACH , uflow , wg , wgk , xgk
      INTEGER j , jtw , jtwm1
      EXTERNAL F
!
      DIMENSION fv1(30) , fv2(30) , xgk(31) , wgk(31) , wg(15)
!
!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
!           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
!           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
!                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
!                   GAUSS RULE
!                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
!                   TO THE 30-POINT GAUSS RULE
!
!           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
!
!           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
!
      SAVE xgk , wgk , wg
      DATA xgk(1) , xgk(2) , xgk(3) , xgk(4) , xgk(5) , xgk(6) , xgk(7) , xgk(8)
     &     , xgk(9) , xgk(10)/0.9994844100504906E+00 , 0.9968934840746495E+00 , 
     &     0.9916309968704046E+00 , 0.9836681232797472E+00 , 
     &     0.9731163225011263E+00 , 0.9600218649683075E+00 , 
     &     0.9443744447485600E+00 , 0.9262000474292743E+00 , 
     &     0.9055733076999078E+00 , 0.8825605357920527E+00/
      DATA xgk(11) , xgk(12) , xgk(13) , xgk(14) , xgk(15) , xgk(16) , xgk(17) , 
     &     xgk(18) , xgk(19) , xgk(20)/0.8572052335460611E+00 , 
     &     0.8295657623827684E+00 , 0.7997278358218391E+00 , 
     &     0.7677774321048262E+00 , 0.7337900624532268E+00 , 
     &     0.6978504947933158E+00 , 0.6600610641266270E+00 , 
     &     0.6205261829892429E+00 , 0.5793452358263617E+00 , 
     &     0.5366241481420199E+00/
      DATA xgk(21) , xgk(22) , xgk(23) , xgk(24) , xgk(25) , xgk(26) , xgk(27) , 
     &     xgk(28) , xgk(29) , xgk(30) , xgk(31)/0.4924804678617786E+00 , 
     &     0.4470337695380892E+00 , 0.4004012548303944E+00 , 
     &     0.3527047255308781E+00 , 0.3040732022736251E+00 , 
     &     0.2546369261678898E+00 , 0.2045251166823099E+00 , 
     &     0.1538699136085835E+00 , 0.1028069379667370E+00 , 
     &     0.5147184255531770E-01 , 0.0E+00/
      DATA wgk(1) , wgk(2) , wgk(3) , wgk(4) , wgk(5) , wgk(6) , wgk(7) , wgk(8)
     &     , wgk(9) , wgk(10)/0.1389013698677008E-02 , 0.3890461127099884E-02 , 
     &     0.6630703915931292E-02 , 0.9273279659517763E-02 , 
     &     0.1182301525349634E-01 , 0.1436972950704580E-01 , 
     &     0.1692088918905327E-01 , 0.1941414119394238E-01 , 
     &     0.2182803582160919E-01 , 0.2419116207808060E-01/
      DATA wgk(11) , wgk(12) , wgk(13) , wgk(14) , wgk(15) , wgk(16) , wgk(17) , 
     &     wgk(18) , wgk(19) , wgk(20)/0.2650995488233310E-01 , 
     &     0.2875404876504129E-01 , 0.3090725756238776E-01 , 
     &     0.3298144705748373E-01 , 0.3497933802806002E-01 , 
     &     0.3688236465182123E-01 , 0.3867894562472759E-01 , 
     &     0.4037453895153596E-01 , 0.4196981021516425E-01 , 
     &     0.4345253970135607E-01/
      DATA wgk(21) , wgk(22) , wgk(23) , wgk(24) , wgk(25) , wgk(26) , wgk(27) , 
     &     wgk(28) , wgk(29) , wgk(30) , wgk(31)/0.4481480013316266E-01 , 
     &     0.4605923827100699E-01 , 0.4718554656929915E-01 , 
     &     0.4818586175708713E-01 , 0.4905543455502978E-01 , 
     &     0.4979568342707421E-01 , 0.5040592140278235E-01 , 
     &     0.5088179589874961E-01 , 0.5122154784925877E-01 , 
     &     0.5142612853745903E-01 , 0.5149472942945157E-01/
      DATA wg(1) , wg(2) , wg(3) , wg(4) , wg(5) , wg(6) , wg(7) , 
     &     wg(8)/0.7968192496166606E-02 , 0.1846646831109096E-01 , 
     &     0.2878470788332337E-01 , 0.3879919256962705E-01 , 
     &     0.4840267283059405E-01 , 0.5749315621761907E-01 , 
     &     0.6597422988218050E-01 , 0.7375597473770521E-01/
      DATA wg(9) , wg(10) , wg(11) , wg(12) , wg(13) , wg(14) , wg(15)
     &     /0.8075589522942022E-01 , 0.8689978720108298E-01 , 
     &     0.9212252223778613E-01 , 0.9636873717464426E-01 , 
     &     0.9959342058679527E-01 , 0.1017623897484055E+00 , 
     &     0.1028526528935588E+00/
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 30-POINT GAUSS RULE
!           RESK   - RESULT OF THE 61-POINT KRONROD RULE
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
!                    OVER (A,B), I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QK61
      epmach = R1MACH(4)
      uflow = R1MACH(1)
!
      centr = 0.5E+00*(B+A)
      hlgth = 0.5E+00*(B-A)
      dhlgth = ABS(hlgth)
!
!           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
!           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
      resg = 0.0E+00
      fc = F(centr)
      resk = wgk(31)*fc
      Resabs = ABS(resk)
      DO j = 1 , 15
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = F(centr-absc)
        fval2 = F(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1 + fval2
        resg = resg + wg(j)*fsum
        resk = resk + wgk(jtw)*fsum
        Resabs = Resabs + wgk(jtw)*(ABS(fval1)+ABS(fval2))
      ENDDO
      DO j = 1 , 15
        jtwm1 = j*2 - 1
        absc = hlgth*xgk(jtwm1)
        fval1 = F(centr-absc)
        fval2 = F(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1 + fval2
        resk = resk + wgk(jtwm1)*fsum
        Resabs = Resabs + wgk(jtwm1)*(ABS(fval1)+ABS(fval2))
      ENDDO
      reskh = resk*0.5E+00
      Resasc = wgk(31)*ABS(fc-reskh)
      DO j = 1 , 30
        Resasc = Resasc + wgk(j)*(ABS(fv1(j)-reskh)+ABS(fv2(j)-reskh))
      ENDDO
      Result = resk*hlgth
      Resabs = Resabs*dhlgth
      Resasc = Resasc*dhlgth
      Abserr = ABS((resk-resg)*hlgth)
      IF ( Resasc/=0.0E+00.AND.Abserr/=0.0E+00 )
     &     Abserr = Resasc*MIN(0.1E+01,(0.2E+03*Abserr/Resasc)**1.5E+00)
      IF ( Resabs>uflow/(0.5E+02*epmach) ) Abserr = MAX((epmach*0.5E+02)*Resabs,
     &     Abserr)
      END SUBROUTINE QK61
