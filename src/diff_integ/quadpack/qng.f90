!** QNG
SUBROUTINE QNG(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier)
  !>
  !  The routine calculates an approximation result to a
  !            given definite integral I = integral of F over (A,B),
  !            hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A1
  !***
  ! **Type:**      SINGLE PRECISION (QNG-S, DQNG-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, GAUSS-KRONROD(PATTERSON) RULES,
  !             NONADAPTIVE, QUADPACK, QUADRATURE, SMOOTH INTEGRAND
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **Description:**
  !
  ! NON-ADAPTIVE INTEGRATION
  ! STANDARD FORTRAN SUBROUTINE
  ! REAL VERSION
  !
  !           F      - Real version
  !                    Function subprogram defining the integrand function
  !                    F(X). The actual name for F needs to be declared
  !                    E X T E R N A L in the driver program.
  !
  !           A      - Real version
  !                    Lower limit of integration
  !
  !           B      - Real version
  !                    Upper limit of integration
  !
  !           EPSABS - Real
  !                    Absolute accuracy requested
  !           EPSREL - Real
  !                    Relative accuracy requested
  !                    If  EPSABS.LE.0
  !                    And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                    The routine will end with IER = 6.
  !
  !         ON RETURN
  !           RESULT - Real
  !                    Approximation to the integral I
  !                    Result is obtained by applying the 21-POINT
  !                    GAUSS-KRONROD RULE (RES21) obtained by optimal
  !                    addition of abscissae to the 10-POINT GAUSS RULE
  !                    (RES10), or by applying the 43-POINT RULE (RES43)
  !                    obtained by optimal addition of abscissae to the
  !                    21-POINT GAUSS-KRONROD RULE, or by applying the
  !                    87-POINT RULE (RES87) obtained by optimal addition
  !                    of abscissae to the 43-POINT RULE.
  !
  !           ABSERR - Real
  !                    Estimate of the modulus of the absolute error,
  !                    which should EQUAL or EXCEED ABS(I-RESULT)
  !
  !           NEVAL  - Integer
  !                    Number of integrand evaluations
  !
  !           IER    - IER = 0 normal and reliable termination of the
  !                            routine. It is assumed that the requested
  !                            accuracy has been achieved.
  !                    IER.GT.0 Abnormal termination of the routine. It is
  !                            assumed that the requested accuracy has
  !                            not been achieved.
  !           ERROR MESSAGES
  !                    IER = 1 The maximum number of steps has been
  !                            executed. The integral is probably too
  !                            difficult to be calculated by DQNG.
  !                        = 6 The input is invalid, because
  !                            EPSABS.LE.0 AND
  !                            EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28).
  !                            RESULT, ABSERR and NEVAL are set to zero.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : XERMSG, R1MACH
  INTERFACE
    REAL(SP) FUNCTION F(X)
      IMPORT SP
      REAL(SP) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Ier, Neval
  REAL(SP) :: A, Abserr, B, Epsabs, Epsrel, Result
  INTEGER :: ipx, k, l
  REAL(SP) :: absc, centr, dhlgth, epmach, fcentr, fval, fval1, fval2, fv1(5), &
    fv2(5), fv3(5), fv4(5), hlgth, res10, res21, res43, res87, resabs, resasc, &
    reskh, savfun(21), uflow
  !
  !           THE FOLLOWING DATA STATEMENTS CONTAIN THE
  !           ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED.
  !
  !           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43-
  !                   AND 87-POINT RULE
  !           X2      ABSCISSAE COMMON TO THE 21-, 43- AND
  !                   87-POINT RULE
  !           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT
  !                   RULE
  !           X4      ABSCISSAE OF THE 87-POINT RULE
  !           W10     WEIGHTS OF THE 10-POINT FORMULA
  !           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR
  !                   ABSCISSAE X1
  !           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR
  !                   ABSCISSAE X2
  !           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR
  !                   ABSCISSAE X1, X3
  !           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR
  !                   ABSCISSAE X3
  !           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR
  !                   ABSCISSAE X1, X2, X3
  !           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR
  !                   ABSCISSAE X4
  !
  REAL(SP), PARAMETER :: x1(5) = [ 0.9739065285171717E+00, 0.8650633666889845E+00, &
    0.6794095682990244E+00, 0.4333953941292472E+00, 0.1488743389816312E+00 ]
  REAL(SP), PARAMETER :: x2(5) = [ 0.9956571630258081E+00, 0.9301574913557082E+00, &
    0.7808177265864169E+00, 0.5627571346686047E+00, 0.2943928627014602E+00 ]
  REAL(SP), PARAMETER :: x3(11) = [ 0.9993333609019321E+00, 0.9874334029080889E+00, &
    0.9548079348142663E+00, 0.9001486957483283E+00, 0.8251983149831142E+00, &
    0.7321483889893050E+00, 0.6228479705377252E+00, 0.4994795740710565E+00, &
    0.3649016613465808E+00, 0.2222549197766013E+00, 0.7465061746138332E-01 ]
  REAL(SP), PARAMETER :: x4(22) = [ 0.9999029772627292E+00, 0.9979898959866787E+00, &
    0.9921754978606872E+00, 0.9813581635727128E+00, 0.9650576238583846E+00, &
    0.9431676131336706E+00, 0.9158064146855072E+00, 0.8832216577713165E+00, &
    0.8457107484624157E+00, 0.8035576580352310E+00, 0.7570057306854956E+00, &
    0.7062732097873218E+00, 0.6515894665011779E+00, 0.5932233740579611E+00, &
    0.5314936059708319E+00, 0.4667636230420228E+00, 0.3994248478592188E+00, &
    0.3298748771061883E+00, 0.2585035592021616E+00, 0.1856953965683467E+00, &
    0.1118422131799075E+00, 0.3735212339461987E-01 ]
  REAL(SP), PARAMETER :: w10(5) = [ 0.6667134430868814E-01, 0.1494513491505806E+00, &
    0.2190863625159820E+00, 0.2692667193099964E+00, 0.2955242247147529E+00 ]
  REAL(SP), PARAMETER :: w21a(5) = [ 0.3255816230796473E-01, 0.7503967481091995E-01, &
    0.1093871588022976E+00, 0.1347092173114733E+00, 0.1477391049013385E+00 ]
  REAL(SP), PARAMETER :: w21b(6) = [ 0.1169463886737187E-01, 0.5475589657435200E-01, &
    0.9312545458369761E-01, 0.1234919762620659E+00, 0.1427759385770601E+00, &
    0.1494455540029169E+00 ]
  REAL(SP), PARAMETER :: w43a(10) = [ 0.1629673428966656E-01, 0.3752287612086950E-01, &
    0.5469490205825544E-01, 0.6735541460947809E-01, 0.7387019963239395E-01, &
    0.5768556059769796E-02, 0.2737189059324884E-01, 0.4656082691042883E-01, &
    0.6174499520144256E-01, 0.7138726726869340E-01 ]
  REAL(SP), PARAMETER :: w43b(12) = [ 0.1844477640212414E-02, 0.1079868958589165E-01, &
    0.2189536386779543E-01, 0.3259746397534569E-01, 0.4216313793519181E-01, &
    0.5074193960018458E-01, 0.5837939554261925E-01, 0.6474640495144589E-01, &
    0.6956619791235648E-01, 0.7282444147183321E-01, 0.7450775101417512E-01, &
    0.7472214751740301E-01 ]
  REAL(SP), PARAMETER :: w87a(21) = [ 0.8148377384149173E-02, 0.1876143820156282E-01, &
    0.2734745105005229E-01, 0.3367770731163793E-01, 0.3693509982042791E-01, &
    0.2884872430211531E-02, 0.1368594602271270E-01, 0.2328041350288831E-01, &
    0.3087249761171336E-01, 0.3569363363941877E-01, 0.9152833452022414E-03, &
    0.5399280219300471E-02, 0.1094767960111893E-01, 0.1629873169678734E-01, &
    0.2108156888920384E-01, 0.2537096976925383E-01, 0.2918969775647575E-01, &
    0.3237320246720279E-01, 0.3478309895036514E-01, 0.3641222073135179E-01, &
    0.3725387550304771E-01 ]
  REAL(SP), PARAMETER :: w87b(23) = [ 0.2741455637620724E-03, 0.1807124155057943E-02, &
    0.4096869282759165E-02, 0.6758290051847379E-02, 0.9549957672201647E-02, &
    0.1232944765224485E-01, 0.1501044734638895E-01, 0.1754896798624319E-01, &
    0.1993803778644089E-01, 0.2219493596101229E-01, 0.2433914712600081E-01, &
    0.2637450541483921E-01, 0.2828691078877120E-01, 0.3005258112809270E-01, &
    0.3164675137143993E-01, 0.3305041341997850E-01, 0.3425509970422606E-01, &
    0.3526241266015668E-01, 0.3607698962288870E-01, 0.3669860449845609E-01, &
    0.3712054926983258E-01, 0.3733422875193504E-01, 0.3736107376267902E-01 ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
  !           FCENTR - FUNCTION VALUE AT MID POINT
  !           ABSC   - ABSCISSA
  !           FVAL   - FUNCTION VALUE
  !           SAVFUN - ARRAY OF FUNCTION VALUES WHICH
  !                    HAVE ALREADY BEEN COMPUTED
  !           RES10  - 10-POINT GAUSS RESULT
  !           RES21  - 21-POINT KRONROD RESULT
  !           RES43  - 43-POINT RESULT
  !           RES87  - 87-POINT RESULT
  !           RESABS - APPROXIMATION TO THE INTEGRAL OF ABS(F)
  !           RESASC - APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  QNG
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  Result = 0.0E+00
  Abserr = 0.0E+00
  Neval = 0
  Ier = 6
  IF ( Epsabs>0.0E+00.OR.Epsrel>=MAX(0.5E-14,0.5E+02*epmach) ) THEN
    hlgth = 0.5E+00*(B-A)
    dhlgth = ABS(hlgth)
    centr = 0.5E+00*(B+A)
    fcentr = F(centr)
    Neval = 21
    Ier = 1
    !
    !          COMPUTE THE INTEGRAL USING THE 10- AND 21-POINT FORMULA.
    !
    DO l = 1, 3
      SELECT CASE (l)
        CASE (2)
          !
          !          COMPUTE THE INTEGRAL USING THE 43-POINT FORMULA.
          !
          res43 = w43b(12)*fcentr
          Neval = 43
          DO k = 1, 10
            res43 = res43 + savfun(k)*w43a(k)
          END DO
          DO k = 1, 11
            ipx = ipx + 1
            absc = hlgth*x3(k)
            fval = F(absc+centr) + F(centr-absc)
            res43 = res43 + fval*w43b(k)
            savfun(ipx) = fval
          END DO
          !
          !          TEST FOR CONVERGENCE.
          !
          Result = res43*hlgth
          Abserr = ABS((res43-res21)*hlgth)
        CASE (3)
          !
          !          COMPUTE THE INTEGRAL USING THE 87-POINT FORMULA.
          !
          res87 = w87b(23)*fcentr
          Neval = 87
          DO k = 1, 21
            res87 = res87 + savfun(k)*w87a(k)
          END DO
          DO k = 1, 22
            absc = hlgth*x4(k)
            res87 = res87 + w87b(k)*(F(absc+centr)+F(centr-absc))
          END DO
          Result = res87*hlgth
          Abserr = ABS((res87-res43)*hlgth)
        CASE DEFAULT
          res10 = 0.0E+00
          res21 = w21b(6)*fcentr
          resabs = w21b(6)*ABS(fcentr)
          DO k = 1, 5
            absc = hlgth*x1(k)
            fval1 = F(centr+absc)
            fval2 = F(centr-absc)
            fval = fval1 + fval2
            res10 = res10 + w10(k)*fval
            res21 = res21 + w21a(k)*fval
            resabs = resabs + w21a(k)*(ABS(fval1)+ABS(fval2))
            savfun(k) = fval
            fv1(k) = fval1
            fv2(k) = fval2
          END DO
          ipx = 5
          DO k = 1, 5
            ipx = ipx + 1
            absc = hlgth*x2(k)
            fval1 = F(centr+absc)
            fval2 = F(centr-absc)
            fval = fval1 + fval2
            res21 = res21 + w21b(k)*fval
            resabs = resabs + w21b(k)*(ABS(fval1)+ABS(fval2))
            savfun(ipx) = fval
            fv3(k) = fval1
            fv4(k) = fval2
          END DO
          !
          !          TEST FOR CONVERGENCE.
          !
          Result = res21*hlgth
          resabs = resabs*dhlgth
          reskh = 0.5E+00*res21
          resasc = w21b(6)*ABS(fcentr-reskh)
          DO k = 1, 5
            resasc = resasc + w21a(k)*(ABS(fv1(k)-reskh)+ABS(fv2(k)-reskh))&
              + w21b(k)*(ABS(fv3(k)-reskh)+ABS(fv4(k)-reskh))
          END DO
          Abserr = ABS((res21-res10)*hlgth)
          resasc = resasc*dhlgth
      END SELECT
      IF ( resasc/=0.0E+00.AND.Abserr/=0.0E+00 )&
        Abserr = resasc*MIN(0.1E+01,(0.2E+03*Abserr/resasc)**1.5E+00)
      IF ( resabs>uflow/(0.5E+02*epmach) )&
        Abserr = MAX((epmach*0.5E+02)*resabs,Abserr)
      IF ( Abserr<=MAX(Epsabs,Epsrel*ABS(Result)) ) Ier = 0
      !- **JUMP OUT OF DO-LOOP
      IF ( Ier==0 ) RETURN
    END DO
  END IF
  CALL XERMSG('QNG','ABNORMAL RETURN',Ier,0)
  RETURN
END SUBROUTINE QNG
