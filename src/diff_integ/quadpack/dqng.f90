!** DQNG
PURE SUBROUTINE DQNG(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier)
  !> The routine calculates an approximation result to a given definite integral
  !  I = integral of F over (A,B), hopefully satisfying following claim for accuracy
  !  ABS(I-RESULT)<=MAX(EPSABS,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A1
  !***
  ! **Type:**      DOUBLE PRECISION (QNG-S, DQNG-D)
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
  ! DOUBLE PRECISION VERSION
  !
  !           F      - Double precision
  !                    Function subprogram defining the integrand function
  !                    F(X). The actual name for F needs to be declared
  !                    E X T E R N A L in the driver program.
  !
  !           A      - Double precision
  !                    Lower limit of integration
  !
  !           B      - Double precision
  !                    Upper limit of integration
  !
  !           EPSABS - Double precision
  !                    Absolute accuracy requested
  !           EPSREL - Double precision
  !                    Relative accuracy requested
  !                    If  EPSABS<=0
  !                    And EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
  !                    The routine will end with IER = 6.
  !
  !         ON RETURN
  !           RESULT - Double precision
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
  !           ABSERR - Double precision
  !                    Estimate of the modulus of the absolute error,
  !                    which should EQUAL or EXCEED ABS(I-RESULT)
  !
  !           NEVAL  - Integer
  !                    Number of integrand evaluations
  !
  !           IER    - IER = 0 normal and reliable termination of the
  !                            routine. It is assumed that the requested
  !                            accuracy has been achieved.
  !                    IER>0 Abnormal termination of the routine. It is
  !                            assumed that the requested accuracy has
  !                            not been achieved.
  !           ERROR MESSAGES
  !                    IER = 1 The maximum number of steps has been
  !                            executed. The integral is probably too
  !                            difficult to be calculated by DQNG.
  !                        = 6 The input is invalid, because
  !                            EPSABS<=0 AND
  !                            EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28).
  !                            RESULT, ABSERR and NEVAL are set to zero.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : tiny_dp, eps_dp
  !
  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(OUT) :: Ier, Neval
  REAL(DP), INTENT(IN) :: A, B, Epsabs, Epsrel
  REAL(DP), INTENT(OUT) :: Abserr, Result
  !
  INTEGER :: ipx, k, l
  REAL(DP) :: absc, centr, dhlgth, epmach, fcentr, fval, fval1, fval2, fv1(5), &
    fv2(5), fv3(5), fv4(5), hlgth, res10, res21, res43, res87, resabs, resasc, &
    reskh, savfun(21), uflow
  !
  !           THE FOLLOWING DATA STATEMENTS CONTAIN THE
  !           ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED.
  !
  !           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43- AND 87-
  !                   POINT RULE
  !           X2      ABSCISSAE COMMON TO THE 21-, 43- AND 87-POINT RULE
  !           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT RULE
  !           X4      ABSCISSAE OF THE 87-POINT RULE
  !           W10     WEIGHTS OF THE 10-POINT FORMULA
  !           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X1
  !           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X2
  !           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X1, X3
  !           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X3
  !           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X1,
  !                   X2, X3
  !           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X4
  !
  !
  ! GAUSS-KRONROD-PATTERSON QUADRATURE COEFFICIENTS FOR USE IN
  ! QUADPACK ROUTINE QNG.  THESE COEFFICIENTS WERE CALCULATED WITH
  ! 101 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, BELL LABS, NOV 1981.
  !
  REAL(DP), PARAMETER :: x1(5) = [ 0.973906528517171720077964012084452_DP, &
    0.865063366688984510732096688423493_DP, 0.679409568299024406234327365114874_DP, &
    0.433395394129247190799265943165784_DP, 0.148874338981631210884826001129720_DP ]
  REAL(DP), PARAMETER :: w10(5) = [ 0.066671344308688137593568809893332_DP, &
    0.149451349150580593145776339657697_DP, 0.219086362515982043995534934228163_DP, &
    0.269266719309996355091226921569469_DP, 0.295524224714752870173892994651338_DP ]
  REAL(DP), PARAMETER :: x2(5) = [ 0.995657163025808080735527280689003_DP, &
    0.930157491355708226001207180059508_DP, 0.780817726586416897063717578345042_DP, &
    0.562757134668604683339000099272694_DP, 0.294392862701460198131126603103866_DP ]
  REAL(DP), PARAMETER :: w21a(5) = [ 0.032558162307964727478818972459390_DP, &
    0.075039674810919952767043140916190_DP, 0.109387158802297641899210590325805_DP, &
    0.134709217311473325928054001771707_DP, 0.147739104901338491374841515972068_DP ]
  REAL(DP), PARAMETER :: w21b(6) = [ 0.011694638867371874278064396062192_DP, &
    0.054755896574351996031381300244580_DP, 0.093125454583697605535065465083366_DP, &
    0.123491976262065851077958109831074_DP, 0.142775938577060080797094273138717_DP, &
    0.149445554002916905664936468389821_DP ]
  REAL(DP), PARAMETER :: x3(11) = [ 0.999333360901932081394099323919911_DP, &
    0.987433402908088869795961478381209_DP, 0.954807934814266299257919200290473_DP, &
    0.900148695748328293625099494069092_DP, 0.825198314983114150847066732588520_DP, &
    0.732148388989304982612354848755461_DP, 0.622847970537725238641159120344323_DP, &
    0.499479574071056499952214885499755_DP, 0.364901661346580768043989548502644_DP, &
    0.222254919776601296498260928066212_DP, 0.074650617461383322043914435796506_DP ]
  REAL(DP), PARAMETER :: w43a(10) = [ 0.016296734289666564924281974617663_DP, &
    0.037522876120869501461613795898115_DP, 0.054694902058255442147212685465005_DP, &
    0.067355414609478086075553166302174_DP, 0.073870199632393953432140695251367_DP, &
    0.005768556059769796184184327908655_DP, 0.027371890593248842081276069289151_DP, &
    0.046560826910428830743339154433824_DP, 0.061744995201442564496240336030883_DP, &
    0.071387267268693397768559114425516_DP ]
  REAL(DP), PARAMETER :: w43b(12) = [ 0.001844477640212414100389106552965_DP, &
    0.010798689585891651740465406741293_DP, 0.021895363867795428102523123075149_DP, &
    0.032597463975345689443882222526137_DP, 0.042163137935191811847627924327955_DP, &
    0.050741939600184577780189020092084_DP, 0.058379395542619248375475369330206_DP, &
    0.064746404951445885544689259517511_DP, 0.069566197912356484528633315038405_DP, &
    0.072824441471833208150939535192842_DP, 0.074507751014175118273571813842889_DP, &
    0.074722147517403005594425168280423_DP ]
  REAL(DP), PARAMETER :: x4(22) = [ 0.999902977262729234490529830591582_DP, &
    0.997989895986678745427496322365960_DP, 0.992175497860687222808523352251425_DP, &
    0.981358163572712773571916941623894_DP, 0.965057623858384619128284110607926_DP, &
    0.943167613133670596816416634507426_DP, 0.915806414685507209591826430720050_DP, &
    0.883221657771316501372117548744163_DP, 0.845710748462415666605902011504855_DP, &
    0.803557658035230982788739474980964_DP, 0.757005730685495558328942793432020_DP, &
    0.706273209787321819824094274740840_DP, 0.651589466501177922534422205016736_DP, &
    0.593223374057961088875273770349144_DP, 0.531493605970831932285268948562671_DP, &
    0.466763623042022844871966781659270_DP, 0.399424847859218804732101665817923_DP, &
    0.329874877106188288265053371824597_DP, 0.258503559202161551802280975429025_DP, &
    0.185695396568346652015917141167606_DP, 0.111842213179907468172398359241362_DP, &
    0.037352123394619870814998165437704_DP ]
  REAL(DP), PARAMETER :: w87a(21) = [ 0.008148377384149172900002878448190_DP, &
    0.018761438201562822243935059003794_DP, 0.027347451050052286161582829741283_DP, &
    0.033677707311637930046581056957588_DP, 0.036935099820427907614589586742499_DP, &
    0.002884872430211530501334156248695_DP, 0.013685946022712701888950035273128_DP, &
    0.023280413502888311123409291030404_DP, 0.030872497611713358675466394126442_DP, &
    0.035693633639418770719351355457044_DP, 0.000915283345202241360843392549948_DP, &
    0.005399280219300471367738743391053_DP, 0.010947679601118931134327826856808_DP, &
    0.016298731696787335262665703223280_DP, 0.021081568889203835112433060188190_DP, &
    0.025370969769253827243467999831710_DP, 0.029189697756475752501446154084920_DP, &
    0.032373202467202789685788194889595_DP, 0.034783098950365142750781997949596_DP, &
    0.036412220731351787562801163687577_DP, 0.037253875503047708539592001191226_DP ]
  REAL(DP), PARAMETER :: w87b(23) = [ 0.000274145563762072350016527092881_DP, &
    0.001807124155057942948341311753254_DP, 0.004096869282759164864458070683480_DP, &
    0.006758290051847378699816577897424_DP, 0.009549957672201646536053581325377_DP, &
    0.012329447652244853694626639963780_DP, 0.015010447346388952376697286041943_DP, &
    0.017548967986243191099665352925900_DP, 0.019938037786440888202278192730714_DP, &
    0.022194935961012286796332102959499_DP, 0.024339147126000805470360647041454_DP, &
    0.026374505414839207241503786552615_DP, 0.028286910788771200659968002987960_DP, &
    0.030052581128092695322521110347341_DP, 0.031646751371439929404586051078883_DP, &
    0.033050413419978503290785944862689_DP, 0.034255099704226061787082821046821_DP, &
    0.035262412660156681033782717998428_DP, 0.036076989622888701185500318003895_DP, &
    0.036698604498456094498018047441094_DP, 0.037120549269832576114119958413599_DP, &
    0.037334228751935040321235449094698_DP, 0.037361073762679023410321241766599_DP ]
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
  !           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
  !           FCENTR - FUNCTION VALUE AT MID POINT
  !           ABSC   - ABSCISSA
  !           FVAL   - FUNCTION VALUE
  !           SAVFUN - ARRAY OF FUNCTION VALUES WHICH HAVE ALREADY BEEN
  !                    COMPUTED
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
  !* FIRST EXECUTABLE STATEMENT  DQNG
  epmach = eps_dp
  uflow = tiny_dp
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  Result = 0._DP
  Abserr = 0._DP
  Neval = 0
  Ier = 6
  IF( Epsabs>0._DP .OR. Epsrel>=MAX(0.5E+02_DP*epmach,0.5E-28_DP) ) THEN
    hlgth = 0.5_DP*(B-A)
    dhlgth = ABS(hlgth)
    centr = 0.5_DP*(B+A)
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
          res10 = 0._DP
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
          reskh = 0.5_DP*res21
          resasc = w21b(6)*ABS(fcentr-reskh)
          DO k = 1, 5
            resasc = resasc + w21a(k)*(ABS(fv1(k)-reskh)+ABS(fv2(k)-reskh))&
              + w21b(k)*(ABS(fv3(k)-reskh)+ABS(fv4(k)-reskh))
          END DO
          Abserr = ABS((res21-res10)*hlgth)
          resasc = resasc*dhlgth
      END SELECT
      IF( resasc/=0._DP .AND. Abserr/=0._DP )&
        Abserr = resasc*MIN(1._DP,(0.2E+03_DP*Abserr/resasc)**1.5_DP)
      IF( resabs>uflow/(0.5E+02_DP*epmach) )&
        Abserr = MAX((epmach*0.5E+02_DP)*resabs,Abserr)
      IF( Abserr<=MAX(Epsabs,Epsrel*ABS(Result)) ) Ier = 0
      !- **JUMP OUT OF DO-LOOP
      IF( Ier==0 ) RETURN
    END DO
  END IF
  ERROR STOP 'DQNG : ABNORMAL RETURN'
  !
  RETURN
END SUBROUTINE DQNG