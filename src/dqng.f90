!*==DQNG.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQNG
SUBROUTINE DQNG(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier)
  IMPLICIT NONE
  !*--DQNG5
  !***BEGIN PROLOGUE  DQNG
  !***PURPOSE  The routine calculates an approximation result to a
  !            given definite integral I = integral of F over (A,B),
  !            hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A1
  !***TYPE      DOUBLE PRECISION (QNG-S, DQNG-D)
  !***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD(PATTERSON) RULES,
  !             NONADAPTIVE, QUADPACK, QUADRATURE, SMOOTH INTEGRAND
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
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
  !                    If  EPSABS.LE.0
  !                    And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DQNG
  !
  REAL(8) :: A, absc, Abserr, B, centr, dhlgth, D1MACH, &
    epmach, Epsabs, Epsrel, F, fcentr, fval, fval1, &
    fval2, fv1, fv2, fv3, fv4, hlgth, Result, res10, &
    res21, res43, res87, resabs, resasc, reskh, &
    savfun, uflow, w10, w21a, w21b, w43a, w43b, w87a, &
    w87b, x1, x2, x3, x4
  INTEGER Ier, ipx, k, l, Neval
  EXTERNAL F
  !
  DIMENSION fv1(5), fv2(5), fv3(5), fv4(5), x1(5), x2(5), x3(11), &
    x4(22), w10(5), w21a(5), w21b(6), w43a(10), w43b(12), &
    w87a(21), w87b(23), savfun(21)
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
  SAVE x1, w10, x2, w21a, w21b, x3, w43a, w43b, x4, w87a, w87b
  DATA x1(1)/0.973906528517171720077964012084452D0/
  DATA x1(2)/0.865063366688984510732096688423493D0/
  DATA x1(3)/0.679409568299024406234327365114874D0/
  DATA x1(4)/0.433395394129247190799265943165784D0/
  DATA x1(5)/0.148874338981631210884826001129720D0/
  DATA w10(1)/0.066671344308688137593568809893332D0/
  DATA w10(2)/0.149451349150580593145776339657697D0/
  DATA w10(3)/0.219086362515982043995534934228163D0/
  DATA w10(4)/0.269266719309996355091226921569469D0/
  DATA w10(5)/0.295524224714752870173892994651338D0/
  !
  DATA x2(1)/0.995657163025808080735527280689003D0/
  DATA x2(2)/0.930157491355708226001207180059508D0/
  DATA x2(3)/0.780817726586416897063717578345042D0/
  DATA x2(4)/0.562757134668604683339000099272694D0/
  DATA x2(5)/0.294392862701460198131126603103866D0/
  DATA w21a(1)/0.032558162307964727478818972459390D0/
  DATA w21a(2)/0.075039674810919952767043140916190D0/
  DATA w21a(3)/0.109387158802297641899210590325805D0/
  DATA w21a(4)/0.134709217311473325928054001771707D0/
  DATA w21a(5)/0.147739104901338491374841515972068D0/
  DATA w21b(1)/0.011694638867371874278064396062192D0/
  DATA w21b(2)/0.054755896574351996031381300244580D0/
  DATA w21b(3)/0.093125454583697605535065465083366D0/
  DATA w21b(4)/0.123491976262065851077958109831074D0/
  DATA w21b(5)/0.142775938577060080797094273138717D0/
  DATA w21b(6)/0.149445554002916905664936468389821D0/
  !
  DATA x3(1)/0.999333360901932081394099323919911D0/
  DATA x3(2)/0.987433402908088869795961478381209D0/
  DATA x3(3)/0.954807934814266299257919200290473D0/
  DATA x3(4)/0.900148695748328293625099494069092D0/
  DATA x3(5)/0.825198314983114150847066732588520D0/
  DATA x3(6)/0.732148388989304982612354848755461D0/
  DATA x3(7)/0.622847970537725238641159120344323D0/
  DATA x3(8)/0.499479574071056499952214885499755D0/
  DATA x3(9)/0.364901661346580768043989548502644D0/
  DATA x3(10)/0.222254919776601296498260928066212D0/
  DATA x3(11)/0.074650617461383322043914435796506D0/
  DATA w43a(1)/0.016296734289666564924281974617663D0/
  DATA w43a(2)/0.037522876120869501461613795898115D0/
  DATA w43a(3)/0.054694902058255442147212685465005D0/
  DATA w43a(4)/0.067355414609478086075553166302174D0/
  DATA w43a(5)/0.073870199632393953432140695251367D0/
  DATA w43a(6)/0.005768556059769796184184327908655D0/
  DATA w43a(7)/0.027371890593248842081276069289151D0/
  DATA w43a(8)/0.046560826910428830743339154433824D0/
  DATA w43a(9)/0.061744995201442564496240336030883D0/
  DATA w43a(10)/0.071387267268693397768559114425516D0/
  DATA w43b(1)/0.001844477640212414100389106552965D0/
  DATA w43b(2)/0.010798689585891651740465406741293D0/
  DATA w43b(3)/0.021895363867795428102523123075149D0/
  DATA w43b(4)/0.032597463975345689443882222526137D0/
  DATA w43b(5)/0.042163137935191811847627924327955D0/
  DATA w43b(6)/0.050741939600184577780189020092084D0/
  DATA w43b(7)/0.058379395542619248375475369330206D0/
  DATA w43b(8)/0.064746404951445885544689259517511D0/
  DATA w43b(9)/0.069566197912356484528633315038405D0/
  DATA w43b(10)/0.072824441471833208150939535192842D0/
  DATA w43b(11)/0.074507751014175118273571813842889D0/
  DATA w43b(12)/0.074722147517403005594425168280423D0/
  !
  DATA x4(1)/0.999902977262729234490529830591582D0/
  DATA x4(2)/0.997989895986678745427496322365960D0/
  DATA x4(3)/0.992175497860687222808523352251425D0/
  DATA x4(4)/0.981358163572712773571916941623894D0/
  DATA x4(5)/0.965057623858384619128284110607926D0/
  DATA x4(6)/0.943167613133670596816416634507426D0/
  DATA x4(7)/0.915806414685507209591826430720050D0/
  DATA x4(8)/0.883221657771316501372117548744163D0/
  DATA x4(9)/0.845710748462415666605902011504855D0/
  DATA x4(10)/0.803557658035230982788739474980964D0/
  DATA x4(11)/0.757005730685495558328942793432020D0/
  DATA x4(12)/0.706273209787321819824094274740840D0/
  DATA x4(13)/0.651589466501177922534422205016736D0/
  DATA x4(14)/0.593223374057961088875273770349144D0/
  DATA x4(15)/0.531493605970831932285268948562671D0/
  DATA x4(16)/0.466763623042022844871966781659270D0/
  DATA x4(17)/0.399424847859218804732101665817923D0/
  DATA x4(18)/0.329874877106188288265053371824597D0/
  DATA x4(19)/0.258503559202161551802280975429025D0/
  DATA x4(20)/0.185695396568346652015917141167606D0/
  DATA x4(21)/0.111842213179907468172398359241362D0/
  DATA x4(22)/0.037352123394619870814998165437704D0/
  DATA w87a(1)/0.008148377384149172900002878448190D0/
  DATA w87a(2)/0.018761438201562822243935059003794D0/
  DATA w87a(3)/0.027347451050052286161582829741283D0/
  DATA w87a(4)/0.033677707311637930046581056957588D0/
  DATA w87a(5)/0.036935099820427907614589586742499D0/
  DATA w87a(6)/0.002884872430211530501334156248695D0/
  DATA w87a(7)/0.013685946022712701888950035273128D0/
  DATA w87a(8)/0.023280413502888311123409291030404D0/
  DATA w87a(9)/0.030872497611713358675466394126442D0/
  DATA w87a(10)/0.035693633639418770719351355457044D0/
  DATA w87a(11)/0.000915283345202241360843392549948D0/
  DATA w87a(12)/0.005399280219300471367738743391053D0/
  DATA w87a(13)/0.010947679601118931134327826856808D0/
  DATA w87a(14)/0.016298731696787335262665703223280D0/
  DATA w87a(15)/0.021081568889203835112433060188190D0/
  DATA w87a(16)/0.025370969769253827243467999831710D0/
  DATA w87a(17)/0.029189697756475752501446154084920D0/
  DATA w87a(18)/0.032373202467202789685788194889595D0/
  DATA w87a(19)/0.034783098950365142750781997949596D0/
  DATA w87a(20)/0.036412220731351787562801163687577D0/
  DATA w87a(21)/0.037253875503047708539592001191226D0/
  DATA w87b(1)/0.000274145563762072350016527092881D0/
  DATA w87b(2)/0.001807124155057942948341311753254D0/
  DATA w87b(3)/0.004096869282759164864458070683480D0/
  DATA w87b(4)/0.006758290051847378699816577897424D0/
  DATA w87b(5)/0.009549957672201646536053581325377D0/
  DATA w87b(6)/0.012329447652244853694626639963780D0/
  DATA w87b(7)/0.015010447346388952376697286041943D0/
  DATA w87b(8)/0.017548967986243191099665352925900D0/
  DATA w87b(9)/0.019938037786440888202278192730714D0/
  DATA w87b(10)/0.022194935961012286796332102959499D0/
  DATA w87b(11)/0.024339147126000805470360647041454D0/
  DATA w87b(12)/0.026374505414839207241503786552615D0/
  DATA w87b(13)/0.028286910788771200659968002987960D0/
  DATA w87b(14)/0.030052581128092695322521110347341D0/
  DATA w87b(15)/0.031646751371439929404586051078883D0/
  DATA w87b(16)/0.033050413419978503290785944862689D0/
  DATA w87b(17)/0.034255099704226061787082821046821D0/
  DATA w87b(18)/0.035262412660156681033782717998428D0/
  DATA w87b(19)/0.036076989622888701185500318003895D0/
  DATA w87b(20)/0.036698604498456094498018047441094D0/
  DATA w87b(21)/0.037120549269832576114119958413599D0/
  DATA w87b(22)/0.037334228751935040321235449094698D0/
  DATA w87b(23)/0.037361073762679023410321241766599D0/
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
  !***FIRST EXECUTABLE STATEMENT  DQNG
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  Result = 0.0D+00
  Abserr = 0.0D+00
  Neval = 0
  Ier = 6
  IF ( Epsabs>0.0D+00.OR.Epsrel>=MAX(0.5D+02*epmach,0.5D-28) ) THEN
    hlgth = 0.5D+00*(B-A)
    dhlgth = ABS(hlgth)
    centr = 0.5D+00*(B+A)
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
          ENDDO
          DO k = 1, 11
            ipx = ipx + 1
            absc = hlgth*x3(k)
            fval = F(absc+centr) + F(centr-absc)
            res43 = res43 + fval*w43b(k)
            savfun(ipx) = fval
          ENDDO
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
          ENDDO
          DO k = 1, 22
            absc = hlgth*x4(k)
            res87 = res87 + w87b(k)*(F(absc+centr)+F(centr-absc))
          ENDDO
          Result = res87*hlgth
          Abserr = ABS((res87-res43)*hlgth)
        CASE DEFAULT
          res10 = 0.0D+00
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
          ENDDO
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
          ENDDO
          !
          !          TEST FOR CONVERGENCE.
          !
          Result = res21*hlgth
          resabs = resabs*dhlgth
          reskh = 0.5D+00*res21
          resasc = w21b(6)*ABS(fcentr-reskh)
          DO k = 1, 5
            resasc = resasc + w21a(k)*(ABS(fv1(k)-reskh)+ABS(fv2(k)-reskh))&
              + w21b(k)*(ABS(fv3(k)-reskh)+ABS(fv4(k)-reskh))
          ENDDO
          Abserr = ABS((res21-res10)*hlgth)
          resasc = resasc*dhlgth
      END SELECT
      IF ( resasc/=0.0D+00.AND.Abserr/=0.0D+00 )&
        Abserr = resasc*MIN(0.1D+01,(0.2D+03*Abserr/resasc)**1.5D+00)
      IF ( resabs>uflow/(0.5D+02*epmach) )&
        Abserr = MAX((epmach*0.5D+02)*resabs,Abserr)
      IF ( Abserr<=MAX(Epsabs,Epsrel*ABS(Result)) ) Ier = 0
      ! ***JUMP OUT OF DO-LOOP
      IF ( Ier==0 ) GOTO 99999
    ENDDO
  ENDIF
  CALL XERMSG('SLATEC','DQNG','ABNORMAL RETURN',Ier,0)
  99999 CONTINUE
  END SUBROUTINE DQNG
