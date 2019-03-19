!** DBIE
REAL(8) FUNCTION DBIE(X)
  IMPLICIT NONE
  !>
  !***
  !  Calculate the Bairy function for a negative argument and an
  !            exponentially scaled Bairy function for a non-negative
  !            argument.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      DOUBLE PRECISION (BIE-S, DBIE-D)
  !***
  ! **Keywords:**  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBIE(X) calculates the double precision Airy function of the
  ! second kind or the double precision exponentially scaled Airy
  ! function of the second kind, depending on the value of the
  ! double precision argument X.
  !
  ! Evaluate BI(X) for X .LE. 0.0  and  BI(X)*EXP(-ZETA)  where
  ! ZETA = 2/3 * X**(3/2)  for X .GE. 0.0
  !
  !
  ! Series for BIF        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   1.45E-32
  !                                         log weighted error  31.84
  !                               significant figures required  30.85
  !                                    decimal places required  32.40
  !
  !
  ! Series for BIG        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   1.29E-33
  !                                         log weighted error  32.89
  !                               significant figures required  31.48
  !                                    decimal places required  33.45
  !
  !
  ! Series for BIF2       on the interval  1.00000E+00 to  8.00000E+00
  !                                        with weighted error   6.08E-32
  !                                         log weighted error  31.22
  !                        approx significant figures required  30.8
  !                                    decimal places required  31.80
  !
  !
  ! Series for BIG2       on the interval  1.00000E+00 to  8.00000E+00
  !                                        with weighted error   4.91E-33
  !                                         log weighted error  32.31
  !                        approx significant figures required  31.6
  !                                    decimal places required  32.90
  !
  !
  ! Series for BIP1       on the interval  1.25000E-01 to  3.53553E-01
  !                                        with weighted error   1.06E-32
  !                                         log weighted error  31.98
  !                               significant figures required  30.61
  !                                    decimal places required  32.81
  !
  !
  ! Series for BIP2       on the interval  0.          to  1.25000E-01
  !                                        with weighted error   4.04E-33
  !                                         log weighted error  32.39
  !                               significant figures required  31.15
  !                                    decimal places required  33.37
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9AIMP, DCSEVL, INITDS

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL eta
  INTEGER INITDS, nbif, nbif2, nbig, nbig2, nbip1, nbip2
  REAL(8) :: X, bifcs(13), bigcs(13), bif2cs(15), big2cs(15), &
    bip1cs(47), bip2cs(88), atr, btr, sqrtx, theta, &
    xbig, xm, x3sml, x32sml, z, D1MACH, DCSEVL
  LOGICAL first
  SAVE bifcs, bigcs, bif2cs, big2cs, bip1cs, bip2cs, atr, btr, &
    nbif, nbig, nbif2, nbig2, nbip1, nbip2, x3sml, x32sml, xbig, &
    first
  DATA bifcs(1)/ - .16730216471986649483537423928176D-1/
  DATA bifcs(2)/ + .10252335834249445611426362777757D+0/
  DATA bifcs(3)/ + .17083092507381516539429650242013D-2/
  DATA bifcs(4)/ + .11862545467744681179216459210040D-4/
  DATA bifcs(5)/ + .44932907017792133694531887927242D-7/
  DATA bifcs(6)/ + .10698207143387889067567767663628D-9/
  DATA bifcs(7)/ + .17480643399771824706010517628573D-12/
  DATA bifcs(8)/ + .20810231071761711025881891834399D-15/
  DATA bifcs(9)/ + .18849814695665416509927971733333D-18/
  DATA bifcs(10)/ + .13425779173097804625882666666666D-21/
  DATA bifcs(11)/ + .77159593429658887893333333333333D-25/
  DATA bifcs(12)/ + .36533879617478566399999999999999D-28/
  DATA bifcs(13)/ + .14497565927953066666666666666666D-31/
  DATA bigcs(1)/ + .22466223248574522283468220139024D-1/
  DATA bigcs(2)/ + .37364775453019545441727561666752D-1/
  DATA bigcs(3)/ + .44476218957212285696215294326639D-3/
  DATA bigcs(4)/ + .24708075636329384245494591948882D-5/
  DATA bigcs(5)/ + .79191353395149635134862426285596D-8/
  DATA bigcs(6)/ + .16498079851827779880887872402706D-10/
  DATA bigcs(7)/ + .24119906664835455909247501122841D-13/
  DATA bigcs(8)/ + .26103736236091436985184781269333D-16/
  DATA bigcs(9)/ + .21753082977160323853123792000000D-19/
  DATA bigcs(10)/ + .14386946400390433219483733333333D-22/
  DATA bigcs(11)/ + .77349125612083468629333333333333D-26/
  DATA bigcs(12)/ + .34469292033849002666666666666666D-29/
  DATA bigcs(13)/ + .12938919273216000000000000000000D-32/
  DATA bif2cs(1)/ + .0998457269381604104468284257993D+0/
  DATA bif2cs(2)/ + .47862497786300553772211467318231D+0/
  DATA bif2cs(3)/ + .25155211960433011771324415436675D-1/
  DATA bif2cs(4)/ + .58206938852326456396515697872216D-3/
  DATA bif2cs(5)/ + .74997659644377865943861457378217D-5/
  DATA bif2cs(6)/ + .61346028703493836681403010356474D-7/
  DATA bif2cs(7)/ + .34627538851480632900434268733359D-9/
  DATA bif2cs(8)/ + .14288910080270254287770846748931D-11/
  DATA bif2cs(9)/ + .44962704298334641895056472179200D-14/
  DATA bif2cs(10)/ + .11142323065833011708428300106666D-16/
  DATA bif2cs(11)/ + .22304791066175002081517866666666D-19/
  DATA bif2cs(12)/ + .36815778736393142842922666666666D-22/
  DATA bif2cs(13)/ + .50960868449338261333333333333333D-25/
  DATA bif2cs(14)/ + .60003386926288554666666666666666D-28/
  DATA bif2cs(15)/ + .60827497446570666666666666666666D-31/
  DATA big2cs(1)/ + .033305662145514340465176188111647D+0/
  DATA big2cs(2)/ + .161309215123197067613287532084943D+0/
  DATA big2cs(3)/ + .631900730961342869121615634921173D-2/
  DATA big2cs(4)/ + .118790456816251736389780192304567D-3/
  DATA big2cs(5)/ + .130453458862002656147116485012843D-5/
  DATA big2cs(6)/ + .937412599553521729546809615508936D-8/
  DATA big2cs(7)/ + .474580188674725153788510169834595D-10/
  DATA big2cs(8)/ + .178310726509481399800065667560946D-12/
  DATA big2cs(9)/ + .516759192784958180374276356640000D-15/
  DATA big2cs(10)/ + .119004508386827125129496251733333D-17/
  DATA big2cs(11)/ + .222982880666403517277063466666666D-20/
  DATA big2cs(12)/ + .346551923027689419722666666666666D-23/
  DATA big2cs(13)/ + .453926336320504514133333333333333D-26/
  DATA big2cs(14)/ + .507884996513522346666666666666666D-29/
  DATA big2cs(15)/ + .491020674696533333333333333333333D-32/
  DATA bip1cs(1)/ - .83220474779434474687471864707973D-1/
  DATA bip1cs(2)/ + .11461189273711742889920226128031D-1/
  DATA bip1cs(3)/ + .42896440718911509494134472566635D-3/
  DATA bip1cs(4)/ - .14906639379950514017847677732954D-3/
  DATA bip1cs(5)/ - .13076597267876290663136340998881D-4/
  DATA bip1cs(6)/ + .63275983961030344754535716032494D-5/
  DATA bip1cs(7)/ - .42226696982681924884778515889433D-6/
  DATA bip1cs(8)/ - .19147186298654689632835494181277D-6/
  DATA bip1cs(9)/ + .64531062845583173611038157880934D-7/
  DATA bip1cs(10)/ - .78448546771397719289748310448628D-8/
  DATA bip1cs(11)/ - .96077216623785085879198533565432D-9/
  DATA bip1cs(12)/ + .70004713316443966339006074402068D-9/
  DATA bip1cs(13)/ - .17731789132814932022083128056698D-9/
  DATA bip1cs(14)/ + .22720894783465236347282126389311D-10/
  DATA bip1cs(15)/ + .16540456313972049847032860681891D-11/
  DATA bip1cs(16)/ - .18517125559292316390755369896693D-11/
  DATA bip1cs(17)/ + .59576312477117290165680715534277D-12/
  DATA bip1cs(18)/ - .12194348147346564781055769498986D-12/
  DATA bip1cs(19)/ + .13347869253513048815386347813597D-13/
  DATA bip1cs(20)/ + .17278311524339746664384792889731D-14/
  DATA bip1cs(21)/ - .14590732013016720735268871713166D-14/
  DATA bip1cs(22)/ + .49010319927115819978994989520104D-15/
  DATA bip1cs(23)/ - .11556545519261548129262972762521D-15/
  DATA bip1cs(24)/ + .19098807367072411430671732441524D-16/
  DATA bip1cs(25)/ - .11768966854492179886913995957862D-17/
  DATA bip1cs(26)/ - .63271925149530064474537459677047D-18/
  DATA bip1cs(27)/ + .33861838880715361614130191322316D-18/
  DATA bip1cs(28)/ - .10725825321758625254992162219622D-18/
  DATA bip1cs(29)/ + .25995709605617169284786933115562D-19/
  DATA bip1cs(30)/ - .48477583571081193660962309494101D-20/
  DATA bip1cs(31)/ + .55298913982121625361505513198933D-21/
  DATA bip1cs(32)/ + .49421660826069471371748197444266D-22/
  DATA bip1cs(33)/ - .55162121924145707458069720814933D-22/
  DATA bip1cs(34)/ + .21437560417632550086631884499626D-22/
  DATA bip1cs(35)/ - .61910313387655605798785061137066D-23/
  DATA bip1cs(36)/ + .14629362707391245659830967336959D-23/
  DATA bip1cs(37)/ - .27918484471059005576177866069333D-24/
  DATA bip1cs(38)/ + .36455703168570246150906795349333D-25/
  DATA bip1cs(39)/ + .58511821906188711839382459733333D-27/
  DATA bip1cs(40)/ - .24946950487566510969745047551999D-26/
  DATA bip1cs(41)/ + .10979323980338380977919579477333D-26/
  DATA bip1cs(42)/ - .34743388345961115015034088106666D-27/
  DATA bip1cs(43)/ + .91373402635349697363171082240000D-28/
  DATA bip1cs(44)/ - .20510352728210629186247720959999D-28/
  DATA bip1cs(45)/ + .37976985698546461748651622399999D-29/
  DATA bip1cs(46)/ - .48479458497755565887848448000000D-30/
  DATA bip1cs(47)/ - .10558306941230714314205866666666D-31/
  DATA bip2cs(1)/ - .11359673758598867913797310895527D+0/
  DATA bip2cs(2)/ + .41381473947881595760052081171444D-2/
  DATA bip2cs(3)/ + .13534706221193329857696921727508D-3/
  DATA bip2cs(4)/ + .10427316653015353405887183456780D-4/
  DATA bip2cs(5)/ + .13474954767849907889589911958925D-5/
  DATA bip2cs(6)/ + .16965374054383983356062511163756D-6/
  DATA bip2cs(7)/ - .10096500865641624301366228396373D-7/
  DATA bip2cs(8)/ - .16729119493778475127836973095943D-7/
  DATA bip2cs(9)/ - .45815364485068383217152795613391D-8/
  DATA bip2cs(10)/ + .37366813665655477274064749384284D-9/
  DATA bip2cs(11)/ + .57669303201452448119584643502111D-9/
  DATA bip2cs(12)/ + .62181265087850324095393408792371D-10/
  DATA bip2cs(13)/ - .63294120282743068241589177281354D-10/
  DATA bip2cs(14)/ - .14915047908598767633999091989487D-10/
  DATA bip2cs(15)/ + .78896213942486771938172394294891D-11/
  DATA bip2cs(16)/ + .24960513721857797984888064000127D-11/
  DATA bip2cs(17)/ - .12130075287291659477746664734814D-11/
  DATA bip2cs(18)/ - .37404939108727277887343460402716D-12/
  DATA bip2cs(19)/ + .22377278140321476798783446931091D-12/
  DATA bip2cs(20)/ + .47490296312192466341986077472514D-13/
  DATA bip2cs(21)/ - .45261607991821224810605655831294D-13/
  DATA bip2cs(22)/ - .30172271841986072645112245876020D-14/
  DATA bip2cs(23)/ + .91058603558754058327592683478908D-14/
  DATA bip2cs(24)/ - .98149238033807062926643864207709D-15/
  DATA bip2cs(25)/ - .16429400647889465253601245251589D-14/
  DATA bip2cs(26)/ + .55334834214274215451182114635164D-15/
  DATA bip2cs(27)/ + .21750479864482655984374381998156D-15/
  DATA bip2cs(28)/ - .17379236200220656971287029558087D-15/
  DATA bip2cs(29)/ - .10470023471443714959283909313604D-17/
  DATA bip2cs(30)/ + .39219145986056386925441403311462D-16/
  DATA bip2cs(31)/ - .11621293686345196925824005665910D-16/
  DATA bip2cs(32)/ - .54027474491754245533735411307773D-17/
  DATA bip2cs(33)/ + .45441582123884610882675428553304D-17/
  DATA bip2cs(34)/ - .28775599625221075729427585480086D-18/
  DATA bip2cs(35)/ - .10017340927225341243596162960440D-17/
  DATA bip2cs(36)/ + .44823931215068369856332561906313D-18/
  DATA bip2cs(37)/ + .76135968654908942328948982366775D-19/
  DATA bip2cs(38)/ - .14448324094881347238956060145422D-18/
  DATA bip2cs(39)/ + .40460859449205362251624847392112D-19/
  DATA bip2cs(40)/ + .20321085700338446891325190707277D-19/
  DATA bip2cs(41)/ - .19602795471446798718272758041962D-19/
  DATA bip2cs(42)/ + .34273038443944824263518958211738D-20/
  DATA bip2cs(43)/ + .37023705853905135480024651593154D-20/
  DATA bip2cs(44)/ - .26879595172041591131400332966712D-20/
  DATA bip2cs(45)/ + .28121678463531712209714454683364D-21/
  DATA bip2cs(46)/ + .60933963636177797173271119680329D-21/
  DATA bip2cs(47)/ - .38666621897150844994172977893413D-21/
  DATA bip2cs(48)/ + .25989331253566943450895651927228D-22/
  DATA bip2cs(49)/ + .97194393622938503767281175216084D-22/
  DATA bip2cs(50)/ - .59392817834375098415630478204591D-22/
  DATA bip2cs(51)/ + .38864949977113015409591960439444D-23/
  DATA bip2cs(52)/ + .15334307393617272869721512868769D-22/
  DATA bip2cs(53)/ - .97513555209762624036336521409724D-23/
  DATA bip2cs(54)/ + .96340644440489471424741339383726D-24/
  DATA bip2cs(55)/ + .23841999400208880109946748792454D-23/
  DATA bip2cs(56)/ - .16896986315019706184848044205207D-23/
  DATA bip2cs(57)/ + .27352715888928361222578444801478D-24/
  DATA bip2cs(58)/ + .35660016185409578960111685025730D-24/
  DATA bip2cs(59)/ - .30234026608258827249534280666954D-24/
  DATA bip2cs(60)/ + .75002041605973930653144204823232D-25/
  DATA bip2cs(61)/ + .48403287575851388827455319838748D-25/
  DATA bip2cs(62)/ - .54364137654447888432698010297766D-25/
  DATA bip2cs(63)/ + .19281214470820962653345978809756D-25/
  DATA bip2cs(64)/ + .50116355020532656659611814172172D-26/
  DATA bip2cs(65)/ - .95040744582693253786034620869972D-26/
  DATA bip2cs(66)/ + .46372646157101975948696332245611D-26/
  DATA bip2cs(67)/ + .21177170704466954163768170577046D-28/
  DATA bip2cs(68)/ - .15404850268168594303692204548726D-26/
  DATA bip2cs(69)/ + .10387944293201213662047889194441D-26/
  DATA bip2cs(70)/ - .19890078156915416751316728235153D-27/
  DATA bip2cs(71)/ - .21022173878658495471177044522532D-27/
  DATA bip2cs(72)/ + .21353099724525793150633356670491D-27/
  DATA bip2cs(73)/ - .79040810747961342319023537632627D-28/
  DATA bip2cs(74)/ - .16575359960435585049973741763592D-28/
  DATA bip2cs(75)/ + .38868342850124112587625586496537D-28/
  DATA bip2cs(76)/ - .22309237330896866182621562424717D-28/
  DATA bip2cs(77)/ + .27777244420176260265625977404382D-29/
  DATA bip2cs(78)/ + .57078543472657725368712433782772D-29/
  DATA bip2cs(79)/ - .51743084445303852800173371555280D-29/
  DATA bip2cs(80)/ + .18413280751095837198450927071569D-29/
  DATA bip2cs(81)/ + .44422562390957094598544071068647D-30/
  DATA bip2cs(82)/ - .98504142639629801547464958226943D-30/
  DATA bip2cs(83)/ + .58857201353585104884754198881995D-30/
  DATA bip2cs(84)/ - .97636075440429787961402312628595D-31/
  DATA bip2cs(85)/ - .13581011996074695047063597884122D-30/
  DATA bip2cs(86)/ + .13999743518492413270568048380345D-30/
  DATA bip2cs(87)/ - .59754904545248477620884562981118D-31/
  DATA bip2cs(88)/ - .40391653875428313641045327529856D-32/
  DATA atr/8.75069057084843450880771988210148D0/
  DATA btr/ - 2.09383632135605431360096498526268D0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  DBIE
  IF ( first ) THEN
    eta = 0.1*REAL(D1MACH(3))
    nbif = INITDS(bifcs,13,eta)
    nbig = INITDS(bigcs,13,eta)
    nbif2 = INITDS(bif2cs,15,eta)
    nbig2 = INITDS(big2cs,15,eta)
    nbip1 = INITDS(bip1cs,47,eta)
    nbip2 = INITDS(bip2cs,88,eta)
    !
    x3sml = eta**0.3333
    x32sml = 1.3104D0*x3sml**2
    xbig = D1MACH(2)**0.6666D0
  ENDIF
  first = .FALSE.
  !
  IF ( X<(-1.0D0) ) THEN
    CALL D9AIMP(X,xm,theta)
    DBIE = xm*SIN(theta)
    RETURN
    !
  ELSEIF ( X<=1.0D0 ) THEN
    z = 0.D0
    IF ( ABS(X)>x3sml ) z = X**3
    DBIE = 0.625D0 + DCSEVL(z,bifcs,nbif)&
      + X*(0.4375D0+DCSEVL(z,bigcs,nbig))
    IF ( X>x32sml ) DBIE = DBIE*EXP(-2.0D0*X*SQRT(X)/3.0D0)
    RETURN
    !
  ELSEIF ( X<=2.0D0 ) THEN
    z = (2.0D0*X**3-9.0D0)/7.0D0
    DBIE = EXP(-2.0D0*X*SQRT(X)/3.0D0)&
      *(1.125D0+DCSEVL(z,bif2cs,nbif2)+X*(0.625D0+&
      DCSEVL(z,big2cs,nbig2)))
    RETURN
    !
  ELSEIF ( X>4.0D0 ) THEN
    !
    sqrtx = SQRT(X)
    z = -1.0D0
    IF ( X<xbig ) z = 16.D0/(X*sqrtx) - 1.0D0
    DBIE = (0.625D0+DCSEVL(z,bip2cs,nbip2))/SQRT(sqrtx)
    RETURN
  ENDIF
  sqrtx = SQRT(X)
  z = atr/(X*sqrtx) + btr
  DBIE = (0.625D0+DCSEVL(z,bip1cs,nbip1))/SQRT(sqrtx)
  RETURN
END FUNCTION DBIE
