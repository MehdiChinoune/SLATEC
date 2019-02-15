!DECK PSIXN
FUNCTION PSIXN(N)
  IMPLICIT NONE
  REAL PSIXN
  !***BEGIN PROLOGUE  PSIXN
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to EXINT
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PSIXN-S, DPSIXN-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     This subroutine returns values of PSI(X)=derivative of log
  !     GAMMA(X), X .GT. 0.0 at integer arguments. A table look-up is
  !     performed for N .LE. 100, and the asymptotic expansion is
  !     evaluated for N .GT. 100.
  !
  !***SEE ALSO  EXINT
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  PSIXN
  !
  INTEGER N, k
  REAL ax, b, c, fn, rfn2, trm, s, wdtol
  REAL R1MACH
  DIMENSION b(6), c(100)
  !-----------------------------------------------------------------------
  !             PSIXN(N), N = 1,100
  !-----------------------------------------------------------------------
  DATA c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), &
    c(10), c(11), c(12), c(13), c(14), c(15), c(16), c(17), &
    c(18), c(19), c(20), c(21), c(22), c(23), &
    c(24)/ - 5.77215664901532861E-01, 4.22784335098467139E-01, &
    9.22784335098467139E-01, 1.25611766843180047E+00, &
    1.50611766843180047E+00, 1.70611766843180047E+00, &
    1.87278433509846714E+00, 2.01564147795561000E+00, &
    2.14064147795561000E+00, 2.25175258906672111E+00, &
    2.35175258906672111E+00, 2.44266167997581202E+00, &
    2.52599501330914535E+00, 2.60291809023222227E+00, &
    2.67434666166079370E+00, 2.74101332832746037E+00, &
    2.80351332832746037E+00, 2.86233685773922507E+00, &
    2.91789241329478063E+00, 2.97052399224214905E+00, &
    3.02052399224214905E+00, 3.06814303986119667E+00, &
    3.11359758531574212E+00, 3.15707584618530734E+00/
  DATA c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32), &
    c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40), &
    c(41), c(42), c(43), c(44), c(45), c(46), c(47), &
    c(48)/3.19874251285197401E+00, 3.23874251285197401E+00, &
    3.27720405131351247E+00, 3.31424108835054951E+00, &
    3.34995537406483522E+00, 3.38443813268552488E+00, &
    3.41777146601885821E+00, 3.45002953053498724E+00, &
    3.48127953053498724E+00, 3.51158256083801755E+00, &
    3.54099432554389990E+00, 3.56956575411532847E+00, &
    3.59734353189310625E+00, 3.62437055892013327E+00, &
    3.65068634839381748E+00, 3.67632737403484313E+00, &
    3.70132737403484313E+00, 3.72571761793728215E+00, &
    3.74952714174680596E+00, 3.77278295570029433E+00, &
    3.79551022842756706E+00, 3.81773245064978928E+00, &
    3.83947158108457189E+00, 3.86074817682925274E+00/
  DATA c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56), &
    c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64), &
    c(65), c(66), c(67), c(68), c(69), c(70), c(71), &
    c(72)/3.88158151016258607E+00, 3.90198967342789220E+00, &
    3.92198967342789220E+00, 3.94159751656514710E+00, &
    3.96082828579591633E+00, 3.97969621032421822E+00, &
    3.99821472884273674E+00, 4.01639654702455492E+00, &
    4.03425368988169777E+00, 4.05179754953082058E+00, &
    4.06903892884116541E+00, 4.08598808138353829E+00, &
    4.10265474805020496E+00, 4.11904819067315578E+00, &
    4.13517722293122029E+00, 4.15105023880423617E+00, &
    4.16667523880423617E+00, 4.18205985418885155E+00, &
    4.19721136934036670E+00, 4.21213674247469506E+00, &
    4.22684262482763624E+00, 4.24133537845082464E+00, &
    4.25562109273653893E+00, 4.26970559977879245E+00/
  DATA c(73), c(74), c(75), c(76), c(77), c(78), c(79), c(80), &
    c(81), c(82), c(83), c(84), c(85), c(86), c(87), c(88), &
    c(89), c(90), c(91), c(92), c(93), c(94), c(95), &
    c(96)/4.28359448866768134E+00, 4.29729311880466764E+00, &
    4.31080663231818115E+00, 4.32413996565151449E+00, &
    4.33729786038835659E+00, 4.35028487337536958E+00, &
    4.36310538619588240E+00, 4.37576361404398366E+00, &
    4.38826361404398366E+00, 4.40060929305632934E+00, &
    4.41280441500754886E+00, 4.42485260777863319E+00, &
    4.43675736968339510E+00, 4.44852207556574804E+00, &
    4.46014998254249223E+00, 4.47164423541605544E+00, &
    4.48300787177969181E+00, 4.49424382683587158E+00, &
    4.50535493794698269E+00, 4.51634394893599368E+00, &
    4.52721351415338499E+00, 4.53796620232542800E+00, &
    4.54860450019776842E+00, 4.55913081598724211E+00/
  DATA c(97), c(98), c(99), c(100)/4.56954748265390877E+00, &
    4.57985676100442424E+00, 4.59006084263707730E+00, &
    4.60016185273808740E+00/
  !-----------------------------------------------------------------------
  !             COEFFICIENTS OF ASYMPTOTIC EXPANSION
  !-----------------------------------------------------------------------
  DATA b(1), b(2), b(3), b(4), b(5), b(6)/8.33333333333333333E-02, &
    -8.33333333333333333E-03, 3.96825396825396825E-03, &
    -4.16666666666666666E-03, 7.57575757575757576E-03, &
    -2.10927960927960928E-02/
  !
  !***FIRST EXECUTABLE STATEMENT  PSIXN
  IF ( N>100 ) THEN
    wdtol = MAX(R1MACH(4),1.0E-18)
    fn = N
    ax = 1.0E0
    s = -0.5E0/fn
    IF ( ABS(s)>wdtol ) THEN
      rfn2 = 1.0E0/(fn*fn)
      DO k = 1, 6
        ax = ax*rfn2
        trm = -b(k)*ax
        IF ( ABS(trm)<wdtol ) EXIT
        s = s + trm
      ENDDO
    ENDIF
    PSIXN = s + LOG(fn)
  ELSE
    PSIXN = c(N)
    RETURN
  ENDIF
END FUNCTION PSIXN
