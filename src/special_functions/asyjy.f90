!** ASYJY
PURE SUBROUTINE ASYJY(FUNJY,X,Fnu,Flgjy,In,Y,Wk,Iflw)
  !> Subsidiary to BESJ and BESY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (ASYJY-S, DASYJY-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !      ASYJY computes Bessel functions J and Y
  !      for arguments X>0.0 and orders FNU>=35.0
  !      on FLGJY = 1 and FLGJY = -1 respectively
  !
  !                                  INPUT
  !
  !      FUNJY - external function JAIRY or YAIRY
  !          X - argument, X>0.0E0
  !        FNU - order of the first Bessel function
  !      FLGJY - selection flag
  !              FLGJY =  1.0E0 gives the J function
  !              FLGJY = -1.0E0 gives the Y function
  !         IN - number of functions desired, IN = 1 or 2
  !
  !                                  OUTPUT
  !
  !         Y  - a vector whose first in components contain the sequence
  !       IFLW - a flag indicating underflow or overflow
  !                    return variables for BESJ only
  !      WK(1) = 1 - (X/FNU)**2 = W**2
  !      WK(2) = SQRT(ABS(WK(1)))
  !      WK(3) = ABS(WK(2) - ATAN(WK(2)))  or
  !              ABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
  !            = ABS((2/3)*ZETA**(3/2))
  !      WK(4) = FNU*WK(3)
  !      WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3)
  !      WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3)
  !      WK(7) = FNU**(1/3)
  !
  !     Abstract
  !         ASYJY implements the uniform asymptotic expansion of
  !         the J and Y Bessel functions for FNU>=35 and real
  !         X>0.0E0. The forms are identical except for a change
  !         in sign of some of the terms. This change in sign is
  !         accomplished by means of the flag FLGJY = 1 or -1. On
  !         FLGJY = 1 the AIRY functions AI(X) and DAI(X) are
  !         supplied by the external function JAIRY, and on
  !         FLGJY = -1 the AIRY functions BI(X) and DBI(X) are
  !         supplied by the external function YAIRY.
  !
  !***
  ! **See also:**  BESJ, BESY
  !***
  ! **Routines called:**  I1MACH, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR section.  (WRB)
  USE service, ONLY : digits_sp, min_exp_sp, eps_2_sp, log10_radix_sp, tiny_sp
  !
  INTERFACE
    PURE SUBROUTINE FUNJY(A,B,C,D,E)
      IMPORT SP
      REAL(SP), INTENT(IN) :: A
      REAL(SP), INTENT(INOUT) :: B, C
      REAL(SP), INTENT(OUT) :: D, E
    END SUBROUTINE FUNJY
  END INTERFACE
  INTEGER, INTENT(IN) :: In
  INTEGER, INTENT(OUT) :: Iflw
  REAL(SP), INTENT(IN) :: Flgjy, Fnu, X
  REAL(SP), INTENT(OUT) :: Wk(7), Y(In)
  INTEGER :: i, j, jn, jr, ju, k, kb, klast, kmax(5), kp1, &
    ks, ksp1, kstemp, l, lr, lrp1, iseta, isetb
  REAL(SP) :: abw2, akm, ap, asum, az, bsum, cr(10), crz32, dfi, elim, dr(10), fi, &
    fn, fn2, phi, rcz, rden, relb, rfn2, rtz, rzden, sa, sb, suma, &
    sumb, s1, ta, tau, tb, tfn, tol, t2, upol(10), xx, z, z32
  REAL(SP), PARAMETER :: tols = -6.90775527898214_SP
  REAL(SP), PARAMETER :: con1 = 6.66666666666667E-01_SP, con2 = 3.33333333333333E-01_SP, &
    con548 = 1.04166666666667E-01_SP
  REAL(SP), PARAMETER :: ar(8) = [ 8.35503472222222E-02_SP, 1.28226574556327E-01_SP, &
    2.91849026464140E-01_SP, 8.81627267443758E-01_SP, 3.32140828186277E+00_SP, &
    1.49957629868626E+01_SP, 7.89230130115865E+01_SP, 4.74451538868264E+02_SP ]
  REAL(SP), PARAMETER :: br(10) = [ -1.45833333333333E-01_SP, -9.87413194444444E-02_SP, &
    -1.43312053915895E-01_SP, -3.17227202678414E-01_SP, -9.42429147957120E-01_SP, &
    -3.51120304082635E+00_SP, -1.57272636203680E+01_SP, -8.22814390971859E+01_SP, &
    -4.92355370523671E+02_SP, -3.31621856854797E+03_SP ]
  REAL(SP), PARAMETER :: c(65) = [ -2.08333333333333E-01_SP, 1.25000000000000E-01_SP, &
    3.34201388888889E-01_SP, -4.01041666666667E-01_SP, 7.03125000000000E-02_SP, &
    -1.02581259645062E+00_SP, 1.84646267361111E+00_SP, -8.91210937500000E-01_SP, &
    7.32421875000000E-02_SP, 4.66958442342625E+00_SP, -1.12070026162230E+01_SP, &
    8.78912353515625E+00_SP, -2.36408691406250E+00_SP, 1.12152099609375E-01_SP, &
    -2.82120725582002E+01_SP, 8.46362176746007E+01_SP, -9.18182415432400E+01_SP, &
    4.25349987453885E+01_SP, -7.36879435947963E+00_SP, 2.27108001708984E-01_SP, &
    2.12570130039217E+02_SP, -7.65252468141182E+02_SP, 1.05999045252800E+03_SP, &
    -6.99579627376133E+02_SP, 2.18190511744212E+02_SP, -2.64914304869516E+01_SP, &
    5.72501420974731E-01_SP, -1.91945766231841E+03_SP, 8.06172218173731E+03_SP, &
    -1.35865500064341E+04_SP, 1.16553933368645E+04_SP, -5.30564697861340E+03_SP, &
    1.20090291321635E+03_SP, -1.08090919788395E+02_SP, 1.72772750258446E+00_SP, &
    2.02042913309661E+04_SP, -9.69805983886375E+04_SP, 1.92547001232532E+05_SP, &
    -2.03400177280416E+05_SP, 1.22200464983017E+05_SP, -4.11926549688976E+04_SP, &
    7.10951430248936E+03_SP, -4.93915304773088E+02_SP, 6.07404200127348E+00_SP, &
    -2.42919187900551E+05_SP, 1.31176361466298E+06_SP, -2.99801591853811E+06_SP, &
    3.76327129765640E+06_SP, -2.81356322658653E+06_SP, 1.26836527332162E+06_SP, &
    -3.31645172484564E+05_SP, 4.52187689813627E+04_SP, -2.49983048181121E+03_SP, &
    2.43805296995561E+01_SP, 3.28446985307204E+06_SP, -1.97068191184322E+07_SP, &
    5.09526024926646E+07_SP, -7.41051482115327E+07_SP, 6.63445122747290E+07_SP, &
    -3.75671766607634E+07_SP, 1.32887671664218E+07_SP, -2.78561812808645E+06_SP, &
    3.08186404612662E+05_SP, -1.38860897537170E+04_SP, 1.10017140269247E+02_SP ]
  REAL(SP), PARAMETER :: alfa1(26,2) = RESHAPE( [ -4.44444444444444E-03_SP, &
    -9.22077922077922E-04_SP, -8.84892884892885E-05_SP, 1.65927687832450E-04_SP, &
    2.46691372741793E-04_SP, 2.65995589346255E-04_SP, 2.61824297061501E-04_SP, &
    2.48730437344656E-04_SP, 2.32721040083232E-04_SP, 2.16362485712365E-04_SP, &
    2.00738858762752E-04_SP, 1.86267636637545E-04_SP, 1.73060775917876E-04_SP, &
    1.61091705929016E-04_SP, 1.50274774160908E-04_SP, 1.40503497391270E-04_SP, &
    1.31668816545923E-04_SP, 1.23667445598253E-04_SP, 1.16405271474738E-04_SP, &
    1.09798298372713E-04_SP, 1.03772410422993E-04_SP, 9.82626078369363E-05_SP, &
    9.32120517249503E-05_SP, 8.85710852478712E-05_SP, 8.42963105715700E-05_SP, &
    8.03497548407791E-05_SP, 6.93735541354589E-04_SP, 2.32241745182922E-04_SP, &
    -1.41986273556691E-05_SP, -1.16444931672049E-04_SP, -1.50803558053049E-04_SP, &
    -1.55121924918096E-04_SP, -1.46809756646466E-04_SP, -1.33815503867491E-04_SP, &
    -1.19744975684254E-04_SP, -1.06184319207974E-04_SP, -9.37699549891194E-05_SP, &
    -8.26923045588193E-05_SP, -7.29374348155221E-05_SP, -6.44042357721016E-05_SP, &
    -5.69611566009369E-05_SP, -5.04731044303562E-05_SP, -4.48134868008883E-05_SP, &
    -3.98688727717599E-05_SP, -3.55400532972042E-05_SP, -3.17414256609022E-05_SP, &
    -2.83996793904175E-05_SP, -2.54522720634871E-05_SP, -2.28459297164725E-05_SP, &
    -2.05352753106481E-05_SP, -1.84816217627666E-05_SP, -1.66519330021394E-05_SP ], [26,2] )
  REAL(SP), PARAMETER :: alfa2(26,2) = RESHAPE( [ -3.54211971457744E-04_SP, &
    -1.56161263945159E-04_SP, 3.04465503594936E-05_SP, 1.30198655773243E-04_SP, &
    1.67471106699712E-04_SP, 1.70222587683593E-04_SP, 1.56501427608595E-04_SP, &
    1.36339170977445E-04_SP, 1.14886692029825E-04_SP, 9.45869093034688E-05_SP, &
    7.64498419250898E-05_SP, 6.07570334965197E-05_SP, 4.74394299290509E-05_SP, &
    3.62757512005344E-05_SP, 2.69939714979225E-05_SP, 1.93210938247939E-05_SP, &
    1.30056674793963E-05_SP, 7.82620866744497E-06_SP, 3.59257485819352E-06_SP, &
    1.44040049814252E-07_SP, -2.65396769697939E-06_SP, -4.91346867098486E-06_SP, &
    -6.72739296091248E-06_SP, -8.17269379678658E-06_SP, -9.31304715093561E-06_SP, &
    -1.02011418798016E-05_SP, 3.78194199201773E-04_SP, 2.02471952761816E-04_SP, &
    -6.37938506318862E-05_SP, -2.38598230603006E-04_SP, -3.10916256027362E-04_SP, &
    -3.13680115247576E-04_SP, -2.78950273791323E-04_SP, -2.28564082619141E-04_SP, &
    -1.75245280340847E-04_SP, -1.25544063060690E-04_SP, -8.22982872820208E-05_SP, &
    -4.62860730588116E-05_SP, -1.72334302366962E-05_SP, 5.60690482304602E-06_SP, &
    2.31395443148287E-05_SP, 3.62642745856794E-05_SP, 4.58006124490189E-05_SP, &
    5.24595294959114E-05_SP, 5.68396208545815E-05_SP, 5.94349820393104E-05_SP, &
    6.06478527578422E-05_SP, 6.08023907788436E-05_SP, 6.01577894539460E-05_SP, &
    5.89199657344698E-05_SP, 5.72515823777593E-05_SP, 5.52804375585853E-05_SP ], [26,2] )
  REAL(SP), PARAMETER :: beta1(26,2) = RESHAPE( [ 1.79988721413553E-02_SP, &
    5.59964911064388E-03_SP, 2.88501402231133E-03_SP, 1.80096606761054E-03_SP, &
    1.24753110589199E-03_SP, 9.22878876572938E-04_SP, 7.14430421727287E-04_SP, &
    5.71787281789705E-04_SP, 4.69431007606482E-04_SP, 3.93232835462917E-04_SP, &
    3.34818889318298E-04_SP, 2.88952148495752E-04_SP, 2.52211615549573E-04_SP, &
    2.22280580798883E-04_SP, 1.97541838033063E-04_SP, 1.76836855019718E-04_SP, &
    1.59316899661821E-04_SP, 1.44347930197334E-04_SP, 1.31448068119965E-04_SP, &
    1.20245444949303E-04_SP, 1.10449144504599E-04_SP, 1.01828770740567E-04_SP, &
    9.41998224204238E-05_SP, 8.74130545753834E-05_SP, 8.13466262162801E-05_SP, &
    7.59002269646219E-05_SP, -1.49282953213429E-03_SP, -8.78204709546389E-04_SP, &
    -5.02916549572035E-04_SP, -2.94822138512746E-04_SP, -1.75463996970783E-04_SP, &
    -1.04008550460816E-04_SP, -5.96141953046458E-05_SP, -3.12038929076098E-05_SP, &
    -1.26089735980230E-05_SP, -2.42892608575730E-07_SP, 8.05996165414274E-06_SP, &
    1.36507009262147E-05_SP, 1.73964125472926E-05_SP, 1.98672978842134E-05_SP, &
    2.14463263790823E-05_SP, 2.23954659232457E-05_SP, 2.28967783814713E-05_SP, &
    2.30785389811178E-05_SP, 2.30321976080909E-05_SP, 2.28236073720349E-05_SP, &
    2.25005881105292E-05_SP, 2.20981015361991E-05_SP, 2.16418427448104E-05_SP, &
    2.11507649256221E-05_SP, 2.06388749782171E-05_SP, 2.01165241997082E-05_SP ], [26,2] )
  REAL(SP), PARAMETER :: beta2(26,2) = RESHAPE( [ 5.52213076721293E-04_SP, &
    4.47932581552385E-04_SP, 2.79520653992021E-04_SP, 1.52468156198447E-04_SP, &
    6.93271105657044E-05_SP, 1.76258683069991E-05_SP, -1.35744996343269E-05_SP, &
    -3.17972413350427E-05_SP, -4.18861861696693E-05_SP, -4.69004889379141E-05_SP, &
    -4.87665447413787E-05_SP, -4.87010031186735E-05_SP, -4.74755620890087E-05_SP, &
    -4.55813058138628E-05_SP, -4.33309644511266E-05_SP, -4.09230193157750E-05_SP, &
    -3.84822638603221E-05_SP, -3.60857167535411E-05_SP, -3.37793306123367E-05_SP, &
    -3.15888560772110E-05_SP, -2.95269561750807E-05_SP, -2.75978914828336E-05_SP, &
    -2.58006174666884E-05_SP, -2.41308356761280E-05_SP, -2.25823509518346E-05_SP, &
    -2.11479656768913E-05_SP, -4.74617796559960E-04_SP, -4.77864567147321E-04_SP, &
    -3.20390228067038E-04_SP, -1.61105016119962E-04_SP, -4.25778101285435E-05_SP, &
    3.44571294294968E-05_SP, 7.97092684075675E-05_SP, 1.03138236708272E-04_SP, &
    1.12466775262204E-04_SP, 1.13103642108481E-04_SP, 1.08651634848774E-04_SP, &
    1.01437951597662E-04_SP, 9.29298396593364E-05_SP, 8.40293133016090E-05_SP, &
    7.52727991349134E-05_SP, 6.69632521975731E-05_SP, 5.92564547323195E-05_SP, &
    5.22169308826976E-05_SP, 4.58539485165361E-05_SP, 4.01445513891487E-05_SP, &
    3.50481730031328E-05_SP, 3.05157995034347E-05_SP, 2.64956119950516E-05_SP, &
    2.29363633690998E-05_SP, 1.97893056664022E-05_SP, 1.70091984636413E-05_SP ], [26,2] )
  REAL(SP), PARAMETER :: beta3(26,1) = RESHAPE( [ 7.36465810572578E-04_SP, &
    8.72790805146194E-04_SP, 6.22614862573135E-04_SP, 2.85998154194304E-04_SP, &
    3.84737672879366E-06_SP, -1.87906003636972E-04_SP, -2.97603646594555E-04_SP, &
    -3.45998126832656E-04_SP, -3.53382470916038E-04_SP, -3.35715635775049E-04_SP, &
    -3.04321124789040E-04_SP, -2.66722723047613E-04_SP, -2.27654214122820E-04_SP, &
    -1.89922611854562E-04_SP, -1.55058918599094E-04_SP, -1.23778240761874E-04_SP, &
    -9.62926147717644E-05_SP, -7.25178327714425E-05_SP, -5.22070028895634E-05_SP, &
    -3.50347750511901E-05_SP, -2.06489761035552E-05_SP, -8.70106096849767E-06_SP, &
    1.13698686675100E-06_SP, 9.16426474122779E-06_SP, 1.56477785428873E-05_SP, &
    2.08223629482467E-05_SP ], [26,1] )
  REAL(SP), PARAMETER :: alfa(26,4) = RESHAPE( [ alfa1, alfa2 ], [26,4] )
  REAL(SP), PARAMETER :: beta(26,5) = RESHAPE( [ beta1, beta2, beta3 ], [26,5] )
  REAL(SP), PARAMETER :: gama(26) = [ 6.29960524947437E-01_SP, 2.51984209978975E-01_SP, &
    1.54790300415656E-01_SP, 1.10713062416159E-01_SP, 8.57309395527395E-02_SP, &
    6.97161316958684E-02_SP, 5.86085671893714E-02_SP, 5.04698873536311E-02_SP, &
    4.42600580689155E-02_SP, 3.93720661543510E-02_SP, 3.54283195924455E-02_SP, &
    3.21818857502098E-02_SP, 2.94646240791158E-02_SP, 2.71581677112934E-02_SP, &
    2.51768272973862E-02_SP, 2.34570755306079E-02_SP, 2.19508390134907E-02_SP, &
    2.06210828235646E-02_SP, 1.94388240897881E-02_SP, 1.83810633800683E-02_SP, &
    1.74293213231963E-02_SP, 1.65685837786612E-02_SP, 1.57865285987918E-02_SP, &
    1.50729501494096E-02_SP, 1.44193250839955E-02_SP, 1.38184805735342E-02_SP ]
  !* FIRST EXECUTABLE STATEMENT  ASYJY
  ta = eps_2_sp
  tol = MAX(ta,1.0E-15_SP)
  tb = log10_radix_sp
  ju = min_exp_sp
  IF( Flgjy==1._SP ) THEN
    elim = -2.303_SP*(tb*ju+3._SP)
  ELSE
    jr = digits_sp
    elim = -2.303_SP*tb*(ju+jr)
  END IF
  fn = Fnu
  Iflw = 0
  DO jn = 1, In
    xx = X/fn
    Wk(1) = 1._SP - xx*xx
    abw2 = ABS(Wk(1))
    Wk(2) = SQRT(abw2)
    Wk(7) = fn**con2
    IF( abw2>0.27750_SP ) THEN
      !
      upol(1) = 1._SP
      tau = 1._SP/Wk(2)
      t2 = 1._SP/Wk(1)
      IF( Wk(1)>=0._SP ) THEN
        !
        !     CASES FOR (X/FN)<SQRT(0.7225)
        !
        Wk(3) = ABS(LOG((1._SP+Wk(2))/xx)-Wk(2))
        Wk(4) = Wk(3)*fn
        rcz = con1/Wk(4)
        IF( Wk(4)<=elim ) THEN
          z32 = 1.5_SP*Wk(3)
          rtz = z32**con2
          Wk(7) = fn**con2
          Wk(5) = rtz*Wk(7)
          Wk(6) = Wk(5)*Wk(5)
          GOTO 100
        END IF
      ELSE
        !
        !     CASES FOR (X/FN)>SQRT(1.2775)
        !
        Wk(3) = ABS(Wk(2)-ATAN(Wk(2)))
        Wk(4) = Wk(3)*fn
        rcz = -con1/Wk(4)
        z32 = 1.5_SP*Wk(3)
        rtz = z32**con2
        Wk(5) = rtz*Wk(7)
        Wk(6) = -Wk(5)*Wk(5)
        GOTO 100
      END IF
    ELSE
      !
      !     ASYMPTOTIC EXPANSION
      !     CASES NEAR X=FN, ABS(1.-(X/FN)**2)<=0.2775
      !     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES
      !
      !     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES
      !
      !     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)
      !
      sa = 0._SP
      IF( abw2/=0._SP ) sa = tols/LOG(abw2)
      sb = sa
      DO i = 1, 5
        akm = MAX(sa,2._SP)
        kmax(i) = INT(akm)
        sa = sa + sb
      END DO
      kb = kmax(5)
      klast = kb - 1
      sa = gama(kb)
      DO k = 1, klast
        kb = kb - 1
        sa = sa*Wk(1) + gama(kb)
      END DO
      z = Wk(1)*sa
      az = ABS(z)
      rtz = SQRT(az)
      Wk(3) = con1*az*rtz
      Wk(4) = Wk(3)*fn
      Wk(5) = rtz*Wk(7)
      Wk(6) = -Wk(5)*Wk(5)
      IF( z>0._SP ) THEN
        IF( Wk(4)>elim ) GOTO 50
        Wk(6) = -Wk(6)
      END IF
      phi = SQRT(SQRT(sa+sa+sa+sa))
      !
      !     B(ZETA) FOR S=0
      !
      kb = kmax(5)
      klast = kb - 1
      sb = beta(kb,1)
      DO k = 1, klast
        kb = kb - 1
        sb = sb*Wk(1) + beta(kb,1)
      END DO
      ksp1 = 1
      fn2 = fn*fn
      rfn2 = 1._SP/fn2
      rden = 1._SP
      asum = 1._SP
      relb = tol*ABS(sb)
      bsum = sb
      DO ks = 1, 4
        ksp1 = ksp1 + 1
        rden = rden*rfn2
        !
        !     A(ZETA) AND B(ZETA) FOR S=1,2,3,4
        !
        kstemp = 5 - ks
        kb = kmax(kstemp)
        klast = kb - 1
        sa = alfa(kb,ks)
        sb = beta(kb,ksp1)
        DO k = 1, klast
          kb = kb - 1
          sa = sa*Wk(1) + alfa(kb,ks)
          sb = sb*Wk(1) + beta(kb,ksp1)
        END DO
        ta = sa*rden
        tb = sb*rden
        asum = asum + ta
        bsum = bsum + tb
        IF( ABS(ta)<=tol .AND. ABS(tb)<=relb ) EXIT
      END DO
      bsum = bsum/(fn*Wk(7))
      GOTO 150
    END IF
    !
    50  Iflw = 1
    RETURN
    100  phi = SQRT((rtz+rtz)*tau)
    tb = 1._SP
    asum = 1._SP
    tfn = tau/fn
    rden = 1._SP/fn
    rfn2 = rden*rden
    rden = 1._SP
    upol(2) = (c(1)*t2+c(2))*tfn
    crz32 = con548*rcz
    bsum = upol(2) + crz32
    relb = tol*ABS(bsum)
    ap = tfn
    ks = 0
    kp1 = 2
    rzden = rcz
    l = 2
    iseta = 0
    isetb = 0
    DO lr = 2, 8, 2
      !
      !     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)
      !
      lrp1 = lr + 1
      DO k = lr, lrp1
        ks = ks + 1
        kp1 = kp1 + 1
        l = l + 1
        s1 = c(l)
        DO j = 2, kp1
          l = l + 1
          s1 = s1*t2 + c(l)
        END DO
        ap = ap*tfn
        upol(kp1) = ap*s1
        cr(ks) = br(ks)*rzden
        rzden = rzden*rcz
        dr(ks) = ar(ks)*rzden
      END DO
      suma = upol(lrp1)
      sumb = upol(lr+2) + upol(lrp1)*crz32
      ju = lrp1
      DO jr = 1, lr
        ju = ju - 1
        suma = suma + cr(jr)*upol(ju)
        sumb = sumb + dr(jr)*upol(ju)
      END DO
      rden = rden*rfn2
      tb = -tb
      IF( Wk(1)>0._SP ) tb = ABS(tb)
      IF( rden<tol ) THEN
        IF( iseta/=1 ) THEN
          IF( ABS(suma)<tol ) iseta = 1
          asum = asum + suma*tb
        END IF
        IF( isetb/=1 ) THEN
          IF( ABS(sumb)<relb ) isetb = 1
          bsum = bsum + sumb*tb
        END IF
        IF( iseta==1 .AND. isetb==1 ) EXIT
      ELSE
        asum = asum + suma*tb
        bsum = bsum + sumb*tb
      END IF
    END DO
    tb = Wk(5)
    IF( Wk(1)>0._SP ) tb = -tb
    bsum = bsum/tb
    !
    150  CALL FUNJY(Wk(6),Wk(5),Wk(4),fi,dfi)
    ta = 1._SP/tol
    tb = tiny_sp*ta*1.0E+3_SP
    IF( ABS(fi)<=tb ) THEN
      fi = fi*ta
      dfi = dfi*ta
      phi = phi*tol
    END IF
    Y(jn) = Flgjy*phi*(fi*asum+dfi*bsum)/Wk(7)
    fn = fn - Flgjy
  END DO
END SUBROUTINE ASYJY
