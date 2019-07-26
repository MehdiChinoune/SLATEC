!** DSPENC
REAL(DP) ELEMENTAL FUNCTION DSPENC(X)
  !> Compute a form of Spence's integral due to K. Mitchell.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      DOUBLE PRECISION (SPENC-S, DSPENC-D)
  !***
  ! **Keywords:**  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DSPENC(X) calculates the double precision Spence's integral
  ! for double precision argument X.  Spence's function defined by
  !        integral from 0 to X of  -LOG(1-Y)/Y  DY.
  ! For ABS(X) <= 1, the uniformly convergent expansion
  !        DSPENC = sum K=1,infinity  X**K / K**2     is valid.
  ! This is a form of Spence's integral due to K. Mitchell which differs
  ! from the definition in the NBS Handbook of Mathematical Functions.
  !
  ! Spence's function can be used to evaluate much more general integral
  ! forms.  For example,
  !        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  =
  !             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C
  !             - DSPENC (A*(C*Z+D)/(A*D-B*C)) / C.
  !
  ! Ref -- K. Mitchell, Philosophical Magazine, 40, p.351 (1949).
  !        Stegun and Abromowitz, AMS 55, p.1004.
  !
  !
  ! Series for SPEN       on the interval  0.          to  5.00000E-01
  !                                        with weighted error   4.74E-32
  !                                         log weighted error  31.32
  !                               significant figures required  30.37
  !                                    decimal places required  32.11
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, INITDS

  !* REVISION HISTORY  (YYMMDD)
  !   780201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891115  Corrected third argument in reference to INITDS.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : eps_2_dp
  !
  REAL(DP), INTENT(IN) :: X
  !
  REAL(DP) :: aln
  INTEGER, PARAMETER :: nspenc = 20
  REAL(DP), PARAMETER :: xbig = 1._DP/eps_2_dp
  REAL(DP), PARAMETER :: spencs(38) = [ +.1527365598892405872946684910028E+0_DP, &
    +.8169658058051014403501838185271E-1_DP, +.5814157140778730872977350641182E-2_DP, &
    +.5371619814541527542247889005319E-3_DP, +.5724704675185826233210603054782E-4_DP, &
    +.6674546121649336343607835438589E-5_DP, +.8276467339715676981584391689011E-6_DP, &
    +.1073315673030678951270005873354E-6_DP, +.1440077294303239402334590331513E-7_DP, &
    +.1984442029965906367898877139608E-8_DP, +.2794005822163638720201994821615E-9_DP, &
    +.4003991310883311823072580445908E-10_DP, +.5823462892044638471368135835757E-11_DP, &
    +.8576708692638689278097914771224E-12_DP, +.1276862586280193045989483033433E-12_DP, &
    +.1918826209042517081162380416062E-13_DP, +.2907319206977138177795799719673E-14_DP, &
    +.4437112685276780462557473641745E-15_DP, +.6815727787414599527867359135607E-16_DP, &
    +.1053017386015574429547019416644E-16_DP, +.1635389806752377100051821734570E-17_DP, &
    +.2551852874940463932310901642581E-18_DP, +.3999020621999360112770470379519E-19_DP, &
    +.6291501645216811876514149171199E-20_DP, +.9933827435675677643803887752533E-21_DP, &
    +.1573679570749964816721763805866E-21_DP, +.2500595316849476129369270954666E-22_DP, &
    +.3984740918383811139210663253333E-23_DP, +.6366473210082843892691326293333E-24_DP, &
    +.1019674287239678367077061973333E-24_DP, +.1636881058913518841111074133333E-25_DP, &
    +.2633310439417650117345279999999E-26_DP, +.4244811560123976817224362666666E-27_DP, &
    +.6855411983680052916824746666666E-28_DP, +.1109122433438056434018986666666E-28_DP, &
    +.1797431304999891457365333333333E-29_DP, +.2917505845976095173290666666666E-30_DP, &
    +.4742646808928671061333333333333E-31_DP ]
  REAL(DP), PARAMETER :: pi26 = +1.644934066848226436472415166646025189219_DP
  !* FIRST EXECUTABLE STATEMENT  DSPENC
  ! nspenc = INITDS(spencs,0.1_SP*eps_2_dp)
  !
  IF( X>2._DP ) THEN
    ! X > 2.0
    DSPENC = 2._DP*pi26 - 0.5_DP*LOG(X)**2
    IF( X<xbig ) DSPENC = DSPENC - (1._DP+DCSEVL(4._DP/X-1._DP,spencs(1:nspenc)))/X
  ELSEIF( X>1._DP ) THEN
    ! 1.0 < X <= 2.0
    DSPENC = pi26 - 0.5_DP*LOG(X)*LOG((X-1._DP)**2/X) + (X-1._DP)&
      *(1._DP+DCSEVL(4._DP*(X-1._DP)/X-1._DP,spencs(1:nspenc)))/X
  ELSEIF( X>0.5_DP ) THEN
    ! 0.5 < X <= 1.0
    DSPENC = pi26
    IF( X/=1._DP ) DSPENC = pi26 - LOG(X)*LOG(1._DP-X) - (1._DP-X)&
      *(1._DP+DCSEVL(4._DP*(1._DP-X)-1._DP,spencs(1:nspenc)))
  ELSEIF( X>=0._DP ) THEN
    ! 0.0 <= X <= 0.5
    DSPENC = X*(1._DP+DCSEVL(4._DP*X-1._DP,spencs(1:nspenc)))
    RETURN
  ELSEIF( X>(-1._DP) ) THEN
    ! -1.0 < X < 0.0
    DSPENC = -0.5_DP*LOG(1._DP-X)**2 &
      - X*(1._DP+DCSEVL(4._DP*X/(X-1._DP)-1._DP,spencs(1:nspenc)))/(X-1._DP)
  ELSE
    ! HERE IF X <= -1.0
    aln = LOG(1._DP-X)
    DSPENC = -pi26 - 0.5_DP*aln*(2._DP*LOG(-X)-aln)
    IF( X>(-xbig) ) DSPENC = DSPENC + &
      (1._DP+DCSEVL(4._DP/(1._DP-X)-1._DP,spencs(1:nspenc)))/(1._DP-X)
  END IF

  RETURN
END FUNCTION DSPENC