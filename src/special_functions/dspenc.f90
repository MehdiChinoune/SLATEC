!** DSPENC
REAL(8) FUNCTION DSPENC(X)
  !>
  !  Compute a form of Spence's integral due to K. Mitchell.
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
  ! For ABS(X) .LE. 1, the uniformly convergent expansion
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
  USE service, ONLY : D1MACH
  REAL(8) :: X, aln
  INTEGER, SAVE :: nspenc
  REAL(8), SAVE :: xbig
  REAL(8), PARAMETER :: spencs(38) = [ +.1527365598892405872946684910028D+0, &
    +.8169658058051014403501838185271D-1, +.5814157140778730872977350641182D-2, &
    +.5371619814541527542247889005319D-3, +.5724704675185826233210603054782D-4, &
    +.6674546121649336343607835438589D-5, +.8276467339715676981584391689011D-6, &
    +.1073315673030678951270005873354D-6, +.1440077294303239402334590331513D-7, &
    +.1984442029965906367898877139608D-8, +.2794005822163638720201994821615D-9, &
    +.4003991310883311823072580445908D-10, +.5823462892044638471368135835757D-11, &
    +.8576708692638689278097914771224D-12, +.1276862586280193045989483033433D-12, &
    +.1918826209042517081162380416062D-13, +.2907319206977138177795799719673D-14, &
    +.4437112685276780462557473641745D-15, +.6815727787414599527867359135607D-16, &
    +.1053017386015574429547019416644D-16, +.1635389806752377100051821734570D-17, &
    +.2551852874940463932310901642581D-18, +.3999020621999360112770470379519D-19, &
    +.6291501645216811876514149171199D-20, +.9933827435675677643803887752533D-21, &
    +.1573679570749964816721763805866D-21, +.2500595316849476129369270954666D-22, &
    +.3984740918383811139210663253333D-23, +.6366473210082843892691326293333D-24, &
    +.1019674287239678367077061973333D-24, +.1636881058913518841111074133333D-25, &
    +.2633310439417650117345279999999D-26, +.4244811560123976817224362666666D-27, &
    +.6855411983680052916824746666666D-28, +.1109122433438056434018986666666D-28, &
    +.1797431304999891457365333333333D-29, +.2917505845976095173290666666666D-30, &
    +.4742646808928671061333333333333D-31 ]
  REAL(8), PARAMETER :: pi26 = +1.644934066848226436472415166646025189219D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DSPENC
  IF ( first ) THEN
    nspenc = INITDS(spencs,38,0.1*REAL(D1MACH(3)))
    xbig = 1.0D0/D1MACH(3)
    first = .FALSE.
  END IF
  !
  IF ( X>2.0D0 ) THEN
    !
    ! X .GT. 2.0
    !
    DSPENC = 2.0D0*pi26 - 0.5D0*LOG(X)**2
    IF ( X<xbig ) DSPENC = DSPENC - (1.D0+DCSEVL(4.D0/X-1.D0,spencs,nspenc))/X
    RETURN
  ELSEIF ( X<=1.0D0 ) THEN
    IF ( X>0.5D0 ) THEN
      !
      ! 0.5 .LT. X .LE. 1.0
      !
      DSPENC = pi26
      IF ( X/=1.D0 ) DSPENC = pi26 - LOG(X)*LOG(1.0D0-X) - (1.D0-X)&
        *(1.D0+DCSEVL(4.D0*(1.D0-X)-1.D0,spencs,nspenc))
      RETURN
    ELSEIF ( X>=0.0D0 ) THEN
      !
      ! 0.0 .LE. X .LE. 0.5
      !
      DSPENC = X*(1.D0+DCSEVL(4.D0*X-1.D0,spencs,nspenc))
      RETURN
    ELSEIF ( X>(-1.D0) ) THEN
      !
      ! -1.0 .LT. X .LT. 0.0
      !
      DSPENC = -0.5D0*LOG(1.0D0-X)&
        **2 - X*(1.D0+DCSEVL(4.D0*X/(X-1.D0)-1.D0,spencs,nspenc))/(X-1.D0)
      RETURN
    ELSE
      !
      ! HERE IF X .LE. -1.0
      !
      aln = LOG(1.0D0-X)
      DSPENC = -pi26 - 0.5D0*aln*(2.0D0*LOG(-X)-aln)
      IF ( X>(-xbig) ) DSPENC = DSPENC + &
        (1.D0+DCSEVL(4.D0/(1.D0-X)-1.D0,spencs,nspenc))/(1.D0-X)
      RETURN
    END IF
  END IF
  !
  ! 1.0 .LT. X .LE. 2.0
  !
  DSPENC = pi26 - 0.5D0*LOG(X)*LOG((X-1.D0)**2/X) + (X-1.D0)&
    *(1.D0+DCSEVL(4.D0*(X-1.D0)/X-1.D0,spencs,nspenc))/X
  RETURN
END FUNCTION DSPENC
