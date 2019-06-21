!** PSI
REAL(SP) FUNCTION PSI(X)
  !> Compute the Psi (or Digamma) function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7C
  !***
  ! **Type:**      SINGLE PRECISION (PSI-S, DPSI-D, CPSI-C)
  !***
  ! **Keywords:**  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! PSI(X) calculates the psi (or digamma) function for real argument X.
  ! PSI(X) is the logarithmic derivative of the gamma function of X.
  !
  ! Series for PSI        on the interval  0.          to  1.00000D+00
  !                                        with weighted error   2.03E-17
  !                                         log weighted error  16.69
  !                               significant figures required  16.39
  !                                    decimal places required  17.37
  !
  ! Series for APSI       on the interval  0.          to  2.50000D-01
  !                                        with weighted error   5.54E-17
  !                                         log weighted error  16.26
  !                               significant figures required  14.42
  !                                    decimal places required  16.86
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  COT, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL(SP) :: X
  INTEGER :: i, n
  REAL(SP) :: aux, y
  INTEGER, SAVE :: ntpsi, ntapsi
  REAL(SP), PARAMETER :: xbig = 1._SP/SQRT(R1MACH(3)), dxrel = SQRT(R1MACH(4))
  REAL(SP), PARAMETER :: psics(23) = [ -.038057080835217922_SP, .49141539302938713_SP, &
    -.056815747821244730_SP, .008357821225914313_SP,-.001333232857994342_SP, &
    .000220313287069308_SP, -.000037040238178456_SP, .000006283793654854_SP, &
    -.000001071263908506_SP, .000000183128394654_SP,-.000000031353509361_SP, &
    .000000005372808776_SP, -.000000000921168141_SP, .000000000157981265_SP, &
    -.000000000027098646_SP, .000000000004648722_SP,-.000000000000797527_SP, &
    .000000000000136827_SP, -.000000000000023475_SP, .000000000000004027_SP, &
    -.000000000000000691_SP, .000000000000000118_SP,-.000000000000000020_SP ]
  REAL(SP), PARAMETER :: apsics(16) = [ -.0204749044678185_SP, -.0101801271534859_SP, &
    .0000559718725387_SP, -.0000012917176570_SP, .0000000572858606_SP, &
    -.0000000038213539_SP, .0000000003397434_SP,-.0000000000374838_SP, &
    .0000000000048990_SP, -.0000000000007344_SP, .0000000000001233_SP, &
    -.0000000000000228_SP, .0000000000000045_SP,-.0000000000000009_SP, &
    .0000000000000002_SP, -.0000000000000000_SP ]
  REAL(SP), PARAMETER :: pi = 3.14159265358979324_SP
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  PSI
  IF( first ) THEN
    ntpsi = INITS(psics,23,0.1_SP*R1MACH(3))
    ntapsi = INITS(apsics,16,0.1_SP*R1MACH(3))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF( y>=2._SP ) THEN
    !
    ! PSI(X) FOR ABS(X) >= 2.
    !
    aux = 0.
    IF( y<xbig ) aux = CSEVL(8._SP/y**2-1._SP,apsics,ntapsi)
    IF( X<0. ) THEN
      PSI = LOG(ABS(X)) - 0.5_SP/X + aux - pi*COT(pi*X)
    ELSE
      PSI = LOG(X) - 0.5_SP/X + aux
    END IF
    RETURN
  END IF
  !
  ! PSI(X) FOR -2. < X < 2.
  !
  n = INT( X )
  IF( X<0. ) n = n - 1
  y = X - n
  n = n - 1
  PSI = CSEVL(2._SP*y-1._SP,psics,ntpsi)
  IF( n==0 ) RETURN
  !
  n = -n
  IF( X==0. ) CALL XERMSG('PSI','X IS 0',2,2)
  IF( X<0. .AND. X+n-2==0. )&
    CALL XERMSG('PSI','X IS A NEGATIVE INTEGER',3,2)
  IF( X<(-0.5_SP) .AND. ABS((X-AINT(X-0.5_SP))/X)<dxrel )&
    CALL XERMSG('PSI',&
    'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,1)
  !
  DO i = 1, n
    PSI = PSI - 1._SP/(X+i-1)
  END DO
  RETURN
END FUNCTION PSI
