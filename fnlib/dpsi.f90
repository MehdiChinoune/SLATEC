!DECK DPSI
REAL(8) FUNCTION DPSI(X)
  IMPLICIT NONE
  INTEGER i, INITDS, n, ntapsi, ntpsi
  !***BEGIN PROLOGUE  DPSI
  !***PURPOSE  Compute the Psi (or Digamma) function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7C
  !***TYPE      DOUBLE PRECISION (PSI-S, DPSI-D, CPSI-C)
  !***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DPSI calculates the double precision Psi (or Digamma) function for
  ! double precision argument X.  PSI(X) is the logarithmic derivative
  ! of the Gamma function of X.
  !
  ! Series for PSI        on the interval  0.          to  1.00000E+00
  !                                        with weighted error   5.79E-32
  !                                         log weighted error  31.24
  !                               significant figures required  30.93
  !                                    decimal places required  32.05
  !
  !
  ! Series for APSI       on the interval  0.          to  1.00000E-02
  !                                        with weighted error   7.75E-33
  !                                         log weighted error  32.11
  !                               significant figures required  28.88
  !                                    decimal places required  32.71
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCOT, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable name.  (RWC, WRB)
  !***END PROLOGUE  DPSI
  REAL(8) :: X, psics(42), apsics(16), aux, dxrel, pi, xbig, &
    y, DCOT, DCSEVL, D1MACH
  LOGICAL first
  EXTERNAL DCOT
  SAVE psics, apsics, pi, ntpsi, ntapsi, xbig, dxrel, first
  DATA psics(1)/ - .38057080835217921520437677667039D-1/
  DATA psics(2)/ + .49141539302938712748204699654277D+0/
  DATA psics(3)/ - .56815747821244730242892064734081D-1/
  DATA psics(4)/ + .83578212259143131362775650747862D-2/
  DATA psics(5)/ - .13332328579943425998079274172393D-2/
  DATA psics(6)/ + .22031328706930824892872397979521D-3/
  DATA psics(7)/ - .37040238178456883592889086949229D-4/
  DATA psics(8)/ + .62837936548549898933651418717690D-5/
  DATA psics(9)/ - .10712639085061849855283541747074D-5/
  DATA psics(10)/ + .18312839465484165805731589810378D-6/
  DATA psics(11)/ - .31353509361808509869005779796885D-7/
  DATA psics(12)/ + .53728087762007766260471919143615D-8/
  DATA psics(13)/ - .92116814159784275717880632624730D-9/
  DATA psics(14)/ + .15798126521481822782252884032823D-9/
  DATA psics(15)/ - .27098646132380443065440589409707D-10/
  DATA psics(16)/ + .46487228599096834872947319529549D-11/
  DATA psics(17)/ - .79752725638303689726504797772737D-12/
  DATA psics(18)/ + .13682723857476992249251053892838D-12/
  DATA psics(19)/ - .23475156060658972717320677980719D-13/
  DATA psics(20)/ + .40276307155603541107907925006281D-14/
  DATA psics(21)/ - .69102518531179037846547422974771D-15/
  DATA psics(22)/ + .11856047138863349552929139525768D-15/
  DATA psics(23)/ - .20341689616261559308154210484223D-16/
  DATA psics(24)/ + .34900749686463043850374232932351D-17/
  DATA psics(25)/ - .59880146934976711003011081393493D-18/
  DATA psics(26)/ + .10273801628080588258398005712213D-18/
  DATA psics(27)/ - .17627049424561071368359260105386D-19/
  DATA psics(28)/ + .30243228018156920457454035490133D-20/
  DATA psics(29)/ - .51889168302092313774286088874666D-21/
  DATA psics(30)/ + .89027730345845713905005887487999D-22/
  DATA psics(31)/ - .15274742899426728392894971904000D-22/
  DATA psics(32)/ + .26207314798962083136358318079999D-23/
  DATA psics(33)/ - .44964642738220696772598388053333D-24/
  DATA psics(34)/ + .77147129596345107028919364266666D-25/
  DATA psics(35)/ - .13236354761887702968102638933333D-25/
  DATA psics(36)/ + .22709994362408300091277311999999D-26/
  DATA psics(37)/ - .38964190215374115954491391999999D-27/
  DATA psics(38)/ + .66851981388855302310679893333333D-28/
  DATA psics(39)/ - .11469986654920864872529919999999D-28/
  DATA psics(40)/ + .19679385886541405920515413333333D-29/
  DATA psics(41)/ - .33764488189750979801907200000000D-30/
  DATA psics(42)/ + .57930703193214159246677333333333D-31/
  DATA apsics(1)/ - .832710791069290760174456932269D-3/
  DATA apsics(2)/ - .416251842192739352821627121990D-3/
  DATA apsics(3)/ + .103431560978741291174463193961D-6/
  DATA apsics(4)/ - .121468184135904152987299556365D-9/
  DATA apsics(5)/ + .311369431998356155521240278178D-12/
  DATA apsics(6)/ - .136461337193177041776516100945D-14/
  DATA apsics(7)/ + .902051751315416565130837974000D-17/
  DATA apsics(8)/ - .831542997421591464829933635466D-19/
  DATA apsics(9)/ + .101224257073907254188479482666D-20/
  DATA apsics(10)/ - .156270249435622507620478933333D-22/
  DATA apsics(11)/ + .296542716808903896133226666666D-24/
  DATA apsics(12)/ - .674686886765702163741866666666D-26/
  DATA apsics(13)/ + .180345311697189904213333333333D-27/
  DATA apsics(14)/ - .556901618245983607466666666666D-29/
  DATA apsics(15)/ + .195867922607736251733333333333D-30/
  DATA apsics(16)/ - .775195892523335680000000000000D-32/
  DATA pi/3.14159265358979323846264338327950D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DPSI
  IF ( first ) THEN
    ntpsi = INITDS(psics,42,0.1*REAL(D1MACH(3)))
    ntapsi = INITDS(apsics,16,0.1*REAL(D1MACH(3)))
    !
    xbig = 1.0D0/SQRT(D1MACH(3))
    dxrel = SQRT(D1MACH(4))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  !
  IF ( y>10.0D0 ) THEN
    !
    ! DPSI(X) FOR ABS(X) .GT. 10.0
    !
    aux = 0.D0
    IF ( y<xbig ) aux = DCSEVL(2.D0*(10.D0/y)**2-1.D0,apsics,ntapsi)
    !
    IF ( X<0.D0 ) THEN
      DPSI = LOG(ABS(X)) - 0.5D0/X + aux - pi*DCOT(pi*X)
    ELSE
      DPSI = LOG(X) - 0.5D0/X + aux
    ENDIF
    RETURN
  ELSE
    !
    ! DPSI(X) FOR ABS(X) .LE. 2
    !
    n = INT( X )
    IF ( X<0.D0 ) n = n - 1
    y = X - n
    n = n - 1
    DPSI = DCSEVL(2.D0*y-1.D0,psics,ntpsi)
    IF ( n==0 ) RETURN
    !
    IF ( n<=0 ) THEN
      !
      n = -n
      IF ( X==0.D0 ) CALL XERMSG('SLATEC','DPSI','X IS 0',2,2)
      IF ( X<0.D0.AND.X+n-2==0.D0 )&
        CALL XERMSG('SLATEC','DPSI','X IS A NEGATIVE INTEGER',3,2)
      IF ( X<(-0.5D0).AND.ABS((X-AINT(X-0.5D0))/X)<dxrel )&
        CALL XERMSG('SLATEC','DPSI',&
        'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,&
        1)
      !
      DO i = 1, n
        DPSI = DPSI - 1.D0/(X+i-1)
      ENDDO
      RETURN
    ENDIF
  ENDIF
  !
  ! DPSI(X) FOR X .GE. 2.0 AND X .LE. 10.0
  !
  DO i = 1, n
    DPSI = DPSI + 1.0D0/(y+i)
  ENDDO
  RETURN
END FUNCTION DPSI
