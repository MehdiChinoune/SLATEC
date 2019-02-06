!*==PSI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK PSI
      FUNCTION PSI(X)
      IMPLICIT NONE
!*--PSI5
!*** Start of declarations inserted by SPAG
      REAL apsics , aux , COT , CSEVL , dxrel , pi , PSI , psics , R1MACH , X , 
     &     xbig , y
      INTEGER i , INITS , n , ntapsi , ntpsi
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  PSI
!***PURPOSE  Compute the Psi (or Digamma) function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7C
!***TYPE      SINGLE PRECISION (PSI-S, DPSI-D, CPSI-C)
!***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
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
!***REFERENCES  (NONE)
!***ROUTINES CALLED  COT, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  PSI
      DIMENSION psics(23) , apsics(16)
      LOGICAL first
      EXTERNAL COT
      SAVE psics , apsics , pi , ntpsi , ntapsi , xbig , dxrel , first
      DATA psics(1)/ - .038057080835217922E0/
      DATA psics(2)/.49141539302938713E0/
      DATA psics(3)/ - .056815747821244730E0/
      DATA psics(4)/.008357821225914313E0/
      DATA psics(5)/ - .001333232857994342E0/
      DATA psics(6)/.000220313287069308E0/
      DATA psics(7)/ - .000037040238178456E0/
      DATA psics(8)/.000006283793654854E0/
      DATA psics(9)/ - .000001071263908506E0/
      DATA psics(10)/.000000183128394654E0/
      DATA psics(11)/ - .000000031353509361E0/
      DATA psics(12)/.000000005372808776E0/
      DATA psics(13)/ - .000000000921168141E0/
      DATA psics(14)/.000000000157981265E0/
      DATA psics(15)/ - .000000000027098646E0/
      DATA psics(16)/.000000000004648722E0/
      DATA psics(17)/ - .000000000000797527E0/
      DATA psics(18)/.000000000000136827E0/
      DATA psics(19)/ - .000000000000023475E0/
      DATA psics(20)/.000000000000004027E0/
      DATA psics(21)/ - .000000000000000691E0/
      DATA psics(22)/.000000000000000118E0/
      DATA psics(23)/ - .000000000000000020E0/
      DATA apsics(1)/ - .0204749044678185E0/
      DATA apsics(2)/ - .0101801271534859E0/
      DATA apsics(3)/.0000559718725387E0/
      DATA apsics(4)/ - .0000012917176570E0/
      DATA apsics(5)/.0000000572858606E0/
      DATA apsics(6)/ - .0000000038213539E0/
      DATA apsics(7)/.0000000003397434E0/
      DATA apsics(8)/ - .0000000000374838E0/
      DATA apsics(9)/.0000000000048990E0/
      DATA apsics(10)/ - .0000000000007344E0/
      DATA apsics(11)/.0000000000001233E0/
      DATA apsics(12)/ - .0000000000000228E0/
      DATA apsics(13)/.0000000000000045E0/
      DATA apsics(14)/ - .0000000000000009E0/
      DATA apsics(15)/.0000000000000002E0/
      DATA apsics(16)/ - .0000000000000000E0/
      DATA pi/3.14159265358979324E0/
      DATA first/.TRUE./
!***FIRST EXECUTABLE STATEMENT  PSI
      IF ( first ) THEN
        ntpsi = INITS(psics,23,0.1*R1MACH(3))
        ntapsi = INITS(apsics,16,0.1*R1MACH(3))
!
        xbig = 1.0/SQRT(R1MACH(3))
        dxrel = SQRT(R1MACH(4))
      ENDIF
      first = .FALSE.
!
      y = ABS(X)
      IF ( y>=2.0 ) THEN
!
! PSI(X) FOR ABS(X) .GE. 2.
!
        aux = 0.
        IF ( y<xbig ) aux = CSEVL(8./y**2-1.,apsics,ntapsi)
        IF ( X<0. ) PSI = LOG(ABS(X)) - 0.5/X + aux - pi*COT(pi*X)
        IF ( X>0. ) PSI = LOG(X) - 0.5/X + aux
        GOTO 99999
      ENDIF
!
! PSI(X) FOR -2. .LT. X .LT. 2.
!
      n = X
      IF ( X<0. ) n = n - 1
      y = X - n
      n = n - 1
      PSI = CSEVL(2.*y-1.,psics,ntpsi)
      IF ( n==0 ) RETURN
!
      n = -n
      IF ( X==0. ) CALL XERMSG('SLATEC','PSI','X IS 0',2,2)
      IF ( X<0..AND.X+n-2==0. )
     &      CALL XERMSG('SLATEC','PSI','X IS A NEGATIVE INTEGER',3,2)
      IF ( X<(-0.5).AND.ABS((X-AINT(X-0.5))/X)<dxrel )
     &      CALL XERMSG('SLATEC','PSI',
     &     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,1)
!
      DO i = 1 , n
        PSI = PSI - 1.0/(X+i-1)
      ENDDO
      RETURN
!
99999 END FUNCTION PSI
