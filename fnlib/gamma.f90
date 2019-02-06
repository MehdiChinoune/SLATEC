!*==GAMMA.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK GAMMA
      FUNCTION GAMMA(X)
      IMPLICIT NONE
!*--GAMMA5
!*** Start of declarations inserted by SPAG
      REAL CSEVL , dxrel , GAMMA , gcs , pi , R1MACH , R9LGMC , sinpiy , 
     &     sq2pil , X , xmax , xmin , y
      INTEGER i , INITS , n , ngcs
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  GAMMA
!***PURPOSE  Compute the complete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      SINGLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! GAMMA computes the gamma function at X, where X is not 0, -1, -2, ....
! GAMMA and X are single precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, GAMLIM, INITS, R1MACH, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  GAMMA
      DIMENSION gcs(23)
      LOGICAL first
      SAVE gcs , pi , sq2pil , ngcs , xmin , xmax , dxrel , first
      DATA gcs(1)/.008571195590989331E0/
      DATA gcs(2)/.004415381324841007E0/
      DATA gcs(3)/.05685043681599363E0/
      DATA gcs(4)/ - .004219835396418561E0/
      DATA gcs(5)/.001326808181212460E0/
      DATA gcs(6)/ - .0001893024529798880E0/
      DATA gcs(7)/.0000360692532744124E0/
      DATA gcs(8)/ - .0000060567619044608E0/
      DATA gcs(9)/.0000010558295463022E0/
      DATA gcs(10)/ - .0000001811967365542E0/
      DATA gcs(11)/.0000000311772496471E0/
      DATA gcs(12)/ - .0000000053542196390E0/
      DATA gcs(13)/.0000000009193275519E0/
      DATA gcs(14)/ - .0000000001577941280E0/
      DATA gcs(15)/.0000000000270798062E0/
      DATA gcs(16)/ - .0000000000046468186E0/
      DATA gcs(17)/.0000000000007973350E0/
      DATA gcs(18)/ - .0000000000001368078E0/
      DATA gcs(19)/.0000000000000234731E0/
      DATA gcs(20)/ - .0000000000000040274E0/
      DATA gcs(21)/.0000000000000006910E0/
      DATA gcs(22)/ - .0000000000000001185E0/
      DATA gcs(23)/.0000000000000000203E0/
      DATA pi/3.14159265358979324E0/
! SQ2PIL IS LOG (SQRT (2.*PI) )
      DATA sq2pil/0.91893853320467274E0/
      DATA first/.TRUE./
!
! LANL DEPENDENT CODE REMOVED 81.02.04
!
!***FIRST EXECUTABLE STATEMENT  GAMMA
      IF ( first ) THEN
!
! ---------------------------------------------------------------------
! INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF
! TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
! THAN MACHINE PRECISION.
!
        ngcs = INITS(gcs,23,0.1*R1MACH(3))
!
        CALL GAMLIM(xmin,xmax)
        dxrel = SQRT(R1MACH(4))
!
! ---------------------------------------------------------------------
! FINISH INITIALIZATION.  START EVALUATING GAMMA(X).
!
      ENDIF
      first = .FALSE.
!
      y = ABS(X)
      IF ( y>10.0 ) THEN
!
! COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
        IF ( X>xmax ) CALL XERMSG('SLATEC','GAMMA','X SO BIG GAMMA OVERFLOWS',3,
     &                            2)
!
        GAMMA = 0.
        IF ( X<xmin ) CALL XERMSG('SLATEC','GAMMA','X SO SMALL GAMMA UNDERFLOWS'
     &                            ,2,1)
        IF ( X<xmin ) RETURN
!
        GAMMA = EXP((y-0.5)*LOG(y)-y+sq2pil+R9LGMC(y))
        IF ( X>0. ) RETURN
!
        IF ( ABS((X-AINT(X-0.5))/X)<dxrel ) CALL XERMSG('SLATEC','GAMMA',
     &       'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER',1,1)
!
        sinpiy = SIN(pi*y)
        IF ( sinpiy==0. ) CALL XERMSG('SLATEC','GAMMA','X IS A NEGATIVE INTEGER'
     &                                ,4,2)
!
        GAMMA = -pi/(y*sinpiy*GAMMA)
        GOTO 99999
      ELSE
!
! COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND
! FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL.
!
        n = X
        IF ( X<0. ) n = n - 1
        y = X - n
        n = n - 1
        GAMMA = 0.9375 + CSEVL(2.*y-1.,gcs,ngcs)
        IF ( n==0 ) RETURN
!
        IF ( n<=0 ) THEN
!
! COMPUTE GAMMA(X) FOR X .LT. 1.
!
          n = -n
          IF ( X==0. ) CALL XERMSG('SLATEC','GAMMA','X IS 0',4,2)
          IF ( X<0..AND.X+n-2==0. )
     &          CALL XERMSG('SLATEC','GAMMA','X IS A NEGATIVE INTEGER',4,2)
          IF ( X<(-0.5).AND.ABS((X-AINT(X-0.5))/X)<dxrel )
     &          CALL XERMSG('SLATEC','GAMMA',
     &         'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,
     &         1)
!
          DO i = 1 , n
            GAMMA = GAMMA/(X+i-1)
          ENDDO
          RETURN
        ENDIF
      ENDIF
!
! GAMMA(X) FOR X .GE. 2.
!
      DO i = 1 , n
        GAMMA = (y+i)*GAMMA
      ENDDO
      RETURN
!
99999 END FUNCTION GAMMA
