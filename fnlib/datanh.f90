!*==DATANH.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DATANH
DOUBLE PRECISION FUNCTION DATANH(X)
  IMPLICIT NONE
  !*--DATANH5
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , nterms
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DATANH
  !***PURPOSE  Compute the arc hyperbolic tangent.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      DOUBLE PRECISION (ATANH-S, DATANH-D, CATANH-C)
  !***KEYWORDS  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS,
  !             FNLIB, INVERSE HYPERBOLIC TANGENT
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DATANH(X) calculates the double precision arc hyperbolic
  ! tangent for double precision argument X.
  !
  ! Series for ATNH       on the interval  0.          to  2.50000E-01
  !                                        with weighted error   6.86E-32
  !                                         log weighted error  31.16
  !                               significant figures required  30.00
  !                                    decimal places required  31.88
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DATANH
  DOUBLE PRECISION X , atnhcs(27) , dxrel , sqeps , y , DCSEVL , D1MACH
  LOGICAL first
  SAVE atnhcs , nterms , dxrel , sqeps , first
  DATA atnhcs(1)/ + .9439510239319549230842892218633D-1/
  DATA atnhcs(2)/ + .4919843705578615947200034576668D-1/
  DATA atnhcs(3)/ + .2102593522455432763479327331752D-2/
  DATA atnhcs(4)/ + .1073554449776116584640731045276D-3/
  DATA atnhcs(5)/ + .5978267249293031478642787517872D-5/
  DATA atnhcs(6)/ + .3505062030889134845966834886200D-6/
  DATA atnhcs(7)/ + .2126374343765340350896219314431D-7/
  DATA atnhcs(8)/ + .1321694535715527192129801723055D-8/
  DATA atnhcs(9)/ + .8365875501178070364623604052959D-10/
  DATA atnhcs(10)/ + .5370503749311002163881434587772D-11/
  DATA atnhcs(11)/ + .3486659470157107922971245784290D-12/
  DATA atnhcs(12)/ + .2284549509603433015524024119722D-13/
  DATA atnhcs(13)/ + .1508407105944793044874229067558D-14/
  DATA atnhcs(14)/ + .1002418816804109126136995722837D-15/
  DATA atnhcs(15)/ + .6698674738165069539715526882986D-17/
  DATA atnhcs(16)/ + .4497954546494931083083327624533D-18/
  DATA atnhcs(17)/ + .3032954474279453541682367146666D-19/
  DATA atnhcs(18)/ + .2052702064190936826463861418666D-20/
  DATA atnhcs(19)/ + .1393848977053837713193014613333D-21/
  DATA atnhcs(20)/ + .9492580637224576971958954666666D-23/
  DATA atnhcs(21)/ + .6481915448242307604982442666666D-24/
  DATA atnhcs(22)/ + .4436730205723615272632320000000D-25/
  DATA atnhcs(23)/ + .3043465618543161638912000000000D-26/
  DATA atnhcs(24)/ + .2091881298792393474047999999999D-27/
  DATA atnhcs(25)/ + .1440445411234050561365333333333D-28/
  DATA atnhcs(26)/ + .9935374683141640465066666666666D-30/
  DATA atnhcs(27)/ + .6863462444358260053333333333333D-31/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DATANH
  IF ( first ) THEN
    nterms = INITDS(atnhcs,27,0.1*REAL(D1MACH(3)))
    dxrel = SQRT(D1MACH(4))
    sqeps = SQRT(3.0D0*D1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>=1.D0 ) CALL XERMSG('SLATEC','DATANH','ABS(X) GE 1',2,2)
  !
  IF ( 1.D0-y<dxrel ) CALL XERMSG('SLATEC','DATANH',&
    'ANSWER LT HALF PRECISION BECAUSE ABS(X) TOO NEAR 1'&
    ,1,1)
  !
  DATANH = X
  IF ( y>sqeps.AND.y<=0.5D0 )&
    DATANH = X*(1.0D0+DCSEVL(8.D0*X*X-1.D0,atnhcs,nterms))
  IF ( y>0.5D0 ) DATANH = 0.5D0*LOG((1.0D0+X)/(1.0D0-X))
  !
END FUNCTION DATANH
