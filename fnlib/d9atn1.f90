!*==D9ATN1.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9ATN1
REAL(8) FUNCTION D9ATN1(X)
  IMPLICIT NONE
  !*--D9ATN15
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , ntatn1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9ATN1
  !***SUBSIDIARY
  !***PURPOSE  Evaluate DATAN(X) from first order relative accuracy so
  !            that DATAN(X) = X + X**3*D9ATN1(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      DOUBLE PRECISION (R9ATN1-S, D9ATN1-D)
  !***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB,
  !             TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate  DATAN(X)  from first order, that is, evaluate
  ! (DATAN(X)-X)/X**3  with relative error accuracy so that
  !        DATAN(X) = X + X**3*D9ATN1(X).
  !
  ! Series for ATN1       on the interval  0.          to  1.00000E+00
  !                                        with weighted error   3.39E-32
  !                                         log weighted error  31.47
  !                               significant figures required  30.26
  !                                    decimal places required  32.27
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891115  Corrected third argument in reference to INITDS.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  D9ATN1
  REAL(8) :: X , xbig , xmax , xsml , y , atn1cs(40) , eps , DCSEVL , &
    D1MACH
  LOGICAL first
  SAVE atn1cs , ntatn1 , xsml , xbig , xmax , first
  DATA atn1cs(1)/ - .3283997535355202356907939922990D-1/
  DATA atn1cs(2)/ + .5833432343172412449951669914907D-1/
  DATA atn1cs(3)/ - .7400369696719646463809011551413D-2/
  DATA atn1cs(4)/ + .1009784199337288083590357511639D-2/
  DATA atn1cs(5)/ - .1439787163565205621471303697700D-3/
  DATA atn1cs(6)/ + .2114512648992107572072112243439D-4/
  DATA atn1cs(7)/ - .3172321074254667167402564996757D-5/
  DATA atn1cs(8)/ + .4836620365460710825377859384800D-6/
  DATA atn1cs(9)/ - .7467746546814112670437614322776D-7/
  DATA atn1cs(10)/ + .1164800896824429830620998641342D-7/
  DATA atn1cs(11)/ - .1832088370847201392699956242452D-8/
  DATA atn1cs(12)/ + .2901908277966063313175351230455D-9/
  DATA atn1cs(13)/ - .4623885312106326738351805721512D-10/
  DATA atn1cs(14)/ + .7405528668775736917992197048286D-11/
  DATA atn1cs(15)/ - .1191354457845136682370820373417D-11/
  DATA atn1cs(16)/ + .1924090144391772599867855692518D-12/
  DATA atn1cs(17)/ - .3118271051076194272254476155327D-13/
  DATA atn1cs(18)/ + .5069240036567731789694520593032D-14/
  DATA atn1cs(19)/ - .8263694719802866053818284405964D-15/
  DATA atn1cs(20)/ + .1350486709817079420526506123029D-15/
  DATA atn1cs(21)/ - .2212023650481746045840137823191D-16/
  DATA atn1cs(22)/ + .3630654747381356783829047647709D-17/
  DATA atn1cs(23)/ - .5970345328847154052451215859165D-18/
  DATA atn1cs(24)/ + .9834816050077133119448329005738D-19/
  DATA atn1cs(25)/ - .1622655075855062336144387604480D-19/
  DATA atn1cs(26)/ + .2681186176945436796301320301226D-20/
  DATA atn1cs(27)/ - .4436309706785255479636243688106D-21/
  DATA atn1cs(28)/ + .7349691897652496945072465510400D-22/
  DATA atn1cs(29)/ - .1219077508350052588289401378133D-22/
  DATA atn1cs(30)/ + .2024298836805215403184540876799D-23/
  DATA atn1cs(31)/ - .3364871555797354579925576362666D-24/
  DATA atn1cs(32)/ + .5598673968346988749492933973333D-25/
  DATA atn1cs(33)/ - .9323939267272320229628532053333D-26/
  DATA atn1cs(34)/ + .1554133116995970222934807893333D-26/
  DATA atn1cs(35)/ - .2592569534179745922757427199999D-27/
  DATA atn1cs(36)/ + .4328193466245734685037909333333D-28/
  DATA atn1cs(37)/ - .7231013125595437471192405333333D-29/
  DATA atn1cs(38)/ + .1208902859830494772942165333333D-29/
  DATA atn1cs(39)/ - .2022404543449897579315199999999D-30/
  DATA atn1cs(40)/ + .3385428713046493843073706666666D-31/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9ATN1
  IF ( first ) THEN
    eps = D1MACH(3)
    ntatn1 = INITDS(atn1cs,40,0.1*REAL(eps))
    !
    xsml = SQRT(0.1D0*eps)
    xbig = 1.571D0/SQRT(eps)
    xmax = 1.571D0/eps
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1.0D0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','D9ATN1',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',&
      2,2)
    IF ( y>xbig ) CALL XERMSG('SLATEC','D9ATN1',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'&
      ,1,1)
    !
    D9ATN1 = (ATAN(X)-X)/X**3
    GOTO 99999
  ENDIF
  !
  IF ( y<=xsml ) D9ATN1 = -1.0D0/3.0D0
  IF ( y<=xsml ) RETURN
  !
  D9ATN1 = -0.25D0 + DCSEVL(2.D0*y*y-1.D0,atn1cs,ntatn1)
  RETURN
  !
  99999 END FUNCTION D9ATN1
