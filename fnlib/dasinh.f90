!*==DASINH.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DASINH
DOUBLE PRECISION FUNCTION DASINH(X)
  IMPLICIT NONE
  !*--DASINH5
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , nterms
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DASINH
  !***PURPOSE  Compute the arc hyperbolic sine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      DOUBLE PRECISION (ASINH-S, DASINH-D, CASINH-C)
  !***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC SINE
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DASINH(X) calculates the double precision arc hyperbolic
  ! sine for double precision argument X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DASINH
  DOUBLE PRECISION X , asnhcs(39) , aln2 , sqeps , xmax , y , DCSEVL , &
    D1MACH
  LOGICAL first
  SAVE asnhcs , aln2 , nterms , xmax , sqeps , first
  DATA asnhcs(1)/ - .12820039911738186343372127359268D+0/
  DATA asnhcs(2)/ - .58811761189951767565211757138362D-1/
  DATA asnhcs(3)/ + .47274654322124815640725249756029D-2/
  DATA asnhcs(4)/ - .49383631626536172101360174790273D-3/
  DATA asnhcs(5)/ + .58506207058557412287494835259321D-4/
  DATA asnhcs(6)/ - .74669983289313681354755069217188D-5/
  DATA asnhcs(7)/ + .10011693583558199265966192015812D-5/
  DATA asnhcs(8)/ - .13903543858708333608616472258886D-6/
  DATA asnhcs(9)/ + .19823169483172793547317360237148D-7/
  DATA asnhcs(10)/ - .28847468417848843612747272800317D-8/
  DATA asnhcs(11)/ + .42672965467159937953457514995907D-9/
  DATA asnhcs(12)/ - .63976084654366357868752632309681D-10/
  DATA asnhcs(13)/ + .96991686089064704147878293131179D-11/
  DATA asnhcs(14)/ - .14844276972043770830246658365696D-11/
  DATA asnhcs(15)/ + .22903737939027447988040184378983D-12/
  DATA asnhcs(16)/ - .35588395132732645159978942651310D-13/
  DATA asnhcs(17)/ + .55639694080056789953374539088554D-14/
  DATA asnhcs(18)/ - .87462509599624678045666593520162D-15/
  DATA asnhcs(19)/ + .13815248844526692155868802298129D-15/
  DATA asnhcs(20)/ - .21916688282900363984955142264149D-16/
  DATA asnhcs(21)/ + .34904658524827565638313923706880D-17/
  DATA asnhcs(22)/ - .55785788400895742439630157032106D-18/
  DATA asnhcs(23)/ + .89445146617134012551050882798933D-19/
  DATA asnhcs(24)/ - .14383426346571317305551845239466D-19/
  DATA asnhcs(25)/ + .23191811872169963036326144682666D-20/
  DATA asnhcs(26)/ - .37487007953314343674570604543999D-21/
  DATA asnhcs(27)/ + .60732109822064279404549242880000D-22/
  DATA asnhcs(28)/ - .98599402764633583177370173440000D-23/
  DATA asnhcs(29)/ + .16039217452788496315232638293333D-23/
  DATA asnhcs(30)/ - .26138847350287686596716134399999D-24/
  DATA asnhcs(31)/ + .42670849606857390833358165333333D-25/
  DATA asnhcs(32)/ - .69770217039185243299730773333333D-26/
  DATA asnhcs(33)/ + .11425088336806858659812693333333D-26/
  DATA asnhcs(34)/ - .18735292078860968933021013333333D-27/
  DATA asnhcs(35)/ + .30763584414464922794065920000000D-28/
  DATA asnhcs(36)/ - .50577364031639824787046399999999D-29/
  DATA asnhcs(37)/ + .83250754712689142224213333333333D-30/
  DATA asnhcs(38)/ - .13718457282501044163925333333333D-30/
  DATA asnhcs(39)/ + .22629868426552784104106666666666D-31/
  DATA aln2/0.69314718055994530941723212145818D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DASINH
  IF ( first ) THEN
    nterms = INITDS(asnhcs,39,0.1*REAL(D1MACH(3)))
    sqeps = SQRT(D1MACH(3))
    xmax = 1.0D0/sqeps
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1.0D0 ) THEN
    IF ( y<xmax ) DASINH = LOG(y+SQRT(y*y+1.D0))
    IF ( y>=xmax ) DASINH = aln2 + LOG(y)
    DASINH = SIGN(DASINH,X)
    GOTO 99999
  ENDIF
  !
  DASINH = X
  IF ( y>sqeps ) DASINH = X*(1.0D0+DCSEVL(2.D0*X*X-1.D0,asnhcs,nterms))
  RETURN
  !
  99999 END FUNCTION DASINH
