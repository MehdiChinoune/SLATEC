!*==CHECK2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CHECK2
SUBROUTINE CHECK2(Sfac,Sdfac,Dfac,Dqfac,Kprint)
  IMPLICIT NONE
  !*--CHECK25
  !*** Start of declarations inserted by SPAG
  INTEGER i, ICAse, INCx, INCy, j, ki, kn, kni, kpar, Kprint, &
    ksize, lenx, leny, MODe, mx, my, N, NPRint
  REAL sa, sb, sc, Sdfac, SDOT, SDSDOT, Sfac, ss
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CHECK2
  !***PURPOSE  (UNKNOWN)
  !***LIBRARY   SLATEC
  !***AUTHOR  Lawson, C. L., (JPL)
  !***DESCRIPTION
  !
  !     THIS SUBPROGRAM TESTS THE BASIC LINEAR ALGEBRA SUBPROGRAMS 1-11,
  !     14-15, AND 18-25. SUBPROGRAMS IN THIS SET EACH REQUIRE TWO ARRAYS
  !     IN THE PARAMETER LIST.
  !
  !     C. L. LAWSON, JPL, 1975 FEB 26, APR 29, MAY 8, MAY 28
  !
  !***ROUTINES CALLED  CAXPY, CCOPY, CDOTC, CDOTU, CSWAP, DAXPY, DCOPY,
  !                    DDOT, DQDOTA, DQDOTI, DROT, DROTM, DSDOT, DSWAP,
  !                    DTEST, SAXPY, SCOPY, SDOT, SDSDOT, SROT, SROTM,
  !                    SSWAP, STEST
  !***COMMON BLOCKS    COMBLA
  !***REVISION HISTORY  (YYMMDD)
  !   750226  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CHECK2
  COMMON /COMBLA/ NPRint, ICAse, N, INCx, INCy, MODe, PASs
  !
  LOGICAL PASs
  INTEGER incxs(4), incys(4), lens(4,2), ns(4)
  REAL sx(7), sy(7), stx(7), sty(7), ssize1(4), ssize2(14,2)
  REAL ssize(7), qc(30), sparam(5), st7b(4,4), ssize3(4)
  REAL(8) :: dx(7), da, dx1(7), dy1(7), dy(7), dt7(4,4), &
    dt8(7,4,4)
  REAL(8) :: dx2(7), dy2(7), dt2(4,4,2), dparam(5), dpar(5,4)
  REAL(8) :: DSDOT, DDOT, DQDOTI, DQDOTA, Dfac, Dqfac
  REAL(8) :: dt10x(7,4,4), dt10y(7,4,4), db
  REAL(8) :: dsize1(4), dsize2(7,2), dsize(7)
  REAL(8) :: dc, ds, dt9x(7,4,4), dt9y(7,4,4), dtx(7), dty(7)
  REAL(8) :: dt19x(7,4,16), dt19xa(7,4,4), dt19xb(7,4,4)
  REAL(8) :: dt19xc(7,4,4), dt19xd(7,4,4), dt19y(7,4,16)
  REAL(8) :: dt19ya(7,4,4), dt19yb(7,4,4), dt19yc(7,4,4)
  REAL(8) :: dt19yd(7,4,4)
  !
  EQUIVALENCE (dt19x(1,1,1),dt19xa(1,1,1))
  EQUIVALENCE (dt19x(1,1,5),dt19xb(1,1,1))
  EQUIVALENCE (dt19x(1,1,9),dt19xc(1,1,1))
  EQUIVALENCE (dt19x(1,1,13),dt19xd(1,1,1))
  EQUIVALENCE (dt19y(1,1,1),dt19ya(1,1,1))
  EQUIVALENCE (dt19y(1,1,5),dt19yb(1,1,1))
  EQUIVALENCE (dt19y(1,1,9),dt19yc(1,1,1))
  EQUIVALENCE (dt19y(1,1,13),dt19yd(1,1,1))
  COMPLEX cx(7), ca, cx1(7), cy1(7), cy(7), ct6(4,4), ct7(4,4)
  COMPLEX ct8(7,4,4), csize1(4), csize2(7,2)
  COMPLEX ct10x(7,4,4), ct10y(7,4,4)
  COMPLEX CDOTC, CDOTU
  DATA sa, da, ca, db, sb/.3, .3D0, (.4,-.7), .25D0, .1/
  DATA incxs/1, 2, -2, -1/
  DATA incys/1, -2, 1, -2/
  DATA lens/1, 1, 2, 4, 1, 1, 3, 7/
  DATA ns/0, 1, 2, 4/
  DATA sc, ss, dc, ds/.8, .6, .8D0, .6D0/
  DATA dx1/.6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0/
  DATA dy1/.5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0/
  DATA dx2/1.D0, .01D0, .02D0, 1.D0, .06D0, 2.D0, 1.D0/
  DATA dy2/1.D0, .04D0, -.03D0, -1.D0, .05D0, 3.D0, -1.D0/
  !            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
  !            COMPUTATION AT RUN TIME.
  DATA cx1/(.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), (-.9,-.4), (.1,.4)&
    , (-.6,.6)/
  DATA cy1/(.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), (-.1,-.2), (-.5,-.3)&
    , (.8,-.7)/
  !
  !                             FOR DQDOTI AND DQDOTA
  !
  DATA dt2/0.25D0, 1.25D0, 1.2504D0, 0.2498D0, 0.25D0, 1.25D0, &
    0.24D0, 0.2492D0, 0.25D0, 1.25D0, 0.31D0, 0.2518D0, 0.25D0, &
    1.25D0, 1.2497D0, 0.2507D0, 0.D0, 2.D0, 2.0008D0, -.0004D0, &
    0.D0, 2.D0, -.02D0, -.0016D0, 0.D0, 2.D0, .12D0, .0036D0, &
    0.D0, 2.D0, 1.9994D0, .0014D0/
  DATA dt7/0.D0, .30D0, .21D0, .62D0, 0.D0, .30D0, -.07D0, .85D0, &
    0.D0, .30D0, -.79D0, -.74D0, 0.D0, .30D0, .33D0, 1.27D0/
  DATA st7b/.1, .4, .31, .72, .1, .4, .03, .95, .1, .4, -.69, &
    -.64, .1, .4, .43, 1.37/
  !
  !                       FOR CDOTU
  !
  DATA ct7/(0.,0.), (-.06,-.90), (.65,-.47), (-.34,-1.22), (0.,0.), &
    (-.06,-.90), (-.59,-1.46), (-1.04,-.04), (0.,0.), (-.06,-.90), &
    (-.83,.59), (.07,-.37), (0.,0.), (-.06,-.90), (-.76,-1.15), &
    (-1.33,-1.82)/
  !
  !                       FOR CDOTC
  !
  DATA ct6/(0.,0.), (.90,0.06), (.91,-.77), (1.80,-.10), (0.,0.), &
    (.90,0.06), (1.45,.74), (.20,.90), (0.,0.), (.90,0.06), &
    (-.55,.23), (.83,-.39), (0.,0.), (.90,0.06), (1.04,0.79), &
    (1.95,1.22)/
  !
  DATA dt8/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.87D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, .68D0, -.87D0, .15D0, .94D0, 0.D0, 0.D0, &
    0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .35D0, -.9D0, .48D0, &
    0.D0, 0.D0, 0.D0, 0.D0, .38D0, -.9D0, .57D0, .7D0, -.75D0, &
    .2D0, .98D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .68D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .35D0, -.72D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .38D0, -.63D0, .15D0, .88D0, &
    0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .68D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.9D0, &
    .33D0, 0.D0, 0.D0, 0.D0, 0.D0, .68D0, -.9D0, .33D0, .7D0, &
    -.75D0, .2D0, 1.04D0/
  !
  DATA ct8/(.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (.32,-1.41), (-1.55,.5), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (-1.55,.5), (.03,-.89), &
    (-.38,-.96), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (-.07,-.89), (-.9,.5), (.42,-1.41), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (.78,.06), (-.9,.5), (.06,-.13), (.1,-.5), (-.77,-.49)&
    , (-.5,-.3), (.52,-1.51), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.07,-.89), (-1.18,-.31), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.78,.06), &
    (-1.54,.97), (.03,-.89), (-.18,-1.31), (0.,0.), (0.,0.), (0.,0.)&
    , (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (.32,-1.41), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (.32,-1.41), (-.9,.5), (.05,-.6), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (.32,-1.41), (-.9,.5), (.05,-.6), &
    (.1,-.5), (-.77,-.49), (-.5,-.3), (.32,-1.16)/
  !
  !
  !                TRUE X VALUES AFTER ROTATION USING SROT OR DROT.
  DATA dt9x/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, -.46D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, .78D0, -.46D0, -.22D0, 1.06D0, 0.D0, 0.D0, &
    0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .66D0, .1D0, -.1D0, &
    0.D0, 0.D0, 0.D0, 0.D0, .96D0, .1D0, -.76D0, .8D0, .90D0, &
    -.3D0, -.02D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.06D0, .1D0, &
    -.1D0, 0.D0, 0.D0, 0.D0, 0.D0, .90D0, .1D0, -.22D0, .8D0, &
    .18D0, -.3D0, -.02D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .78D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, &
    .26D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .78D0, .26D0, -.76D0, &
    1.12D0, 0.D0, 0.D0, 0.D0/
  !
  !                TRUE Y VALUES AFTER ROTATION USING SROT OR DROT.
  !
  DATA dt9y/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.78D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, .04D0, -.78D0, .54D0, .08D0, 0.D0, 0.D0, &
    0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -.9D0, -.12D0, &
    0.D0, 0.D0, 0.D0, 0.D0, .64D0, -.9D0, -.30D0, .7D0, -.18D0, &
    .2D0, .28D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -1.08D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .64D0, -1.26D0, .54D0, .20D0, &
    0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .04D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.9D0, &
    .18D0, 0.D0, 0.D0, 0.D0, 0.D0, .04D0, -.9D0, .18D0, .7D0, &
    -.18D0, .2D0, .16D0/
  !
  DATA dt10x/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, -.9D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, &
    0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .3D0, .1D0, .5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .8D0, .1D0, -.6D0, .8D0, .3D0, -.3D0, .5D0, &
    .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.9D0, .1D0, .5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .7D0, .1D0, .3D0, .8D0, -.9D0, -.3D0, .5D0, &
    .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, .5D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .5D0, .3D0, -.6D0, .8D0, 0.D0, 0.D0, 0.D0/
  !
  DATA dt10y/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, .1D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, &
    .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.5D0, -.9D0, .6D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -.4D0, -.9D0, .9D0, .7D0, -.5D0, .2D0, .6D0, &
    .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.5D0, .6D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -.4D0, .9D0, -.5D0, .6D0, 0.D0, 0.D0, 0.D0, &
    .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, .6D0, -.9D0, .1D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, -.9D0, .1D0, .7D0, -.5D0, .2D0, .8D0/
  !
  DATA ct10x/(.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
    , (0.,0.), (.6,-.6), (-.9,.5), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (.6,-.6), (-.9,.5), (.7,-.6), (.1,-.5), &
    (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), (0.,0.)&
    , (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.6), (-.4,-.7), &
    (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.8,-.7), &
    (-.4,-.7), (-.1,-.2), (.2,-.8), (.7,-.6), (.1,.4), (.6,-.6), &
    (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
    , (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (-.9,.5), (-.4,-.7), (.6,-.6), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (.1,-.5), (-.4,-.7), (.7,-.6), (.2,-.8), &
    (-.9,.5), (.1,.4), (.6,-.6), (.7,-.8), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.)&
    , (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (.7,-.6), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (.7,-.6)&
    , (-.1,-.2), (.8,-.7), (0.,0.), (0.,0.), (0.,0.)/
  !
  DATA ct10y/(.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
    , (0.,0.), (.7,-.8), (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (.7,-.8), (-.4,-.7), (-.1,-.9), (.2,-.8), &
    (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.)&
    , (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.1,-.9), (-.9,.5), &
    (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (-.6,.6), &
    (-.9,.5), (-.9,-.4), (.1,-.5), (-.1,-.9), (-.5,-.3), (.7,-.8), &
    (.6,-.6), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.)&
    , (.7,-.8), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (-.1,-.9), (.7,-.8), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (-.6,.6), (-.9,-.4), (-.1,-.9), (.7,-.8), &
    (0.,0.), (0.,0.), (0.,0.), (.6,-.6), (0.,0.), (0.,0.), (0.,0.)&
    , (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (0.,0.), (0.,0.), &
    (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.8), (-.9,.5), &
    (-.4,-.7), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (.7,-.8), &
    (-.9,.5), (-.4,-.7), (.1,-.5), (-.1,-.9), (-.5,-.3), (.2,-.8)/
  !                        TRUE X RESULTS F0R ROTATIONS SROTM AND DROTM
  DATA dt19xa/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 3.8D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 2.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 3.5D0, -.4D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, -.8D0, &
    3.8D0, -2.2D0, -1.2D0, 0.D0, 0.D0, 0.D0, -.9D0, 2.8D0, &
    -1.4D0, -1.3D0, 0.D0, 0.D0, 0.D0, 3.5D0, -.4D0, -2.2D0, &
    4.7D0, 0.D0, 0.D0, 0.D0/
  !
  DATA dt19xb/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .1D0, &
    -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, -.3D0, .1D0, -2.0D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 3.3D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, -2.0D0, &
    .1D0, 1.4D0, .8D0, .6D0, -.3D0, -2.8D0, -1.8D0, .1D0, 1.3D0, &
    .8D0, 0.D0, -.3D0, -1.9D0, 3.8D0, .1D0, -3.1D0, .8D0, 4.8D0, &
    -.3D0, -1.5D0/
  !
  DATA dt19xc/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .6D0, .1D0, -.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.8D0, .1D0, &
    -3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 3.3D0, .1D0, -2.0D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 2.1D0, .1D0, -2.0D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .6D0, .1D0, -.5D0, .8D0, .9D0, -.3D0, -.4D0, -1.6D0, &
    .1D0, -2.2D0, .8D0, 5.4D0, -.3D0, -2.8D0, -1.5D0, .1D0, &
    -1.4D0, .8D0, 3.6D0, -.3D0, -1.9D0, 3.7D0, .1D0, -2.2D0, &
    .8D0, 3.6D0, -.3D0, -1.5D0/
  !
  DATA dt19xd/.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .6D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .6D0, .1D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.8D0, -1.0D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, -.9D0, -.8D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 3.5D0, .8D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, .6D0, .1D0, -.5D0, .8D0, 0.D0, 0.D0, 0.D0, -.8D0, &
    -1.0D0, 1.4D0, -1.6D0, 0.D0, 0.D0, 0.D0, -.9D0, -.8D0, &
    1.3D0, -1.6D0, 0.D0, 0.D0, 0.D0, 3.5D0, .8D0, -3.1D0, 4.8D0, &
    0.D0, 0.D0, 0.D0/
  !                        TRUE Y RESULTS FOR ROTATIONS SROTM AND DROTM
  DATA dt19ya/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -4.8D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.7D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -2.6D0, 3.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, .7D0, -4.8D0, &
    3.0D0, 1.1D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.7D0, -.7D0, 2.3D0, &
    0.D0, 0.D0, 0.D0, -2.6D0, 3.5D0, -.7D0, -3.6D0, 0.D0, 0.D0, &
    0.D0/
  !
  DATA dt19yb/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.0D0, -.9D0, -.3D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.5D0, -.9D0, 1.5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -1.5D0, -.9D0, -1.8D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, 3.7D0, -.9D0, &
    -1.2D0, .7D0, -1.5D0, .2D0, 2.2D0, -.3D0, -.9D0, 2.1D0, &
    .7D0, -1.6D0, .2D0, 2.0D0, -1.6D0, -.9D0, -2.1D0, .7D0, &
    2.9D0, .2D0, -3.8D0/
  !
  DATA dt19yc/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    -.9D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 4.0D0, -6.3D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, -.5D0, .3D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -1.5D0, 3.0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .5D0, -.9D0, .3D0, .7D0, 0.D0, 0.D0, 0.D0, 3.7D0, -7.2D0, &
    3.0D0, 1.7D0, 0.D0, 0.D0, 0.D0, -.3D0, .9D0, -.7D0, 1.9D0, &
    0.D0, 0.D0, 0.D0, -1.6D0, 2.7D0, -.7D0, -3.4D0, 0.D0, 0.D0, &
    0.D0/
  !
  DATA dt19yd/.5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, .5D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, 0.D0, 0.D0, 0.D0, &
    0.D0, 0.D0, 0.D0, 1.7D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    0.D0, -2.6D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, .5D0, &
    -.9D0, .3D0, 0.D0, 0.D0, 0.D0, 0.D0, .7D0, -.9D0, 1.2D0, &
    0.D0, 0.D0, 0.D0, 0.D0, 1.7D0, -.9D0, .5D0, 0.D0, 0.D0, &
    0.D0, 0.D0, -2.6D0, -.9D0, -1.3D0, 0.D0, 0.D0, 0.D0, 0.D0, &
    .5D0, -.9D0, .3D0, .7D0, -.6D0, .2D0, .8D0, .7D0, -.9D0, &
    1.2D0, .7D0, -1.5D0, .2D0, 1.6D0, 1.7D0, -.9D0, .5D0, .7D0, &
    -1.6D0, .2D0, 2.4D0, -2.6D0, -.9D0, -1.3D0, .7D0, 2.9D0, &
    .2D0, -4.0D0/
  !
  DATA ssize1/0., .3, 1.6, 3.2/
  DATA dsize1/0.D0, .3D0, 1.6D0, 3.2D0/
  DATA ssize3/.1, .4, 1.7, 3.3/
  !
  !                         FOR CDOTC AND CDOTU
  !
  DATA csize1/(0.,0.), (.9,.9), (1.63,1.73), (2.90,2.78)/
  DATA ssize2/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
    0., 0., 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, &
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17/
  DATA dsize2/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.17D0, &
    1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0, 1.17D0/
  !
  !                         FOR CAXPY
  !
  DATA csize2/(0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), (0.,0.), &
    (0.,0.), (1.54,1.54), (1.54,1.54), (1.54,1.54), (1.54,1.54), &
    (1.54,1.54), (1.54,1.54), (1.54,1.54)/
  !
  !                         FOR SROTM AND DROTM
  !
  DATA dpar/ - 2.D0, 0.D0, 0.D0, 0.D0, 0.D0, -1.D0, 2.D0, -3.D0, &
    -4.D0, 5.D0, 0.D0, 0.D0, 2.D0, -3.D0, 0.D0, 1.D0, 5.D0, &
    2.D0, 0.D0, -4.D0/
  !***FIRST EXECUTABLE STATEMENT  CHECK2
  DO ki = 1, 4
    INCx = incxs(ki)
    INCy = incys(ki)
    mx = ABS(INCx)
    my = ABS(INCy)
    !
    DO kn = 1, 4
      N = ns(kn)
      ksize = MIN(2,kn)
      lenx = lens(kn,mx)
      leny = lens(kn,my)
      !                                       INITIALIZE ALL ARGUMENT ARRAYS.
      DO i = 1, 7
        sx(i) = dx1(i)
        sy(i) = dy1(i)
        dx(i) = dx1(i)
        dy(i) = dy1(i)
        cx(i) = cx1(i)
        cy(i) = cy1(i)
      ENDDO
      !
      !                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
      !
      SELECT CASE (ICAse)
        CASE (2)
          !                                                              2. DSDOT
          CALL STEST(1,REAL(DSDOT(N,sx,INCx,sy,INCy)),REAL(dt7(kn,ki)),&
            ssize1(kn),Sfac,Kprint)
        CASE (3)
          !                                                              3. SDSDOT
          CALL STEST(1,SDSDOT(N,sb,sx,INCx,sy,INCy),st7b(kn,ki),ssize3(kn),&
            Sfac,Kprint)
        CASE (4)
          !                                                              4. DDOT
          CALL DTEST(1,DDOT(N,dx,INCx,dy,INCy),dt7(kn,ki),dsize1(kn),Dfac,&
            Kprint)
        CASE (5)
          !                                                              5. DQDOTI
          !                        DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
          !                        PRECISION ARITHMETIC INTERNALLY.
          !     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
          !     IN THE DIAGNOSTIC OUTPUT.
          !
          MODe = 1
          CALL DTEST(1,DQDOTI(N,db,qc,dx2,INCx,dy2,INCy),dt2(kn,ki,1),&
            dt2(kn,ki,1),Dqfac,Kprint)
        CASE (6)
          !                                                              6. DQDOTA
          !     TO TEST DQDOTA WE ACTUALLY TEST BOTH DQDOTI AND DQDOTA.
          !     THE OUTPUT VALUE OF QX FROM DQDOTI WILL BE USED AS INPUT
          !     TO DQDOTA.  QX IS SUPPOSED TO BE IN A MACHINE-DEPENDENT
          !     EXTENDED PRECISION FORM.
          !     MODE IS SET TO 1 OR 2 TO DISTINGUISH TESTS OF
          !     DQDOTI OR DQDOTA IN THE DIAGNOSTIC OUTPUT.
          !
          MODe = 1
          CALL DTEST(1,DQDOTI(N,db,qc,dx2,INCx,dy2,INCy),dt2(kn,ki,1),&
            dt2(kn,ki,1),Dqfac,Kprint)
          MODe = 2
          CALL DTEST(1,DQDOTA(N,-db,qc,dx2,INCx,dy2,INCy),dt2(kn,ki,2),&
            dt2(kn,ki,2),Dqfac,Kprint)
        CASE (7)
          !                                                              7. CDOTC
          CALL STEST(2,CDOTC(N,cx,INCx,cy,INCy),ct6(kn,ki),csize1(kn),Sfac,&
            Kprint)
        CASE (8)
          !                                                              8. CDOTU
          CALL STEST(2,CDOTU(N,cx,INCx,cy,INCy),ct7(kn,ki),csize1(kn),Sfac,&
            Kprint)
        CASE (9)
          !                                                              9. SAXPY
          CALL SAXPY(N,sa,sx,INCx,sy,INCy)
          DO j = 1, leny
            sty(j) = dt8(j,kn,ki)
          ENDDO
          CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac,Kprint)
        CASE (10)
          !                                                              10. DAXPY
          CALL DAXPY(N,da,dx,INCx,dy,INCy)
          CALL DTEST(leny,dy,dt8(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
        CASE (11)
          !                                                              11. CAXPY
          CALL CAXPY(N,ca,cx,INCx,cy,INCy)
          CALL STEST(2*leny,cy,ct8(1,kn,ki),csize2(1,ksize),Sfac,Kprint)
        CASE (12,13,16,17)
          GOTO 100
        CASE (14)
          !                                                              14. SROT
          DO i = 1, 7
            sx(i) = dx1(i)
            sy(i) = dy1(i)
            stx(i) = dt9x(i,kn,ki)
            sty(i) = dt9y(i,kn,ki)
          ENDDO
          CALL SROT(N,sx,INCx,sy,INCy,sc,ss)
          CALL STEST(lenx,sx,stx,ssize2(1,ksize),Sfac,Kprint)
          CALL STEST(leny,sy,sty,ssize2(1,ksize),Sfac,Kprint)
        CASE (15)
          !                                                             15. DROT
          DO i = 1, 7
            dx(i) = dx1(i)
            dy(i) = dy1(i)
          ENDDO
          CALL DROT(N,dx,INCx,dy,INCy,dc,ds)
          CALL DTEST(lenx,dx,dt9x(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
          CALL DTEST(leny,dy,dt9y(1,kn,ki),dsize2(1,ksize),Dfac,Kprint)
        CASE (18)
          !                                                             18. SROTM
          kni = kn + 4*(ki-1)
          DO kpar = 1, 4
            DO i = 1, 7
              sx(i) = dx1(i)
              sy(i) = dy1(i)
              stx(i) = dt19x(i,kpar,kni)
              sty(i) = dt19y(i,kpar,kni)
            ENDDO
            !
            DO i = 1, 5
              sparam(i) = dpar(i,kpar)
            ENDDO
            !                          SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
            !                          IF ANY
            MODe = INT(sparam(1))
            !
            DO i = 1, lenx
              ssize(i) = stx(i)
            ENDDO
            !                         THE TRUE RESULTS DT19X(1,2,7) AND
            !                         DT19X(5,3,8) ARE ZERO DUE TO CANCELLATION.
            !                         DT19X(1,2,7) = 2.*.6 - 4.*.3 = 0
            !                         DT19X(5,3,8) = .9 - 3.*.3 = 0
            !                         FOR THESE CASES RESPECTIVELY SET SIZE( )
            !                         EQUAL TO 2.4 AND 1.8
            IF ( (kpar==2).AND.(kni==7) ) ssize(1) = 2.4E0
            IF ( (kpar==3).AND.(kni==8) ) ssize(5) = 1.8E0
            !
            CALL SROTM(N,sx,INCx,sy,INCy,sparam)
            CALL STEST(lenx,sx,stx,ssize,Sfac,Kprint)
            CALL STEST(leny,sy,sty,sty,Sfac,Kprint)
          ENDDO
        CASE (19)
          !                                                             19. DROTM
          kni = kn + 4*(ki-1)
          DO kpar = 1, 4
            DO i = 1, 7
              dx(i) = dx1(i)
              dy(i) = dy1(i)
              dtx(i) = dt19x(i,kpar,kni)
              dty(i) = dt19y(i,kpar,kni)
            ENDDO
            !
            DO i = 1, 5
              dparam(i) = dpar(i,kpar)
            ENDDO
            !                            SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
            !                            IF ANY
            MODe = INT(dparam(1))
            !
            DO i = 1, lenx
              dsize(i) = dtx(i)
            ENDDO
            !                             SEE REMARK ABOVE ABOUT DT11X(1,2,7)
            !                             AND DT11X(5,3,8).
            IF ( (kpar==2).AND.(kni==7) ) dsize(1) = 2.4D0
            IF ( (kpar==3).AND.(kni==8) ) dsize(5) = 1.8D0
            !
            CALL DROTM(N,dx,INCx,dy,INCy,dparam)
            CALL DTEST(lenx,dx,dtx,dsize,Dfac,Kprint)
            CALL DTEST(leny,dy,dty,dty,Dfac,Kprint)
          ENDDO
        CASE (20)
          !                                                             20. SCOPY
          DO i = 1, 7
            sty(i) = dt10y(i,kn,ki)
          ENDDO
          CALL SCOPY(N,sx,INCx,sy,INCy)
          CALL STEST(leny,sy,sty,ssize2(1,1),1.,Kprint)
        CASE (21)
          !                                                             21. DCOPY
          CALL DCOPY(N,dx,INCx,dy,INCy)
          CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2(1,1),1.D0,Kprint)
        CASE (22)
          !                                                             22. CCOPY
          CALL CCOPY(N,cx,INCx,cy,INCy)
          CALL STEST(2*leny,cy,ct10y(1,kn,ki),ssize2(1,1),1.,Kprint)
        CASE (23)
          !                                                             23. SSWAP
          CALL SSWAP(N,sx,INCx,sy,INCy)
          DO i = 1, 7
            stx(i) = dt10x(i,kn,ki)
            sty(i) = dt10y(i,kn,ki)
          ENDDO
          CALL STEST(lenx,sx,stx,ssize2(1,1),1.,Kprint)
          CALL STEST(leny,sy,sty,ssize2(1,1),1.,Kprint)
        CASE (24)
          !                                                             24. DSWAP
          CALL DSWAP(N,dx,INCx,dy,INCy)
          CALL DTEST(lenx,dx,dt10x(1,kn,ki),dsize2(1,1),1.D0,Kprint)
          CALL DTEST(leny,dy,dt10y(1,kn,ki),dsize2(1,1),1.D0,Kprint)
        CASE (25)
          !                                                             25. CSWAP
          CALL CSWAP(N,cx,INCx,cy,INCy)
          CALL STEST(2*lenx,cx,ct10x(1,kn,ki),ssize2(1,1),1.,Kprint)
          CALL STEST(2*leny,cy,ct10y(1,kn,ki),ssize2(1,1),1.,Kprint)
        CASE DEFAULT
          !                                                              1. SDOT
          CALL STEST(1,SDOT(N,sx,INCx,sy,INCy),REAL(dt7(kn,ki)),ssize1(kn),&
            Sfac,Kprint)
      END SELECT
      !
      !
      !
    ENDDO
  ENDDO
  RETURN
  !                 THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  100  STOP
END SUBROUTINE CHECK2
