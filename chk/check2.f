*DECK CHECK2
      SUBROUTINE CHECK2 (SFAC, SDFAC, DFAC, DQFAC, KPRINT)
C***BEGIN PROLOGUE  CHECK2
C***PURPOSE  (UNKNOWN)
C***LIBRARY   SLATEC
C***AUTHOR  Lawson, C. L., (JPL)
C***DESCRIPTION
C
C     THIS SUBPROGRAM TESTS THE BASIC LINEAR ALGEBRA SUBPROGRAMS 1-11,
C     14-15, AND 18-25. SUBPROGRAMS IN THIS SET EACH REQUIRE TWO ARRAYS
C     IN THE PARAMETER LIST.
C
C     C. L. LAWSON, JPL, 1975 FEB 26, APR 29, MAY 8, MAY 28
C
C***ROUTINES CALLED  CAXPY, CCOPY, CDOTC, CDOTU, CSWAP, DAXPY, DCOPY,
C                    DDOT, DQDOTA, DQDOTI, DROT, DROTM, DSDOT, DSWAP,
C                    DTEST, SAXPY, SCOPY, SDOT, SDSDOT, SROT, SROTM,
C                    SSWAP, STEST
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   750226  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CHECK2
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
C
      LOGICAL          PASS
      INTEGER          INCXS(4),INCYS(4),LENS(4,2),NS(4)
      REAL             SX(7),SY(7),STX(7),STY(7),SSIZE1(4),SSIZE2(14,2)
      REAL             SSIZE(7),QC(30),SPARAM(5),ST7B(4,4),SSIZE3(4)
      DOUBLE PRECISION DX(7),DA,DX1(7),DY1(7),DY(7),DT7(4,4),DT8(7,4,4)
      DOUBLE PRECISION DX2(7), DY2(7), DT2(4,4,2), DPARAM(5), DPAR(5,4)
      DOUBLE PRECISION DSDOT,DDOT,DQDOTI,DQDOTA,DFAC,DQFAC
      DOUBLE PRECISION DT10X(7,4,4),DT10Y(7,4,4),DB
      DOUBLE PRECISION DSIZE1(4),DSIZE2(7,2),DSIZE(7)
      DOUBLE PRECISION DC,DS,DT9X(7,4,4),DT9Y(7,4,4),DTX(7),DTY(7)
      DOUBLE PRECISION DT19X(7,4,16),DT19XA(7,4,4),DT19XB(7,4,4)
      DOUBLE PRECISION DT19XC(7,4,4),DT19XD(7,4,4),DT19Y(7,4,16)
      DOUBLE PRECISION DT19YA(7,4,4),DT19YB(7,4,4),DT19YC(7,4,4)
      DOUBLE PRECISION DT19YD(7,4,4)
C
      EQUIVALENCE (DT19X(1,1,1),DT19XA(1,1,1)),(DT19X(1,1,5),
     A   DT19XB(1,1,1)),(DT19X(1,1,9),DT19XC(1,1,1)),
     B   (DT19X(1,1,13),DT19XD(1,1,1))
      EQUIVALENCE (DT19Y(1,1,1),DT19YA(1,1,1)),(DT19Y(1,1,5),
     A   DT19YB(1,1,1)),(DT19Y(1,1,9),DT19YC(1,1,1)),
     B   (DT19Y(1,1,13),DT19YD(1,1,1))
      COMPLEX          CX(7),CA,CX1(7),CY1(7),CY(7),CT6(4,4),CT7(4,4)
      COMPLEX          CT8(7,4,4),CSIZE1(4),CSIZE2(7,2)
      COMPLEX          CT10X(7,4,4), CT10Y(7,4,4)
      COMPLEX          CDOTC,CDOTU
      DATA SA,DA,CA,DB,SB/.3,.3D0,(.4,-.7),.25D0,.1/
      DATA INCXS/   1,   2,  -2,  -1 /
      DATA INCYS/   1,  -2,   1,  -2 /
      DATA LENS/1, 1, 2, 4,   1, 1, 3, 7/
      DATA NS   /   0,   1,   2,   4 /
      DATA SC,SS,DC,DS/ .8,.6,.8D0,.6D0/
      DATA DX1/ .6D0, .1D0,-.5D0, .8D0, .9D0,-.3D0,-.4D0/
      DATA DY1/ .5D0,-.9D0, .3D0, .7D0,-.6D0, .2D0, .8D0/
      DATA DX2/ 1.D0,.01D0, .02D0,1.D0,.06D0, 2.D0, 1.D0/
      DATA DY2/ 1.D0,.04D0,-.03D0,-1.D0,.05D0,3.D0,-1.D0/
C            THE TERMS D11(3,2) AND D11(4,2) WILL BE SET BY
C            COMPUTATION AT RUN TIME.
      DATA CX1/(.7,-.8),(-.4,-.7),(-.1,-.9),(.2,-.8),(-.9,-.4),(.1,.4),
     *                                                        (-.6,.6)/
      DATA CY1/(.6,-.6),(-.9,.5),(.7,-.6),(.1,-.5),(-.1,-.2),(-.5,-.3),
     *                                                       (.8,-.7) /
C
C                             FOR DQDOTI AND DQDOTA
C
      DATA DT2/0.25D0,1.25D0,1.2504D0,0.2498D0,
     A         0.25D0,1.25D0,0.24D0,0.2492D0,
     B         0.25D0,1.25D0,0.31D0,0.2518D0,
     C         0.25D0,1.25D0,1.2497D0,0.2507D0,
     D         0.D0,2.D0,2.0008D0,-.0004D0,
     E         0.D0,2.D0,-.02D0,-.0016D0,
     F         0.D0,2.D0,.12D0,.0036D0,
     G         0.D0,2.D0,1.9994D0,.0014D0/
      DATA DT7/ 0.D0,.30D0,.21D0,.62D0,      0.D0,.30D0,-.07D0,.85D0,
     *          0.D0,.30D0,-.79D0,-.74D0,    0.D0,.30D0,.33D0,1.27D0/
      DATA ST7B/ .1, .4, .31, .72,     .1, .4, .03, .95,
     *           .1, .4, -.69, -.64,   .1, .4, .43, 1.37/
C
C                       FOR CDOTU
C
      DATA CT7/(0.,0.),(-.06,-.90),(.65,-.47),(-.34,-1.22),
     1         (0.,0.),(-.06,-.90),(-.59,-1.46),(-1.04,-.04),
     2         (0.,0.),(-.06,-.90),(-.83,.59),  (  .07,-.37),
     3         (0.,0.),(-.06,-.90),(-.76,-1.15),(-1.33,-1.82)/
C
C                       FOR CDOTC
C
      DATA CT6/(0.,0.),(.90,0.06), (.91,-.77),    (1.80,-.10),
     A         (0.,0.),(.90,0.06), (1.45,.74),    (.20,.90),
     B         (0.,0.),(.90,0.06), (-.55,.23),    (.83,-.39),
     C         (0.,0.),(.90,0.06), (1.04,0.79),    (1.95,1.22)/
C
      DATA DT8/.5D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     1         .68D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     2         .68D0,-.87D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     3         .68D0,-.87D0,.15D0,.94D0,          0.D0,0.D0,0.D0,
     4         .5D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     5         .68D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     6         .35D0,-.9D0,.48D0,                   0.D0,0.D0,0.D0,0.D0,
     7         .38D0,-.9D0,.57D0,.7D0,-.75D0,.2D0,.98D0,
     8         .5D0,                      0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     9         .68D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A         .35D0,-.72D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     B         .38D0,-.63D0,.15D0,.88D0,                 0.D0,0.D0,0.D0,
     C         .5D0,                      0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D         .68D0,                     0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E         .68D0,-.9D0,.33D0,                0.D0,0.D0,0.D0,0.D0,
     F         .68D0,-.9D0,.33D0,.7D0,-.75D0,.2D0,1.04D0/
C
      DATA CT8/
     A(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.32,-1.41),(-1.55,.5),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.32,-1.41),(-1.55,.5),(.03,-.89),(-.38,-.96),(0.,0.),(0.,0.),
     E                                                         (0.,0.),
     F(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(-.07,-.89),(-.9,.5),(.42,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     I(.78,.06),(-.9,.5),(.06,-.13),(.1,-.5),(-.77,-.49),(-.5,-.3),
     J                                                     (.52,-1.51),
     K(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     M(-.07,-.89),(-1.18,-.31),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.78,.06),(-1.54,.97),(.03,-.89),(-.18,-1.31),(0.,0.),(0.,0.),
     O(0.,0.),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.32,-1.41),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     Q(.32,-1.41),(-.9,.5),(.05,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     R(.32,-1.41),(-.9,.5),(.05,-.6),(.1,-.5),(-.77,-.49),(-.5,-.3),
     S                                                     (.32,-1.16) /
C
C
C                TRUE X VALUES AFTER ROTATION USING SROT OR DROT.
      DATA DT9X/.6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B          .78D0,-.46D0,               0.D0,0.D0,0.D0,0.D0,0.D0,
     C          .78D0,-.46D0,-.22D0,1.06D0,              0.D0,0.D0,0.D0,
     D          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F          .66D0,.1D0,-.1D0,                   0.D0,0.D0,0.D0,0.D0,
     G          .96D0,.1D0,-.76D0,.8D0,.90D0,-.3D0,-.02D0,
     H          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J          -.06D0,.1D0,-.1D0,                  0.D0,0.D0,0.D0,0.D0,
     K          .90D0,.1D0,-.22D0,.8D0,.18D0,-.3D0,-.02D0,
     L          .6D0,                    0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M          .78D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N          .78D0,.26D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     O          .78D0,.26D0,-.76D0,1.12D0,               0.D0,0.D0,0.D0/
C
C                TRUE Y VALUES AFTER ROTATION USING SROT OR DROT.
C
      DATA DT9Y/ .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .04D0,-.78D0,              0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .04D0,-.78D0, .54D0, .08D0,             0.D0,0.D0,0.D0,
     D           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           .7D0,-.9D0,-.12D0,                 0.D0,0.D0,0.D0,0.D0,
     G           .64D0,-.9D0,-.30D0, .7D0,-.18D0, .2D0, .28D0,
     H           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .04D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           .7D0,-1.08D0,              0.D0,0.D0,0.D0,0.D0,0.D0,
     K           .64D0,-1.26D0,.54D0, .20D0,             0.D0,0.D0,0.D0,
     L           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M          .04D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .04D0,-.9D0, .18D0,                0.D0,0.D0,0.D0,0.D0,
     O           .04D0,-.9D0, .18D0, .7D0,-.18D0, .2D0, .16D0/
C
      DATA DT10X/.6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .5D0,-.9D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .5D0,-.9D0,.3D0,.7D0,                   0.D0,0.D0,0.D0,
     D           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           .3D0,.1D0 ,.5D0,                   0.D0,0.D0,0.D0,0.D0,
     G           .8D0,.1D0 ,-.6D0,.8D0 ,.3D0,-.3D0,.5D0,
     H           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,.1D0,.5D0,                   0.D0,0.D0,0.D0,0.D0,
     K           .7D0, .1D0,.3D0, .8D0,-.9D0,-.3D0,.5D0,
     L           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .5D0,.3D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     O           .5D0,.3D0,-.6D0,.8D0,                   0.D0,0.D0,0.D0/
C
      DATA DT10Y/.5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B           .6D0,.1D0,                 0.D0,0.D0,0.D0,0.D0,0.D0,
     C           .6D0,.1D0,-.5D0,.8D0,                   0.D0,0.D0,0.D0,
     D           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.5D0,-.9D0,.6D0,                  0.D0,0.D0,0.D0,0.D0,
     G           -.4D0,-.9D0,.9D0, .7D0,-.5D0, .2D0,.6D0,
     H           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     I           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,.6D0,                0.D0,0.D0,0.D0,0.D0,0.D0,
     K           -.4D0,.9D0,-.5D0,.6D0,                  0.D0,0.D0,0.D0,
     L           .5D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     M           .6D0,                   0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     N           .6D0,-.9D0,.1D0,                   0.D0,0.D0,0.D0,0.D0,
     O           .6D0,-.9D0,.1D0, .7D0,-.5D0, .2D0, .8D0/
C
      DATA CT10X/
     A(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.6,-.6),(-.9,.5),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.6,-.6),(-.9,.5),(.7,-.6),(.1,-.5),(0.,0.),(0.,0.),(0.,0.),
     E(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     F(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(.7,-.6),(-.4,-.7),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(.8,-.7),(-.4,-.7),(-.1,-.2),(.2,-.8),(.7,-.6),(.1,.4),(.6,-.6),
     I(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     J(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     K(-.9,.5),(-.4,-.7),(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(.1,-.5),(-.4,-.7),(.7,-.6),(.2,-.8),(-.9,.5),(.1,.4),(.6,-.6),
     M(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     O(.6,-.6),(.7,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.6,-.6),(.7,-.6),(-.1,-.2),(.8,-.7),(0.,0.),(0.,0.),(0.,0.)   /
C
      DATA CT10Y/
     A(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     C(.7,-.8),(-.4,-.7),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     D(.7,-.8),(-.4,-.7),(-.1,-.9),(.2,-.8),(0.,0.),(0.,0.),(0.,0.),
     E(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     F(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     G(-.1,-.9),(-.9,.5),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     H(-.6,.6),(-.9,.5),(-.9,-.4),(.1,-.5),(-.1,-.9),(-.5,-.3),(.7,-.8),
     I(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     J(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     K(-.1,-.9),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     L(-.6,.6),(-.9,-.4),(-.1,-.9),(.7,-.8),(0.,0.),(0.,0.),(0.,0.),
     M(.6,-.6),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     N(.7,-.8),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     O(.7,-.8),(-.9,.5),(-.4,-.7),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     P(.7,-.8),(-.9,.5),(-.4,-.7),(.1,-.5),(-.1,-.9),(-.5,-.3),(.2,-.8)/
C                        TRUE X RESULTS F0R ROTATIONS SROTM AND DROTM
      DATA DT19XA/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           -.8D0,  3.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,  2.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K           3.5D0,  -.4D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0,
     M           -.8D0,  3.8D0, -2.2D0, -1.2D0,          0.D0,0.D0,0.D0,
     N           -.9D0,  2.8D0, -1.4D0, -1.3D0,          0.D0,0.D0,0.D0,
     O           3.5D0,  -.4D0, -2.2D0,  4.7D0,          0.D0,0.D0,0.D0/
C
      DATA DT19XB/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0,
     I           0.D0,    .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0,
     J           -.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     K           3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0,
     M          -2.0D0,   .1D0,  1.4D0,   .8D0,   .6D0,  -.3D0, -2.8D0,
     N          -1.8D0,   .1D0,  1.3D0,   .8D0,  0.D0,   -.3D0, -1.9D0,
     O           3.8D0,   .1D0, -3.1D0,   .8D0,  4.8D0,  -.3D0, -1.5D0 /
C
      DATA DT19XC/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,  -.5D0,             0.D0,0.D0,0.D0,0.D0,
     I           4.8D0,   .1D0, -3.0D0,             0.D0,0.D0,0.D0,0.D0,
     J           3.3D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     K           2.1D0,   .1D0, -2.0D0,             0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,   .9D0,  -.3D0,  -.4D0,
     M          -1.6D0,   .1D0, -2.2D0,   .8D0,  5.4D0,  -.3D0, -2.8D0,
     N          -1.5D0,   .1D0, -1.4D0,   .8D0,  3.6D0,  -.3D0, -1.9D0,
     O           3.7D0,   .1D0, -2.2D0,   .8D0,  3.6D0,  -.3D0, -1.5D0 /
C
      DATA DT19XD/.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E           -.8D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           -.9D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G           3.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .6D0,   .1D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           -.8D0, -1.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.9D0,  -.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K           3.5D0,   .8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .6D0,   .1D0,  -.5D0,   .8D0,          0.D0,0.D0,0.D0,
     M           -.8D0, -1.0D0,  1.4D0, -1.6D0,          0.D0,0.D0,0.D0,
     N           -.9D0,  -.8D0,  1.3D0, -1.6D0,          0.D0,0.D0,0.D0,
     O           3.5D0,   .8D0, -3.1D0,  4.8D0,          0.D0,0.D0,0.D0/
C                        TRUE Y RESULTS FOR ROTATIONS SROTM AND DROTM
      DATA DT19YA/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I            .7D0, -4.8D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           1.7D0,  -.7D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K          -2.6D0,  3.5D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0,
     M            .7D0, -4.8D0,  3.0D0,  1.1D0,          0.D0,0.D0,0.D0,
     N           1.7D0,  -.7D0,  -.7D0,  2.3D0,          0.D0,0.D0,0.D0,
     O          -2.6D0,  3.5D0,  -.7D0, -3.6D0,          0.D0,0.D0,0.D0/
C
      DATA DT19YB/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,
     I           4.0D0,  -.9D0,  -.3D0,             0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,  -.9D0,  1.5D0,             0.D0,0.D0,0.D0,0.D0,
     K          -1.5D0,  -.9D0, -1.8D0,             0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0,
     M           3.7D0,  -.9D0, -1.2D0,   .7D0, -1.5D0,   .2D0,  2.2D0,
     N           -.3D0,  -.9D0,  2.1D0,   .7D0, -1.6D0,   .2D0,  2.0D0,
     O          -1.6D0,  -.9D0, -2.1D0,   .7D0,  2.9D0,   .2D0, -3.8D0 /
C
      DATA DT19YC/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     I           4.0D0, -6.3D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     J           -.5D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     K          -1.5D0,  3.0D0,             0.D0,0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,          0.D0,0.D0,0.D0,
     M           3.7D0, -7.2D0,  3.0D0,  1.7D0,          0.D0,0.D0,0.D0,
     N           -.3D0,   .9D0,  -.7D0,  1.9D0,          0.D0,0.D0,0.D0,
     O          -1.6D0,  2.7D0,  -.7D0, -3.4D0,          0.D0,0.D0,0.D0/
C
      DATA DT19YD/.5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     B            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     C            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     D            .5D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     E            .7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     F           1.7D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     G          -2.6D0,                  0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     H            .5D0,  -.9D0,   .3D0,             0.D0,0.D0,0.D0,0.D0,
     I            .7D0,  -.9D0,  1.2D0,             0.D0,0.D0,0.D0,0.D0,
     J           1.7D0,  -.9D0,   .5D0,             0.D0,0.D0,0.D0,0.D0,
     K          -2.6D0,  -.9D0, -1.3D0,             0.D0,0.D0,0.D0,0.D0,
     L            .5D0,  -.9D0,   .3D0,   .7D0,  -.6D0,   .2D0,   .8D0,
     M            .7D0,  -.9D0,  1.2D0,   .7D0, -1.5D0,   .2D0,  1.6D0,
     N           1.7D0,  -.9D0,   .5D0,   .7D0, -1.6D0,   .2D0,  2.4D0,
     O          -2.6D0,  -.9D0, -1.3D0,   .7D0,  2.9D0,   .2D0, -4.0D0 /
C
      DATA SSIZE1/ 0.  , .3  , 1.6  , 3.2   /
      DATA DSIZE1/ 0.D0, .3D0, 1.6D0, 3.2D0 /
      DATA SSIZE3/ .1, .4, 1.7, 3.3 /
C
C                         FOR CDOTC AND CDOTU
C
      DATA CSIZE1/ (0.,0.), (.9,.9), (1.63,1.73), (2.90,2.78) /
      DATA SSIZE2/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     A  1.17,1.17,1.17,1.17,1.17,1.17,1.17,
     B  1.17,1.17,1.17,1.17,1.17,1.17,1.17/
      DATA DSIZE2/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,
     A  1.17D0,1.17D0,1.17D0,1.17D0,1.17D0,1.17D0,1.17D0/
C
C                         FOR CAXPY
C
      DATA CSIZE2/
     A (0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),
     B (1.54,1.54),(1.54,1.54),(1.54,1.54),(1.54,1.54),(1.54,1.54),
     C                                     (1.54,1.54),(1.54,1.54) /
C
C                         FOR SROTM AND DROTM
C
      DATA DPAR/-2.D0,  0.D0,0.D0,0.D0,0.D0,
     A          -1.D0,  2.D0, -3.D0, -4.D0,  5.D0,
     B           0.D0,  0.D0,  2.D0, -3.D0,  0.D0,
     C           1.D0,  5.D0,  2.D0,  0.D0, -4.D0/
C***FIRST EXECUTABLE STATEMENT  CHECK2
        DO 520 KI = 1, 4
        INCX = INCXS(KI)
        INCY = INCYS(KI)
        MX   = ABS(INCX)
        MY   = ABS(INCY)
C
          DO 500 KN=1,4
          N= NS(KN)
          KSIZE=MIN(2,KN)
          LENX = LENS(KN,MX)
          LENY = LENS(KN,MY)
C                                       INITIALIZE ALL ARGUMENT ARRAYS.
               DO 5 I = 1, 7
               SX(I) = DX1(I)
               SY(I) = DY1(I)
               DX(I) = DX1(I)
               DY(I) = DY1(I)
               CX(I) = CX1(I)
    5          CY(I) = CY1(I)
C
C                             BRANCH TO SELECT SUBPROGRAM TO BE TESTED.
C
          GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90,100,
     A           110,999,999,140,150,999,999,180,190,200,
     B           210,220,230,240,250), ICASE
C                                                              1. SDOT
   10     CALL STEST(1,SDOT(N,SX,INCX,SY,INCY),REAL(DT7(KN,KI)),
     *                                         SSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              2. DSDOT
   20     CALL STEST(1,REAL(DSDOT(N,SX,INCX,SY,INCY)),
     *               REAL(DT7(KN,KI)),SSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              3. SDSDOT
   30     CALL STEST(1,SDSDOT(N,SB,SX,INCX,SY,INCY),
     *               ST7B(KN,KI),SSIZE3(KN),SFAC,KPRINT)
          GO TO 500
C                                                              4. DDOT
   40     CALL DTEST(1,DDOT(N,DX,INCX,DY,INCY),DT7(KN,KI),
     *               DSIZE1(KN),DFAC,KPRINT)
          GO TO 500
C                                                              5. DQDOTI
   50 CONTINUE
C                        DQDOTI AND DQDOTA ARE SUPPOSED TO USE EXTENDED
C                        PRECISION ARITHMETIC INTERNALLY.
C     SET MODE = 1 OR 2 TO DISTINGUISH TESTS OF DQDOTI OR DQDOTA
C     IN THE DIAGNOSTIC OUTPUT.
C
          MODE = 1
          CALL DTEST(1,DQDOTI(N,DB,QC,DX2,INCX,DY2,INCY),
     *               DT2(KN,KI,1),DT2(KN,KI,1),DQFAC,KPRINT)
      GO TO 500
C                                                              6. DQDOTA
   60 CONTINUE
C     TO TEST DQDOTA WE ACTUALLY TEST BOTH DQDOTI AND DQDOTA.
C     THE OUTPUT VALUE OF QX FROM DQDOTI WILL BE USED AS INPUT
C     TO DQDOTA.  QX IS SUPPOSED TO BE IN A MACHINE-DEPENDENT
C     EXTENDED PRECISION FORM.
C     MODE IS SET TO 1 OR 2 TO DISTINGUISH TESTS OF
C     DQDOTI OR DQDOTA IN THE DIAGNOSTIC OUTPUT.
C
          MODE = 1
          CALL DTEST(1,DQDOTI(N,DB,QC,DX2,INCX,DY2,INCY),
     *               DT2(KN,KI,1),DT2(KN,KI,1),DQFAC,KPRINT)
          MODE = 2
          CALL DTEST(1,DQDOTA(N,-DB,QC,DX2,INCX,DY2,INCY),
     *               DT2(KN,KI,2),DT2(KN,KI,2),DQFAC,KPRINT)
          GO TO 500
C                                                              7. CDOTC
   70     CALL STEST(2, CDOTC(N,CX,INCX,CY,INCY),
     *               CT6(KN,KI),CSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              8. CDOTU
   80     CALL STEST(2,CDOTU(N,CX,INCX,CY,INCY),
     *               CT7(KN,KI),CSIZE1(KN),SFAC,KPRINT)
          GO TO 500
C                                                              9. SAXPY
   90     CALL SAXPY(N,SA,SX,INCX,SY,INCY)
               DO 95 J = 1, LENY
   95          STY(J) = DT8(J,KN,KI)
          CALL STEST(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                              10. DAXPY
  100      CALL DAXPY(N,DA,DX,INCX,DY,INCY)
          CALL DTEST(LENY,DY,DT8(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          GO TO 500
C                                                              11. CAXPY
  110     CALL CAXPY(N,CA,CX,INCX,CY,INCY)
          CALL STEST(2*LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                              14. SROT
  140     CONTINUE
               DO 144 I = 1, 7
               SX(I) = DX1(I)
               SY(I) = DY1(I)
               STX(I) = DT9X(I,KN,KI)
               STY(I) = DT9Y(I,KN,KI)
  144         CONTINUE
          CALL SROT   (N,SX,INCX,SY,INCY,SC,SS)
          CALL STEST(LENX,SX,STX,SSIZE2(1,KSIZE),SFAC,KPRINT)
          CALL STEST(LENY,SY,STY,SSIZE2(1,KSIZE),SFAC,KPRINT)
          GO TO 500
C                                                             15. DROT
  150     CONTINUE
               DO 154 I = 1, 7
               DX(I) = DX1(I)
               DY(I) = DY1(I)
  154          CONTINUE
          CALL DROT   (N,DX,INCX,DY,INCY,DC,DS)
          CALL DTEST(LENX,DX,DT9X(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          CALL DTEST(LENY,DY,DT9Y(1,KN,KI),DSIZE2(1,KSIZE),DFAC,KPRINT)
          GO TO 500
C                                                             18. SROTM
  180     KNI = KN + 4*(KI-1)
          DO 189 KPAR=1,4
          DO 182 I = 1, 7
          SX(I) = DX1(I)
          SY(I) = DY1(I)
          STX(I) = DT19X(I,KPAR,KNI)
  182     STY(I) = DT19Y(I,KPAR,KNI)
C
          DO 186 I = 1, 5
  186     SPARAM(I) = DPAR(I,KPAR)
C                          SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
C                          IF ANY
          MODE = INT(SPARAM(1))
C
          DO 187 I = 1, LENX
  187     SSIZE(I) = STX(I)
C                         THE TRUE RESULTS DT19X(1,2,7) AND
C                         DT19X(5,3,8) ARE ZERO DUE TO CANCELLATION.
C                         DT19X(1,2,7) = 2.*.6 - 4.*.3 = 0
C                         DT19X(5,3,8) = .9 - 3.*.3 = 0
C                         FOR THESE CASES RESPECTIVELY SET SIZE( )
C                         EQUAL TO 2.4 AND 1.8
          IF ((KPAR .EQ. 2) .AND. (KNI .EQ. 7))
     1           SSIZE(1) = 2.4E0
          IF ((KPAR .EQ. 3) .AND. (KNI .EQ. 8))
     1           SSIZE(5) = 1.8E0
C
          CALL SROTM(N,SX,INCX,SY,INCY,SPARAM)
          CALL STEST(LENX,SX,STX,SSIZE,SFAC,KPRINT)
          CALL STEST(LENY,SY,STY,STY,SFAC,KPRINT)
  189     CONTINUE
          GO TO 500
C                                                             19. DROTM
  190     KNI = KN + 4*(KI-1)
          DO 199 KPAR=1,4
            DO 192 I = 1, 7
            DX(I) = DX1(I)
            DY(I) = DY1(I)
            DTX(I) = DT19X(I,KPAR,KNI)
  192       DTY(I) = DT19Y(I,KPAR,KNI)
C
            DO 196 I = 1, 5
  196       DPARAM(I) = DPAR(I,KPAR)
C                            SET MODE TO IDENTIFY DIAGNOSTIC OUTPUT,
C                            IF ANY
          MODE = INT(DPARAM(1))
C
            DO 197 I = 1, LENX
  197       DSIZE(I) = DTX(I)
C                             SEE REMARK ABOVE ABOUT DT11X(1,2,7)
C                             AND DT11X(5,3,8).
          IF ((KPAR .EQ. 2) .AND. (KNI .EQ. 7))
     1               DSIZE(1) = 2.4D0
          IF ((KPAR .EQ. 3) .AND. (KNI .EQ. 8))
     1               DSIZE(5) = 1.8D0
C
          CALL   DROTM(N,DX,INCX,DY,INCY,DPARAM)
          CALL DTEST(LENX,DX,DTX,DSIZE,DFAC,KPRINT)
          CALL DTEST(LENY,DY,DTY,DTY,DFAC,KPRINT)
  199     CONTINUE
          GO TO 500
C                                                             20. SCOPY
  200     DO 205 I = 1, 7
  205     STY(I) = DT10Y(I,KN,KI)
          CALL SCOPY(N,SX,INCX,SY,INCY)
          CALL STEST(LENY,SY,STY,SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             21. DCOPY
  210     CALL DCOPY(N,DX,INCX,DY,INCY)
          CALL DTEST(LENY,DY,DT10Y(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          GO TO 500
C                                                             22. CCOPY
  220     CALL CCOPY(N,CX,INCX,CY,INCY)
          CALL STEST(2*LENY,CY,CT10Y(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             23. SSWAP
  230     CALL SSWAP(N,SX,INCX,SY,INCY)
               DO 235 I = 1, 7
               STX(I) = DT10X(I,KN,KI)
  235          STY(I) = DT10Y(I,KN,KI)
          CALL STEST(LENX,SX,STX,SSIZE2(1,1),1.,KPRINT)
          CALL STEST(LENY,SY,STY,SSIZE2(1,1),1.,KPRINT)
          GO TO 500
C                                                             24. DSWAP
  240     CALL DSWAP(N,DX,INCX,DY,INCY)
          CALL DTEST(LENX,DX,DT10X(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          CALL DTEST(LENY,DY,DT10Y(1,KN,KI),DSIZE2(1,1),1.D0,KPRINT)
          GO TO 500
C                                                             25. CSWAP
  250     CALL CSWAP(N,CX,INCX,CY,INCY)
          CALL STEST(2*LENX,CX,CT10X(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
          CALL STEST(2*LENY,CY,CT10Y(1,KN,KI),SSIZE2(1,1),1.,KPRINT)
C
C
C
  500     CONTINUE
  520   CONTINUE
      RETURN
C                 THE FOLLOWING STOP SHOULD NEVER BE REACHED.
  999 STOP
      END
