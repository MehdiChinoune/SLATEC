!DECK FUNDOC
SUBROUTINE FUNDOC
  IMPLICIT NONE
  !***BEGIN PROLOGUE  FUNDOC
  !***PURPOSE  Documentation for FNLIB, a collection of routines for
  !            evaluating elementary and special functions.
  !***LIBRARY   SLATEC
  !***CATEGORY  C, Z
  !***TYPE      ALL (FUNDOC-A)
  !***KEYWORDS  DOCUMENTATION, ELEMENTARY FUNCTIONS, SPECIAL FUNCTIONS
  !***AUTHOR  Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  ! The SLATEC Library --  Elementary And Special Functions
  !
  ! This describes the elementary and special function routines available
  ! in the SLATEC library.  Most of the these routines were written by
  ! Wayne Fullerton while at LANL.  Some were written by Don Amos of SNLA.
  ! There are approximately 63 single precision, 63 double precision and
  ! 25 complex user callable elementary and special function routines.
  !
  ! The table below gives a breakdown of routines according to their
  ! function.  Unless otherwise indicated all routines are function
  ! subprograms.
  !                                             Sngl.      Dble.
  ! Description              Notation           Prec.      Prec.   Complex
  !
  !         ***Intrinsic Functions and Fundamental Functions***
  ! Unpack floating point              Call R9UPAK(X,Y,N)  D9UPAK    --
  !  number
  ! Pack floating point                        R9PAK(Y,N)  D9PAK     --
  !  number
  ! Initialize orthogonal               INITS(OS,NOS,ETA)  INITDS    --
  !  polynomial series
  ! Evaluate Chebyshev       summation for  CSEVL(X,CS,N)  DCSEVL    --
  ! series                  i = 1 to n of
  !                          cs(i)*(2*x)**(i-1)
  !
  !                  ***Elementary Functions***
  ! Argument = theta in      z = \ z \ *          --         --    CARG(Z)
  !  radians                 e**(i * theta)
  ! Cube root                                   CBRT(X)    DCBRT   CCBRT
  ! Relative error exponen-  ((e**x) -1) / x    EXPREL(X)  DEXPRL  CEXPRL
  !  tial from first order
  ! Common logarithm         log to the base 10   --         --  CLOG10(Z)
  !                          of z
  ! Relative error logarithm ln(1 + x)          ALNREL(X)  DLNREL  CLNREL
  ! Relative error logarithm (ln(1 + x) - x     R9LN2R(X)  D9LN2R  C9LN2R
  ! from second order        + x**2/2) / x**3
  !               ***Trigonometric and Hyperbolic Functions***
  ! Tangent                  tan z                --         --    CTAN(Z)
  ! Cotangent                cot x              COT(X)     DCOT    CCOT
  ! Sine x in degrees        sin((2*pi*x)/360)  SINDG(X)   DSINDG    --
  ! Cosine x in degrees      cos((2*pi*x)/360)  COSDG(X)   DCOSDG    --
  ! Arc sine                 arcsin (z)           --         --   CASIN(Z)
  ! Arc cosine               arccos (z)           --         --   CACOS(Z)
  ! Arc tangent              arctan (z)           --         --   CATAN(Z)
  ! Quadrant correct         arctan (z1/z2)       --         -- CATAN2(Z1,
  !  arc tangent                                                       Z2)
  ! Hyperbolic sine          sinh z               --         --   CSINH(Z)
  ! Hyperbolic cosine        cosh z               --         --   CCOSH(Z)
  ! Hyperbolic tangent       tanh z               --         --   CTANH(Z)
  ! Arc hyperbolic sine      arcsinh (x)        ASINH(X)   DASINH  CASINH
  ! Arc hyperbolic cosine    arccosh (x)        ACOSH(X)   DACOSH  CACOSH
  ! Arc hyperbolic tangent   arctanh (x)        ATANH(X)   DATANH  CATANH
  ! Relative error arc       (arctan (x) - x)   R9ATN1(X)  D9ATN1    --
  !  tangent from first order   / x**3
  !              ***Exponential Integrals and Related Functions***
  ! Exponential integral     Ei(x) = (minus)    EI(X)      DEI       --
  !                          the integral from
  !                          -x to infinity of
  !                            (e**-t / t)dt
  ! Exponential integral     E sub 1 (x) =      E1(X)      DE1       --
  !                          the integral from x
  !                            to infinity of
  !                          (e**-t / t) dt
  ! Logarithmic integral     li(x) = the        ALI(X)     DLI       --
  !                          integral from 0 to
  !                          x of (1 / ln t) dt
  !   Sequences of exponential integrals.
  !   M values are computed where
  !   k=0,1,...M-1 and n>=1
  ! Exponential integral     E sub n+k (x) Call EXINT(X,   DEXINT    --
  !                        =the integral from   N,KODE,M,TOL,
  !                         1 to infinity of    EN,IERR)
  !                       (e**(-x*t)/t**(n+k))dt
  !                 ***Gamma Functions and Related Functions***
  ! Factorial                n!                 FAC(N)     DFAC      --
  ! Binomial                 n!/(m!*(n-m)!)     BINOM(N,M) DBINOM    --
  ! Gamma                    gamma(x)           GAMMA(X)   DGAMMA  CGAMMA
  ! Gamma(x) under and                     Call GAMLIM(    DGAMLM    --
  !  overflow limits                           XMIN,XMAX)
  ! Reciprocal gamma         1 / gamma(x)       GAMR(X)    DGAMR   CGAMR
  ! Log abs gamma            ln \gamma(x)\      ALNGAM(X)  DLNGAM    --
  ! Log gamma                ln gamma(z)          --         --    CLNGAM
  ! Log abs gamma       g = ln \gamma(x)\  Call ALGAMS(X,  DLGAMS    --
  ! with sign           s = sign gamma(x)      G,S)
  ! Incomplete gamma         gamma(a,x) =       GAMI(A,X)  DGAMI     --
  !                          the integral from
  !                          0 to x of
  !                         (t**(a-1) * e**-t)dt
  ! Complementary            gamma(a,x) =       GAMIC(A,X) DGAMIC    --
  !  incomplete gamma        the integral from
  !                          x to infinity of
  !                         (t**(a-1) * e**-t)dt
  ! Tricomi's             gamma super star(a,x) GAMIT(A,X) DGAMIT    --
  !  incomplete gamma        = x**-a *
  !                         incomplete gamma(a,x)
  !                          / gamma(a)
  ! Psi (Digamma)            psi(x) = gamma'(x) PSI(X)     DPSI    CPSI
  !                          / gamma(x)
  ! Pochhammer's         (a) sub x = gamma(a+x) POCH(A,X)  DPOCH     --
  !  generalized symbol      / gamma(a)
  ! Pochhammer's symbol    ((a) sub x -1) / x   POCH1(A,X) DPOCH1    --
  !  from first order
  ! Beta                     b(a,b) = (gamma(a) BETA(A,B)  DBETA   CBETA
  !                          * gamma(b))
  !                          / gamma(a+b)
  !                           = the integral
  !                           from 0 to 1 of
  !                           (t**(a-1) *
  !                           (1-t)**(b-1))dt
  ! Log beta                 ln b(a,b)         ALBETA(A,B) DLBETA  CLBETA
  ! Incomplete beta          i sub x (a,b) =  BETAI(X,A,B) DBETAI    __
  !                          b sub x (a,b) / b(a,b)
  !                           = 1 / b(a,b) *
  !                          the integral
  !                          from 0 to x of
  !                          (t**(a-1) *
  !                          (1-t)**(b-1))dt
  ! Log gamma correction     ln gamma(x) -      R9LGMC(X)  D9LGMC  C9LGMC
  !  term when Stirling's    (ln(2 * pi))/2 -
  !  approximation is valid  (x - 1/2) * ln(x) + x
  !                ***Error Functions and Fresnel Integrals***
  ! Error function           erf x = (2 /       ERF(X)     DERF      --
  !                          square root of pi) *
  !                          the integral from
  !                          0 to x of
  !                          e**(-t**2)dt
  ! Complementary            erfc x = (2 /      ERFC(X)    DERFC     --
  !  error function          square root of pi) *
  !                          the integral from
  !                          x to infinity of
  !                          e**(-t**2)dt
  ! Dawson's function        F(x) = e**(-x**2)  DAWS(X)    DDAWS     --
  !                          * the integral from
  !                          from 0 to x of
  !                          e**(t**2)dt
  !                         ***Bessel Functions***
  !   Bessel functions of special integer order
  ! First kind, order zero   J sub 0 (x)        BESJ0(X)   DBESJ0    --
  ! First kind, order one    J sub 1 (x)        BESJ1(X)   DBESJ1    --
  ! Second kind, order zero  Y sub 0 (x)        BESY0(X)   DBESY0    --
  ! Second kind, order one   Y sub 1 (x)        BESY1(X)   DBESY1    --
  !   Modified (hyperbolic) Bessel functions of special integer order
  ! First kind, order zero   I sub 0 (x)        BESI0(X)   DBESI0    --
  ! First kind, order one    I sub 1 (x)        BESI1(X)   DBESI1    --
  ! Third kind, order zero   K sub 0 (x)        BESK0(X)   DBESK0    --
  ! Third kind, order one    K sub 1 (x)        BESK1(X)   DBESK1    --
  !   Modified (hyperbolic) Bessel functions of special integer order
  !   scaled by an exponential
  ! First kind, order zero   e**-\x\ * I sub 0(x) BESI0E(X) DBSI0E   --
  ! First kind, order one    e**-\x\ * I sub 1(x) BESI1E(X) DBSI1E   --
  ! Third kind, order zero   e**x * K sub 0 (x)   BESK0E(X) DBSK0E   --
  ! Third kind, order one    e**x * K sub 1 (x)   BESK1E(X) DBSK1E   --
  !   Sequences of Bessel functions of general order.
  !   N values are computed where  k = 1,2,...N and v .ge. 0.
  ! Modified first kind      I sub v+k-1 (x) Call BESI(X,   DBESI    --
  !                          optional scaling  ALPHA,KODE,N,
  !                          by e**(-x)        Y,NZ)
  ! First kind               J sub v+k-1 (x) Call BESJ(X,   DBESJ    --
  !                                            ALPHA,N,Y,NZ)
  ! Second kind              Y sub v+k-1 (x) Call BESY(X,   DBESY    --
  !                                            FNU,N,Y)
  ! Modified third kind      K sub v+k-1 (x) Call BESK(X,   DBESK    --
  !                          optional scaling  FNU,KODE,N,Y,
  !                          by e**(x)         NZ)
  !   Sequences of Bessel functions.  \N\ values are computed where
  !   I = 0, 1, 2, ..., N-1  for N > 0  or I = 0, -1, -2, ..., N+1
  !   for N < 0.
  ! Modified third kind      K sub v+i (x)   Call BESKS(    DBESKS   --
  !                                           XNU,X,N,BK)
  !   Sequences of Bessel functions scaled by an exponential.
  !   \N\ values are computed where  I = 0, 1, 2, ..., N-1
  !   for N > 0  or  I = 0, -1, -2, ..., N+1  for N < 0.
  ! Modified third kind      e**x *         Call BESKES(    DBSKES   --
  !                          K sub v+i (x)     XNU,X,N,BK)
  !                ***Bessel Functions of Fractional Order***
  !   Airy functions
  ! Airy                     Ai(x)              AI(X)      DAI       --
  ! Bairy                    Bi(x)              BI(X)      DBI       --
  !   Exponentially scaled Airy functions
  ! Airy                     Ai(x), x <= 0      AIE(X)     DAIE      --
  !                          exp(2/3 * x**(3/2))
  !                          * Ai(x), x >= 0
  ! Bairy                    Bi(x), x <= 0      BIE(X)     DBIE      --
  !                          exp(-2/3 * x**(3/2))
  !                          * Bi(x), x >= 0
  !                 ***Confluent Hypergeometric Functions***
  ! Confluent                U(a,b,x)           CHU(A,B,X) DCHU      --
  !  hypergeometric
  !                     ***Miscellaneous Functions***
  ! Spence                   s(x) = - the       SPENC(X)   DSPENC    --
  !  dilogarithm             integral from
  !                          0 to x of
  !                          ((ln \1-y\) / y)dy
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801015  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Routine name changed from FNLIBD to FUNDOC.  (WRB)
  !   900723  PURPOSE section revised.  (WRB)
  !***END PROLOGUE  FUNDOC
  !***FIRST EXECUTABLE STATEMENT  FUNDOC
END SUBROUTINE FUNDOC
