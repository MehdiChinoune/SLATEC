!** ERF
REAL FUNCTION ERF(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the error function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C8A, L5A1E
  !***
  ! **Type:**      SINGLE PRECISION (ERF-S, DERF-D)
  !***
  ! **Keywords:**  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ERF(X) calculates the single precision error function for
  ! single precision argument X.
  !
  ! Series for ERF        on the interval  0.          to  1.00000D+00
  !                                        with weighted error   7.10E-18
  !                                         log weighted error  17.15
  !                               significant figures required  16.31
  !                                    decimal places required  17.71
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, ERFC, INITS, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable name.  (RWC, WRB)

  REAL erfcs, sqeps, sqrtpi, X, xbig, y
  INTEGER nterf
  DIMENSION erfcs(13)
  LOGICAL first
  INTEGER, EXTERNAL ::  INITS
  REAL, EXTERNAL :: CSEVL, ERFC, R1MACH
  SAVE erfcs, sqrtpi, nterf, xbig, sqeps, first
  DATA erfcs(1)/ - .049046121234691808E0/
  DATA erfcs(2)/ - .14226120510371364E0/
  DATA erfcs(3)/.010035582187599796E0/
  DATA erfcs(4)/ - .000576876469976748E0/
  DATA erfcs(5)/.000027419931252196E0/
  DATA erfcs(6)/ - .000001104317550734E0/
  DATA erfcs(7)/.000000038488755420E0/
  DATA erfcs(8)/ - .000000001180858253E0/
  DATA erfcs(9)/.000000000032334215E0/
  DATA erfcs(10)/ - .000000000000799101E0/
  DATA erfcs(11)/.000000000000017990E0/
  DATA erfcs(12)/ - .000000000000000371E0/
  DATA erfcs(13)/.000000000000000007E0/
  DATA sqrtpi/1.7724538509055160E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  ERF
  IF ( first ) THEN
    nterf = INITS(erfcs,13,0.1*R1MACH(3))
    xbig = SQRT(-LOG(sqrtpi*R1MACH(3)))
    sqeps = SQRT(2.0*R1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1. ) THEN
    !
    ! ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1.
    !
    IF ( y<=xbig ) THEN
      ERF = SIGN(1.0-ERFC(y),X)
    ELSE
      ERF = SIGN(1.0,X)
    ENDIF
    RETURN
  ENDIF
  !
  ! ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1.
  !
  IF ( y<=sqeps ) ERF = 2.0*X/sqrtpi
  IF ( y>sqeps ) ERF = X*(1.0+CSEVL(2.*X**2-1.,erfcs,nterf))
  RETURN
END FUNCTION ERF
