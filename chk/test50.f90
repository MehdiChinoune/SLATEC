MODULE TEST50_MOD
  IMPLICIT NONE

CONTAINS
  !DECK QXBLKT
  SUBROUTINE QXBLKT(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL am, an, bm, bn, cm, cn, deltas, deltat, ermax, err, hds, &
      hdt, s, t, tds, tdt, temp1, temp2, temp3, w
    REAL y, z
    INTEGER i, idimy, ierror, iflg, Ipass, j, Kprint, Lun, m, mp, &
      n, np
    !***BEGIN PROLOGUE  QXBLKT
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !     PROGRAM TO ILLUSTRATE THE USE OF BLKTRI
    !
    !***ROUTINES CALLED  BLKTRI
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXBLKT
    DIMENSION y(75,105), am(75), bm(75), cm(75), an(105), bn(105), &
      cn(105), w(1952), s(75), t(105)
    !***FIRST EXECUTABLE STATEMENT  QXBLKT
    ermax = 1.E-3
    iflg = 0
    np = 1
    n = 63
    mp = 1
    m = 50
    idimy = 75
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
    !     COEFFICIENTS AND THE ARRAY Y.
    !
    deltas = 1.0E0/(m+1)
    DO i = 1, m
      s(i) = i*deltas
    ENDDO
    deltat = 1.0E0/(n+1)
    DO j = 1, n
      t(j) = j*deltat
    ENDDO
    !
    !     COMPUTE THE COEFFICIENTS AM, BM AND CM CORRESPONDING TO THE S
    !     DIRECTION.
    !
    hds = deltas/2.
    tds = deltas + deltas
    DO i = 1, m
      temp1 = 1./(s(i)*tds)
      temp2 = 1./((s(i)-hds)*tds)
      temp3 = 1./((s(i)+hds)*tds)
      am(i) = temp1*temp2
      cm(i) = temp1*temp3
      bm(i) = -(am(i)+cm(i))
    ENDDO
    !
    !     COMPUTE THE COEFFICIENTS AN, BN AND CN CORRESPONDING TO THE T
    !     DIRECTION.
    !
    hdt = deltat/2.
    tdt = deltat + deltat
    DO j = 1, n
      temp1 = 1./(t(j)*tdt)
      temp2 = 1./((t(j)-hdt)*tdt)
      temp3 = 1./((t(j)+hdt)*tdt)
      an(j) = temp1*temp2
      cn(j) = temp1*temp3
      bn(j) = -(an(j)+cn(j))
    ENDDO
    !
    !     COMPUTE RIGHT SIDE OF EQUATION
    !
    DO j = 1, n
      DO i = 1, m
        y(i,j) = 3.75*s(i)*t(j)*(s(i)**4.+t(j)**4.)
      ENDDO
    ENDDO
    !
    !     INCLUDE NONHOMOGENEOUS BOUNDARY INTO RIGHT SIDE. NOTE THAT THE
    !     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
    !
    DO j = 1, n
      y(m,j) = y(m,j) - cm(m)*t(j)**5.
    ENDDO
    DO i = 1, m
      y(i,n) = y(i,n) - cn(n)*s(i)**5.
    ENDDO
    DO
      !
      CALL BLKTRI(iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y,ierror,w)
      iflg = iflg + 1
      IF ( iflg>1 ) THEN
        !
        !     COMPUTE DISCRETIZATION ERROR
        !
        err = 0.
        DO j = 1, n
          DO i = 1, m
            z = ABS(y(i,j)-(s(i)*t(j))**5.)
            IF ( z>err ) err = z
          ENDDO
        ENDDO
        !
        Ipass = 1
        IF ( err>ermax ) Ipass = 0
        IF ( Kprint==0 ) RETURN
        IF ( Kprint>=2.OR.Ipass==0 ) THEN
          WRITE (Lun,99001) ierror, err, INT(w(1))
          !
          99001       FORMAT ('1',20X,'SUBROUTINE BLKTRI EXAMPLE'///10X,&
            'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
            'IERROR = 0'/18X,'DISCRETIZATION ERROR = 1.6478E-05'/12X,&
            'REQUIRED LENGTH OF W ARRAY = 823'//10X,&
            'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
            'DISCRETIZATION ERROR =',1PE12.5/12X,&
            'REQUIRED LENGTH OF W ARRAY =',I4)
          IF ( Ipass==1 ) THEN
            WRITE (Lun,99002)
            99002         FORMAT (60X,'PASS'/)
          ELSE
            WRITE (Lun,99003)
            99003         FORMAT (60X,'FAIL'/)
          ENDIF
        ENDIF
        RETURN
      ENDIF
    ENDDO
  END SUBROUTINE QXBLKT
  !DECK QXCRT
  SUBROUTINE QXCRT(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL a, b, bda, bdb, bdc, bdd, c, d, dum, elmbda, ermax, err, &
      f, pertrb, pi, piby2, PIMACH, pisq, w, x
    REAL y, z
    INTEGER i, idimf, ierror, Ipass, j, Kprint, Lun, m, mbdcnd, mp1, &
      n, nbdcnd, np1
    !***BEGIN PROLOGUE  QXCRT
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCRT TO SOLVE
    !     THE EQUATION
    !
    !     (D/DX)(DU/DX) + (D/DY)(DU/DY) - 4*U
    !
    !     = (2 - (4 + PI**2/4)*X**2)*COS((Y+1)*PI/2)
    !
    !     WITH THE BOUNDARY CONDITIONS
    !     ON THE RECTANGLE 0 .LT. X .LT. 2, -1 .LT. Y .LT. 3 WITH THE
    !
    !     U(0,Y) = 0
    !                                          -1 .LE. Y .LE. 3
    !     (DU/DX)(2,Y) = 4*COS((Y+1)*PI/2)
    !
    !     AND WITH U PERIODIC IN Y.
    !          THE X-INTERVAL WILL BE DIVIDED INTO 40 PANELS AND THE
    !     Y-INTERVAL WILL BE DIVIDED INTO 80 PANELS.
    !
    !***ROUTINES CALLED  HWSCRT, PIMACH
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXCRT
    DIMENSION f(45,82), bdb(81), w(1200), x(41), y(81)
    !***FIRST EXECUTABLE STATEMENT  QXCRT
    !
    !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
    !     IS DIMENSIONED 6*(N+1) + 8*(M+1).
    !
    idimf = 45
    ermax = 1.E-3
    a = 0.
    b = 2.
    m = 40
    mbdcnd = 2
    c = -1.
    d = 3.
    n = 80
    nbdcnd = 0
    elmbda = -4.
    !
    !     AUXILIARY QUANTITIES.
    !
    pi = PIMACH(dum)
    piby2 = pi/2.
    pisq = pi**2
    mp1 = m + 1
    np1 = n + 1
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
    !     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
    !
    DO i = 1, mp1
      x(i) = (i-1)/20.0E0
    ENDDO
    DO j = 1, np1
      y(j) = -1.0E0 + (j-1)/20.0E0
    ENDDO
    !
    !     GENERATE BOUNDARY DATA.
    !
    DO j = 1, np1
      bdb(j) = 4.*COS((y(j)+1.)*piby2)
    ENDDO
    !
    !     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
    !
    DO j = 1, np1
      f(1,j) = 0.
    ENDDO
    !
    !     GENERATE RIGHT SIDE OF EQUATION.
    !
    DO i = 2, mp1
      DO j = 1, np1
        f(i,j) = (2.-(4.+pisq/4.)*x(i)**2)*COS((y(j)+1.)*piby2)
      ENDDO
    ENDDO
    CALL HWSCRT(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,idimf,&
      pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
    !                U(X,Y) = X**2*COS((Y+1)*PIBY2)
    !
    err = 0.
    DO i = 1, mp1
      DO j = 1, np1
        z = ABS(f(i,j)-x(i)**2*COS((y(j)+1.)*piby2))
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint==0 ) RETURN
    IF ( Kprint>=2.OR.Ipass==0 ) THEN
      WRITE (Lun,99001) ierror, err, INT(w(1))
      !
      99001   FORMAT ('1',20X,'SUBROUTINE HWSCRT EXAMPLE'///10X,&
        'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
        'IERROR = 0'/18X,'DISCRETIZATION ERROR = 5.36508E-04'/12X,&
        'REQUIRED LENGTH OF W ARRAY = 880'//10X,&
        'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
        'DISCRETIZATION ERROR =',1PE12.5/12X,&
        'REQUIRED LENGTH OF W ARRAY =',I4)
      IF ( Ipass==1 ) THEN
        WRITE (Lun,99002)
        99002     FORMAT (60X,'PASS'/)
      ELSE
        WRITE (Lun,99003)
        99003     FORMAT (60X,'FAIL'/)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE QXCRT
  !DECK QXCSP
  SUBROUTINE QXCSP(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL bdrf, bdrs, bdtf, bdts, ci4, dphi, dr, dtheta, dum, elmbda, &
      ermax, err, f, pertrb, pi, PIMACH, r, rf, rs, si
    REAL tf, theta, ts, w, z
    INTEGER i, idimf, ierror, intl, Ipass, j, Kprint, Lun, m, &
      mbdcnd, mp1, n, nbdcnd, np1
    !***BEGIN PROLOGUE  QXCSP
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !     PROGRAM TO ILLUSTRATE THE USE OF HWSCSP
    !
    !***ROUTINES CALLED  HWSCSP, PIMACH
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXCSP
    DIMENSION f(48,33), bdtf(33), w(1200), r(33), theta(48)
    !***FIRST EXECUTABLE STATEMENT  QXCSP
    !
    !     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  SINCE M=36, N=32,
    !     L=N THEREFORE K=5 AND W IS DIMENSIONED 2*(L+1)*(K-1) + 6*(M+N)
    !     + MAX(4*N,6*M) + 14 = 902.
    !
    ermax = 1.E-3
    pi = PIMACH(dum)
    intl = 0
    ts = 0.
    tf = pi/2.
    m = 36
    mbdcnd = 6
    rs = 0.
    rf = 1.
    n = 32
    nbdcnd = 5
    elmbda = 0.
    idimf = 48
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
    !     BOUNDARY DATA AND THE RIGHT SIDE OF THE EQUATION.
    !
    mp1 = m + 1
    dtheta = tf/m
    DO i = 1, mp1
      theta(i) = (i-1)*dtheta
    ENDDO
    np1 = n + 1
    dr = 1.0E0/n
    DO j = 1, np1
      r(j) = (j-1)*dr
    ENDDO
    !
    !     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
    !
    DO j = 1, np1
      bdtf(j) = 0.
    ENDDO
    !
    !     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
    !
    DO i = 1, mp1
      f(i,n+1) = COS(theta(i))**4
    ENDDO
    !
    !     COMPUTE RIGHT SIDE OF EQUATION
    !
    DO i = 1, mp1
      ci4 = 12.0E0*COS(theta(i))**2
      DO j = 1, n
        f(i,j) = ci4*r(j)**2
      ENDDO
    ENDDO
    !
    CALL HWSCSP(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,bdrf,elmbda,&
      f,idimf,pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR
    !
    err = 0.
    DO i = 1, mp1
      ci4 = COS(theta(i))**4
      DO j = 1, n
        z = ABS(f(i,j)-ci4*r(j)**4)
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint/=0 ) THEN
      IF ( Kprint>=2.OR.Ipass==0 ) THEN
        WRITE (Lun,99001) ierror, err, INT(w(1))
        !
        99001     FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 1'///10X,&
          'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
          'IERROR = 0'/18X,'DISCRETIZATION ERROR = 7.99842E-04'/12X,&
          'REQUIRED LENGTH OF W ARRAY = 775'//10X,&
          'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
          'DISCRETIZATION ERROR =',1PE12.5/12X,&
          'REQUIRED LENGTH OF W ARRAY =',I4)
        IF ( Ipass==1 ) THEN
          WRITE (Lun,99003)
        ELSE
          WRITE (Lun,99004)
        ENDIF
      ENDIF
    ENDIF
    !
    !     THE FOLLOWING PROGRAM ILLUSTRATES THE USE OF HWSCSP TO SOLVE
    !     A THREE DIMENSIONAL PROBLEM WHICH HAS LONGITUDINAL DEPENDENCE
    !
    mbdcnd = 2
    nbdcnd = 1
    dphi = pi/72.
    elmbda = -2.0E0*(1.0E0-COS(dphi))/dphi**2
    !
    !     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
    !
    DO i = 1, mp1
      f(i,n+1) = SIN(theta(i))
    ENDDO
    !
    !     COMPUTE RIGHT SIDE OF THE EQUATION
    !
    DO j = 1, n
      DO i = 1, mp1
        f(i,j) = 0.
      ENDDO
    ENDDO
    !
    CALL HWSCSP(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,bdrf,elmbda,&
      f,idimf,pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
    !
    err = 0.
    DO i = 1, mp1
      si = SIN(theta(i))
      DO j = 1, np1
        z = ABS(f(i,j)-r(j)*si)
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint==0 ) RETURN
    IF ( Kprint>=2.OR.Ipass==0 ) THEN
      WRITE (Lun,99002) ierror, err, INT(w(1))
      99002   FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 2'///10X,&
        'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
        'IERROR = 0'/18X,'DISCRETIZATION ERROR = 5.86824E-05'/12X,&
        'REQUIRED LENGTH OF W ARRAY = 775'//10X,&
        'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
        'DISCRETIZATION ERROR =',1PE12.5/12X,&
        'REQUIRED LENGTH OF W ARRAY =',I4)
      IF ( Ipass==1 ) THEN
        WRITE (Lun,99003)
      ELSE
        WRITE (Lun,99004)
      ENDIF
    ENDIF
    RETURN
    99003 FORMAT (60X,'PASS'/)
    99004 FORMAT (60X,'FAIL'/)
  END SUBROUTINE QXCSP
  !DECK QXCYL
  SUBROUTINE QXCYL(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL a, b, bda, bdb, bdc, bdd, c, d, elmbda, ermax, err, f, &
      pertrb, r, w, x, z
    INTEGER i, idimf, ierror, Ipass, j, Kprint, Lun, m, mbdcnd, mp1, &
      n, nbdcnd, np1
    !***BEGIN PROLOGUE  QXCYL
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSCYL TO SOLVE
    !     THE EQUATION
    !
    !     (1/R)(D/DR)(R*(DU/DR)) + (D/DZ)(DU/DZ)
    !
    !     = (2*R*Z)**2*(4*Z**2 + 3*R**2)
    !
    !     ON THE RECTANGLE 0 .LT. R .LT. 1, 0 .LT. Z .LT. 1 WITH THE
    !     BOUNDARY CONDITIONS
    !
    !     U(0,Z) UNSPECIFIED
    !                                            0 .LE. Z .LE. 1
    !     (DU/DR)(1,Z) = 4*Z**4
    !
    !     AND
    !
    !     (DU/DZ)(R,0) = 0
    !                                            0 .LE. R .LE. 1
    !     (DU/DZ)(R,1) = 4*R**4 .
    !
    !          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
    !     Z-INTERVAL WILL BE DIVIDED INTO 100 PANELS.
    !
    !***ROUTINES CALLED  HWSCYL
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !   930415  Test modified to use a 64 by 128 grid.  (WRB)
    !***END PROLOGUE  QXCYL
    DIMENSION f(65,129), bda(129), bdb(129), bdc(65), bdd(65), w(1400), &
      r(65), z(129)
    !***FIRST EXECUTABLE STATEMENT  QXCYL
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1',20X,'SUBROUTINE HWSCYL EXAMPLE'//)
    !
    !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
    !
    idimf = 65
    ermax = 1.0E-3
    a = 0.0
    b = 1.0
    m = 64
    mbdcnd = 6
    c = 0.0
    d = 1.0
    n = 128
    nbdcnd = 3
    elmbda = 0.0
    !
    !     AUXILIARY QUANTITIES.
    !
    mp1 = m + 1
    np1 = n + 1
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
    !     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
    !
    DO i = 1, mp1
      r(i) = (i-1)/64.0E0
    ENDDO
    DO j = 1, np1
      z(j) = (j-1)/128.0E0
    ENDDO
    !
    !     GENERATE BOUNDARY DATA.
    !
    DO j = 1, np1
      bdb(j) = 4.0*z(j)**4
    ENDDO
    DO i = 1, mp1
      bdc(i) = 0.0
      bdd(i) = 4.0*r(i)**4
    ENDDO
    !
    !     BDA IS A DUMMY VARIABLE.
    !
    !     GENERATE RIGHT SIDE OF EQUATION.
    !
    DO i = 1, mp1
      DO j = 1, np1
        f(i,j) = 4.0*r(i)**2*z(j)**2*(4.0*z(j)**2+3.0*r(i)**2)
      ENDDO
    ENDDO
    CALL HWSCYL(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,idimf,&
      pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
    !     NORM(F(I,J) - A - U(R(I),Z(J))).  THE EXACT SOLUTION IS
    !     U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
    !
    x = 0.0
    DO i = 1, mp1
      DO j = 1, np1
        x = x + f(i,j) - (r(i)*z(j))**4
      ENDDO
    ENDDO
    x = x/(np1*mp1)
    DO i = 1, mp1
      DO j = 1, np1
        f(i,j) = f(i,j) - x
      ENDDO
    ENDDO
    err = 0.0
    DO i = 1, mp1
      DO j = 1, np1
        x = ABS(f(i,j)-(r(i)*z(j))**4)
        IF ( x>err ) err = x
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint>=3.OR.(Kprint>=2.AND.Ipass==0) ) WRITE (Lun,99002) ierror, &
      pertrb, err, INT(w(1))
    99002 FORMAT (10X,'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/32X,&
      'PERTRB =',E12.5/18X,'DISCRETIZATION ERROR =',1PE12.5/12X,&
      'REQUIRED LENGTH OF W ARRAY = ',I4)
    !
    !     Print PASS/FAIL message.
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99003)
    99003 FORMAT (25X,'HWSCYL TEST PASSED'/)
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99004)
    99004 FORMAT (25X,'HWSCYL TEST FAILED'/)
    RETURN
  END SUBROUTINE QXCYL
  !DECK QXGBUN
  SUBROUTINE QXGBUN(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL a, b, c, deltax, deltay, dum, dysq, ermax, err, f, pi, &
      PIMACH, s, t, w, x, y, z
    INTEGER i, idimy, ierror, Ipass, j, Kprint, Lun, m, mm1, mperod, &
      n, nperod
    !***BEGIN PROLOGUE  QXGBUN
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE GENBUN
    !
    !***ROUTINES CALLED  GENBUN, PIMACH
    !***REVISION HISTORY  (YYMMDD)
    !   750701  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891009  Removed unreferenced variable.  (WRB)
    !   891009  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXGBUN
    DIMENSION f(25,130), a(20), b(20), c(20), w(1200), x(20), y(120)
    !***FIRST EXECUTABLE STATEMENT  QXGBUN
    !
    !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMY.  ALSO NOTE THAT
    !     W(.) IS DIMENSIONED 6*N + 5*M.
    !
    ermax = 1.E-2
    idimy = 25
    mperod = 1
    m = 20
    deltax = 1.0E0/m
    nperod = 0
    n = 120
    pi = PIMACH(dum)
    deltay = 2.0E0*pi/n
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
    !     COEFFICIENTS AND RIGHT SIDE OF EQUATION.
    !
    DO i = 1, m
      x(i) = (i-1)*deltax
    ENDDO
    DO j = 1, n
      y(j) = -pi + (j-1)*deltay
    ENDDO
    !
    !     GENERATE COEFFICIENTS.
    !
    s = (deltay/deltax)**2
    t = s*deltax
    a(1) = 0.
    b(1) = -2.0E0*s
    c(1) = 2.0E0*s
    DO i = 2, m
      a(i) = (1.+x(i))**2*s + (1.+x(i))*t
      c(i) = (1.+x(i))**2*s - (1.+x(i))*t
      b(i) = -2.0E0*(1.0E0+x(i))**2*s
    ENDDO
    c(m) = 0.
    !
    !     GENERATE RIGHT SIDE OF EQUATION FOR I = 1 SHOWING INTRODUCTION OF
    !     BOUNDARY DATA.
    !
    dysq = deltay**2
    DO j = 1, n
      f(1,j) = dysq*(11.+8./deltax)*SIN(y(j))
    ENDDO
    !
    !     GENERATE RIGHT SIDE.
    !
    mm1 = m - 1
    DO i = 2, mm1
      DO j = 1, n
        f(i,j) = dysq*3.*(1.+x(i))**4*SIN(y(j))
      ENDDO
    ENDDO
    !
    !     GENERATE RIGHT SIDE FOR I = M SHOWING INTRODUCTION OF
    !     BOUNDARY DATA.
    !
    DO j = 1, n
      f(m,j) = dysq*(3.*(1.+x(m))**4-16.*((1.+x(m))/deltax)**2+16.*(1.+x(m))&
        /deltax)*SIN(y(j))
    ENDDO
    CALL GENBUN(nperod,n,mperod,m,a,b,c,idimy,f,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
    !                   U(X,Y) = (1+X)**4*SIN(Y)
    !
    err = 0.
    DO i = 1, m
      DO j = 1, n
        z = ABS(f(i,j)-(1.+x(i))**4*SIN(y(j)))
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint==0 ) RETURN
    IF ( Kprint>=2.OR.Ipass==0 ) THEN
      WRITE (Lun,99001) ierror, err, INT(w(1))
      !
      99001   FORMAT ('1',20X,'SUBROUTINE GENBUN EXAMPLE'///10X,&
        'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
        'IERROR = 0'/18X,'DISCRETIZATION ERROR = 7.94113E-03'/12X,&
        'REQUIRED LENGTH OF W ARRAY = 740'//10X,&
        'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
        'DISCRETIZATION ERROR =',1PE12.5/12X,&
        'REQUIRED LENGTH OF W ARRAY =',I4)
      IF ( Ipass==1 ) THEN
        WRITE (Lun,99002)
        99002     FORMAT (60X,'PASS'/)
      ELSE
        WRITE (Lun,99003)
        99003     FORMAT (60X,'FAIL'/)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE QXGBUN
  !DECK QXPLR
  SUBROUTINE QXPLR(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL a, b, bda, bdb, bdc, bdd, c, d, dum, elmbda, ermax, err, &
      f, pertrb, pi, PIMACH, r, theta, w, z
    INTEGER i, idimf, ierror, Ipass, j, Kprint, Lun, m, mbdcnd, mp1, &
      n, nbdcnd, np1
    !***BEGIN PROLOGUE  QXPLR
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !          PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE HWSPLR TO SOLVE
    !     THE EQUATION
    !
    !     (1/R)(D/DR)(R*(DU/DR)) + (1/R**2)(D/DTHETA)(DU/DTHETA) = 16*R**2
    !
    !     ON THE QUARTER-DISK 0 .LT. R .LT. 1, 0 .LT. THETA .LT. PI/2 WITH
    !     WITH THE BOUNDARY CONDITIONS
    !
    !     U(1,THETA) = 1 - COS(4*THETA), 0 .LE. THETA .LE. 1
    !
    !     AND
    !
    !     (DU/DTHETA)(R,0) = (DU/DTHETA)(R,PI/2) = 0,  0 .LE. R .LE. 1.
    !
    !     (NOTE THAT THE SOLUTION U IS UNSPECIFIED AT R = 0.)
    !          THE R-INTERVAL WILL BE DIVIDED INTO 50 PANELS AND THE
    !     THETA-INTERVAL WILL BE DIVIDED INTO 48 PANELS.
    !
    !***ROUTINES CALLED  HWSPLR, PIMACH
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXPLR
    DIMENSION f(100,50), bdc(51), bdd(51), w(1200), r(51), theta(49)
    !***FIRST EXECUTABLE STATEMENT  QXPLR
    !
    !     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.  ALSO NOTE THAT W
    !     IS DIMENSIONED 6*(N+1) + 8*(M+1).
    !
    idimf = 100
    ermax = 1.E-3
    a = 0.
    b = 1.
    m = 50
    mbdcnd = 5
    c = 0.
    pi = PIMACH(dum)
    d = pi/2.
    n = 48
    nbdcnd = 3
    elmbda = 0.
    !
    !     AUXILIARY QUANTITIES.
    !
    mp1 = m + 1
    np1 = n + 1
    !
    !     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
    !     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
    !
    DO i = 1, mp1
      r(i) = (i-1)/50.0E0
    ENDDO
    DO j = 1, np1
      theta(j) = (j-1)*pi/96.0E0
    ENDDO
    !
    !     GENERATE BOUNDARY DATA.
    !
    DO i = 1, mp1
      bdc(i) = 0.
      bdd(i) = 0.
    ENDDO
    !
    !     BDA AND BDB ARE DUMMY VARIABLES.
    !
    DO j = 1, np1
      f(mp1,j) = 1. - COS(4.*theta(j))
    ENDDO
    !
    !     GENERATE RIGHT SIDE OF EQUATION.
    !
    DO i = 1, m
      DO j = 1, np1
        f(i,j) = 16.*r(i)**2
      ENDDO
    ENDDO
    CALL HWSPLR(a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,elmbda,f,idimf,&
      pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
    !                U(R,THETA) = R**4*(1 - COS(4*THETA))
    !
    err = 0.
    DO i = 1, mp1
      DO j = 1, np1
        z = ABS(f(i,j)-r(i)**4*(1.-COS(4.*theta(j))))
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint==0 ) RETURN
    IF ( Kprint>=2.OR.Ipass==0 ) THEN
      WRITE (Lun,99001) ierror, err, INT(w(1))
      !
      99001   FORMAT ('1',20X,'SUBROUTINE HWSPLR EXAMPLE'///10X,&
        'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
        'IERROR = 0'/18X,'DISCRETIZATION ERROR = 6.19134E-04'/12X,&
        'REQUIRED LENGTH OF W ARRAY = 882'//10X,&
        'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
        'DISCRETIZATION ERROR =',1PE12.5/12X,&
        'REQUIRED LENGTH OF W ARRAY =',I4)
      IF ( Ipass==1 ) THEN
        WRITE (Lun,99002)
        99002     FORMAT (60X,'PASS'/)
      ELSE
        WRITE (Lun,99003)
        99003     FORMAT (60X,'FAIL'/)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE QXPLR
  !DECK QXSSP
  SUBROUTINE QXSSP(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL bdpf, bdps, bdtf, bdts, dphi, dtheta, dum, elmbda, ermax, &
      err, f, pertrb, pf, pi, PIMACH, ps, sinp, sint, tf, ts
    REAL w, z
    INTEGER i, idimf, ierror, Ipass, j, Kprint, Lun, m, mbdcnd, mp1, &
      n, nbdcnd, np1
    !***BEGIN PROLOGUE  QXSSP
    !***PURPOSE
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  (UNKNOWN)
    !***DESCRIPTION
    !
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !     *                                                               *
    !     *                        F I S H P A K                          *
    !     *                                                               *
    !     *                                                               *
    !     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
    !     *                                                               *
    !     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
    !     *                                                               *
    !     *                  (VERSION  3, JUNE 1979)                     *
    !     *                                                               *
    !     *                             BY                                *
    !     *                                                               *
    !     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
    !     *                                                               *
    !     *                             OF                                *
    !     *                                                               *
    !     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
    !     *                                                               *
    !     *                BOULDER, COLORADO  (80307)  U.S.A.             *
    !     *                                                               *
    !     *                   WHICH IS SPONSORED BY                       *
    !     *                                                               *
    !     *              THE NATIONAL SCIENCE FOUNDATION                  *
    !     *                                                               *
    !     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !
    !
    !     PROGRAM TO ILLUSTRATE THE USE OF HWSSSP
    !
    !***ROUTINES CALLED  HWSSSP, PIMACH
    !***REVISION HISTORY  (YYMMDD)
    !   800103  DATE WRITTEN
    !   890718  Changed computation of PI to use PIMACH.  (WRB)
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   890911  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Added PASS/FAIL message and cleaned up FORMATs.  (RWC)
    !***END PROLOGUE  QXSSP
    DIMENSION f(19,73), bdtf(73), sint(19), sinp(73), w(1200)
    !***FIRST EXECUTABLE STATEMENT  QXSSP
    !
    !     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.  W IS
    !     DIMENSIONED 11*(M+1)+6*(N+1)=647 SINCE M=18 AND N=72.
    !
    pi = PIMACH(dum)
    ermax = 5.E-3
    ts = 0.0
    tf = pi/2.
    m = 18
    mbdcnd = 6
    ps = 0.0
    pf = pi + pi
    n = 72
    nbdcnd = 0
    elmbda = 0.
    idimf = 19
    !
    !     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
    !
    dtheta = tf/m
    mp1 = m + 1
    DO i = 1, mp1
      sint(i) = SIN((i-1)*dtheta)
    ENDDO
    dphi = (pi+pi)/n
    np1 = n + 1
    DO j = 1, np1
      sinp(j) = SIN((j-1)*dphi)
    ENDDO
    !
    !     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
    !
    DO j = 1, np1
      DO i = 1, mp1
        f(i,j) = 2. - 6.*(sint(i)*sinp(j))**2
      ENDDO
    ENDDO
    !
    !     STORE DERIVATIVE DATA AT THE EQUATOR
    !
    DO j = 1, np1
      bdtf(j) = 0.
    ENDDO
    !
    CALL HWSSSP(ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,bdpf,elmbda,f,&
      idimf,pertrb,ierror,w)
    !
    !     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
    !     SOLUTION MUST BE NORMALIZED.
    !
    err = 0.0
    DO j = 1, np1
      DO i = 1, mp1
        z = ABS(f(i,j)-(sint(i)*sinp(j))**2-f(1,1))
        IF ( z>err ) err = z
      ENDDO
    ENDDO
    !
    Ipass = 1
    IF ( err>ermax ) Ipass = 0
    IF ( Kprint==0 ) RETURN
    IF ( Kprint>=2.OR.Ipass==0 ) THEN
      WRITE (Lun,99001) ierror, err, INT(w(1))
      !
      99001   FORMAT ('1',20X,'SUBROUTINE HWSSSP EXAMPLE'///10X,&
        'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,&
        'IERROR = 0'/18X,'DISCRETIZATION ERROR = 3.38107E-03'/12X,&
        'REQUIRED LENGTH OF W ARRAY = 600'//10X,&
        'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,&
        'DISCRETIZATION ERROR =',1PE12.5/12X,&
        'REQUIRED LENGTH OF W ARRAY =',I4)
      IF ( Ipass==1 ) THEN
        WRITE (Lun,99002)
        99002     FORMAT (60X,'PASS'/)
      ELSE
        WRITE (Lun,99003)
        99003     FORMAT (60X,'FAIL'/)
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE QXSSP
END MODULE TEST50_MOD
!DECK TEST50
PROGRAM TEST50
  USE TEST50_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST50
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  I2
  !***TYPE      SINGLE PRECISION (TEST50-S)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        HWSCRT
  !        HWSPLR
  !        HWSCYL
  !        HWSSSP
  !        HWSCSP
  !        GENBUN
  !        BLKTRI
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, QXBLKT, QXCRT, QXCSP, QXCYL, QXGBUN, QXPLR,
  !                    QXSSP, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST50
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST50
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test HWSCRT
  !
  CALL QXCRT(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test HWSPLR
  !
  CALL QXPLR(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test HWSCYL
  !
  CALL QXCYL(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test HWSSSP
  !
  CALL QXSSP(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test HWSCSP
  !
  CALL QXCSP(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test GENBUN
  !
  CALL QXGBUN(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test BLKTRI
  !
  CALL QXBLKT(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST50 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST50 *************')
  ENDIF
  STOP
END PROGRAM TEST50
