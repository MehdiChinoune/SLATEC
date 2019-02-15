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
