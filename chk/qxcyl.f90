!*==QXCYL.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXCYL
SUBROUTINE QXCYL(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXCYL5
  !*** Start of declarations inserted by SPAG
  REAL a , b , bda , bdb , bdc , bdd , c , d , elmbda , ermax , err , f , &
    pertrb , r , w , x , z
  INTEGER i , idimf , ierror , Ipass , j , Kprint , Lun , m , mbdcnd , mp1 , &
    n , nbdcnd , np1
  !*** End of declarations inserted by SPAG
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
  !     *                  (VERSION  3 , JUNE 1979)                     *
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
  DIMENSION f(65,129) , bda(129) , bdb(129) , bdc(65) , bdd(65) , w(1400) , &
    r(65) , z(129)
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
  DO i = 1 , mp1
    r(i) = (i-1)/64.0E0
  ENDDO
  DO j = 1 , np1
    z(j) = (j-1)/128.0E0
  ENDDO
  !
  !     GENERATE BOUNDARY DATA.
  !
  DO j = 1 , np1
    bdb(j) = 4.0*z(j)**4
  ENDDO
  DO i = 1 , mp1
    bdc(i) = 0.0
    bdd(i) = 4.0*r(i)**4
  ENDDO
  !
  !     BDA IS A DUMMY VARIABLE.
  !
  !     GENERATE RIGHT SIDE OF EQUATION.
  !
  DO i = 1 , mp1
    DO j = 1 , np1
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
  DO i = 1 , mp1
    DO j = 1 , np1
      x = x + f(i,j) - (r(i)*z(j))**4
    ENDDO
  ENDDO
  x = x/(np1*mp1)
  DO i = 1 , mp1
    DO j = 1 , np1
      f(i,j) = f(i,j) - x
    ENDDO
  ENDDO
  err = 0.0
  DO i = 1 , mp1
    DO j = 1 , np1
      x = ABS(f(i,j)-(r(i)*z(j))**4)
      IF ( x>err ) err = x
    ENDDO
  ENDDO
  !
  Ipass = 1
  IF ( err>ermax ) Ipass = 0
  IF ( Kprint>=3.OR.(Kprint>=2.AND.Ipass==0) ) WRITE (Lun,99002) ierror , &
    pertrb , err , INT(w(1))
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
