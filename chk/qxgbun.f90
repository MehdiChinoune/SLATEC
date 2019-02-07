!*==QXGBUN.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXGBUN
SUBROUTINE QXGBUN(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXGBUN5
  !*** Start of declarations inserted by SPAG
  REAL a , b , c , deltax , deltay , dum , dysq , ermax , err , f , pi , &
    PIMACH , s , t , w , x , y , z
  INTEGER i , idimy , ierror , Ipass , j , Kprint , Lun , m , mm1 , mperod , &
    n , nperod
  !*** End of declarations inserted by SPAG
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
  DIMENSION f(25,130) , a(20) , b(20) , c(20) , w(1200) , x(20) , y(120)
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
  DO i = 1 , m
    x(i) = (i-1)*deltax
  ENDDO
  DO j = 1 , n
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
  DO i = 2 , m
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
  DO j = 1 , n
    f(1,j) = dysq*(11.+8./deltax)*SIN(y(j))
  ENDDO
  !
  !     GENERATE RIGHT SIDE.
  !
  mm1 = m - 1
  DO i = 2 , mm1
    DO j = 1 , n
      f(i,j) = dysq*3.*(1.+x(i))**4*SIN(y(j))
    ENDDO
  ENDDO
  !
  !     GENERATE RIGHT SIDE FOR I = M SHOWING INTRODUCTION OF
  !     BOUNDARY DATA.
  !
  DO j = 1 , n
    f(m,j) = dysq*(3.*(1.+x(m))**4-16.*((1.+x(m))/deltax)**2+16.*(1.+x(m))&
      /deltax)*SIN(y(j))
  ENDDO
  CALL GENBUN(nperod,n,mperod,m,a,b,c,idimy,f,ierror,w)
  !
  !     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
  !                   U(X,Y) = (1+X)**4*SIN(Y)
  !
  err = 0.
  DO i = 1 , m
    DO j = 1 , n
      z = ABS(f(i,j)-(1.+x(i))**4*SIN(y(j)))
      IF ( z>err ) err = z
    ENDDO
  ENDDO
  !
  Ipass = 1
  IF ( err>ermax ) Ipass = 0
  IF ( Kprint==0 ) RETURN
  IF ( Kprint>=2.OR.Ipass==0 ) THEN
    WRITE (Lun,99001) ierror , err , INT(w(1))
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
