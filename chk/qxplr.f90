!*==QXPLR.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXPLR
SUBROUTINE QXPLR(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QXPLR5
  !*** Start of declarations inserted by SPAG
  REAL a , b , bda , bdb , bdc , bdd , c , d , dum , elmbda , ermax , err , &
    f , pertrb , pi , PIMACH , r , theta , w , z
  INTEGER i , idimf , ierror , Ipass , j , Kprint , Lun , m , mbdcnd , mp1 , &
    n , nbdcnd , np1
  !*** End of declarations inserted by SPAG
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
  DIMENSION f(100,50) , bdc(51) , bdd(51) , w(1200) , r(51) , theta(49)
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
  DO i = 1 , mp1
    r(i) = (i-1)/50.0E0
  ENDDO
  DO j = 1 , np1
    theta(j) = (j-1)*pi/96.0E0
  ENDDO
  !
  !     GENERATE BOUNDARY DATA.
  !
  DO i = 1 , mp1
    bdc(i) = 0.
    bdd(i) = 0.
  ENDDO
  !
  !     BDA AND BDB ARE DUMMY VARIABLES.
  !
  DO j = 1 , np1
    f(mp1,j) = 1. - COS(4.*theta(j))
  ENDDO
  !
  !     GENERATE RIGHT SIDE OF EQUATION.
  !
  DO i = 1 , m
    DO j = 1 , np1
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
  DO i = 1 , mp1
    DO j = 1 , np1
      z = ABS(f(i,j)-r(i)**4*(1.-COS(4.*theta(j))))
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
