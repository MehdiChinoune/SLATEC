!*==QXCSP.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXCSP
      SUBROUTINE QXCSP(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--QXCSP5
!*** Start of declarations inserted by SPAG
      REAL bdrf , bdrs , bdtf , bdts , ci4 , dphi , dr , dtheta , dum , elmbda , 
     &     ermax , err , f , pertrb , pi , PIMACH , r , rf , rs , si
      REAL tf , theta , ts , w , z
      INTEGER i , idimf , ierror , intl , Ipass , j , Kprint , Lun , m , 
     &        mbdcnd , mp1 , n , nbdcnd , np1
!*** End of declarations inserted by SPAG
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
      DIMENSION f(48,33) , bdtf(33) , w(1200) , r(33) , theta(48)
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
      DO i = 1 , mp1
        theta(i) = (i-1)*dtheta
      ENDDO
      np1 = n + 1
      dr = 1.0E0/n
      DO j = 1 , np1
        r(j) = (j-1)*dr
      ENDDO
!
!     GENERATE NORMAL DERIVATIVE DATA AT EQUATOR
!
      DO j = 1 , np1
        bdtf(j) = 0.
      ENDDO
!
!     COMPUTE BOUNDARY DATA ON THE SURFACE OF THE SPHERE
!
      DO i = 1 , mp1
        f(i,n+1) = COS(theta(i))**4
      ENDDO
!
!     COMPUTE RIGHT SIDE OF EQUATION
!
      DO i = 1 , mp1
        ci4 = 12.0E0*COS(theta(i))**2
        DO j = 1 , n
          f(i,j) = ci4*r(j)**2
        ENDDO
      ENDDO
!
      CALL HWSCSP(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,bdrf,elmbda,
     &            f,idimf,pertrb,ierror,w)
!
!     COMPUTE DISCRETIZATION ERROR
!
      err = 0.
      DO i = 1 , mp1
        ci4 = COS(theta(i))**4
        DO j = 1 , n
          z = ABS(f(i,j)-ci4*r(j)**4)
          IF ( z>err ) err = z
        ENDDO
      ENDDO
!
      Ipass = 1
      IF ( err>ermax ) Ipass = 0
      IF ( Kprint/=0 ) THEN
        IF ( Kprint>=2.OR.Ipass==0 ) THEN
          WRITE (Lun,99001) ierror , err , INT(w(1))
!
99001     FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 1'///10X,
     &            'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,
     &            'IERROR = 0'/18X,'DISCRETIZATION ERROR = 7.99842E-04'/12X,
     &            'REQUIRED LENGTH OF W ARRAY = 775'//10X,
     &            'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,
     &            'DISCRETIZATION ERROR =',1PE12.5/12X,
     &            'REQUIRED LENGTH OF W ARRAY =',I4)
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
      DO i = 1 , mp1
        f(i,n+1) = SIN(theta(i))
      ENDDO
!
!     COMPUTE RIGHT SIDE OF THE EQUATION
!
      DO j = 1 , n
        DO i = 1 , mp1
          f(i,j) = 0.
        ENDDO
      ENDDO
!
      CALL HWSCSP(intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,bdrf,elmbda,
     &            f,idimf,pertrb,ierror,w)
!
!     COMPUTE DISCRETIZATION ERROR   (FOURIER COEFFICIENTS)
!
      err = 0.
      DO i = 1 , mp1
        si = SIN(theta(i))
        DO j = 1 , np1
          z = ABS(f(i,j)-r(j)*si)
          IF ( z>err ) err = z
        ENDDO
      ENDDO
!
      IF ( err>ermax ) Ipass = 0
      IF ( Kprint==0 ) RETURN
      IF ( Kprint>=2.OR.Ipass==0 ) THEN
        WRITE (Lun,99002) ierror , err , INT(w(1))
99002   FORMAT ('1',20X,'SUBROUTINE HWSCSP EXAMPLE 2'///10X,
     &          'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,
     &          'IERROR = 0'/18X,'DISCRETIZATION ERROR = 5.86824E-05'/12X,
     &          'REQUIRED LENGTH OF W ARRAY = 775'//10X,
     &          'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,
     &          'DISCRETIZATION ERROR =',1PE12.5/12X,
     &          'REQUIRED LENGTH OF W ARRAY =',I4)
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
