!*==QXSSP.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXSSP
      SUBROUTINE QXSSP(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--QXSSP5
!*** Start of declarations inserted by SPAG
      REAL bdpf , bdps , bdtf , bdts , dphi , dtheta , dum , elmbda , ermax , 
     &     err , f , pertrb , pf , pi , PIMACH , ps , sinp , sint , tf , ts
      REAL w , z
      INTEGER i , idimf , ierror , Ipass , j , Kprint , Lun , m , mbdcnd , mp1 , 
     &        n , nbdcnd , np1
!*** End of declarations inserted by SPAG
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
      DIMENSION f(19,73) , bdtf(73) , sint(19) , sinp(73) , w(1200)
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
      DO i = 1 , mp1
        sint(i) = SIN((i-1)*dtheta)
      ENDDO
      dphi = (pi+pi)/n
      np1 = n + 1
      DO j = 1 , np1
        sinp(j) = SIN((j-1)*dphi)
      ENDDO
!
!     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
!
      DO j = 1 , np1
        DO i = 1 , mp1
          f(i,j) = 2. - 6.*(sint(i)*sinp(j))**2
        ENDDO
      ENDDO
!
!     STORE DERIVATIVE DATA AT THE EQUATOR
!
      DO j = 1 , np1
        bdtf(j) = 0.
      ENDDO
!
      CALL HWSSSP(ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,bdpf,elmbda,f,
     &            idimf,pertrb,ierror,w)
!
!     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
!     SOLUTION MUST BE NORMALIZED.
!
      err = 0.0
      DO j = 1 , np1
        DO i = 1 , mp1
          z = ABS(f(i,j)-(sint(i)*sinp(j))**2-f(1,1))
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
99001   FORMAT ('1',20X,'SUBROUTINE HWSSSP EXAMPLE'///10X,
     &          'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,
     &          'IERROR = 0'/18X,'DISCRETIZATION ERROR = 3.38107E-03'/12X,
     &          'REQUIRED LENGTH OF W ARRAY = 600'//10X,
     &          'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,
     &          'DISCRETIZATION ERROR =',1PE12.5/12X,
     &          'REQUIRED LENGTH OF W ARRAY =',I4)
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
