!*==QXBLKT.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QXBLKT
      SUBROUTINE QXBLKT(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--QXBLKT5
!*** Start of declarations inserted by SPAG
      REAL am , an , bm , bn , cm , cn , deltas , deltat , ermax , err , hds , 
     &     hdt , s , t , tds , tdt , temp1 , temp2 , temp3 , w
      REAL y , z
      INTEGER i , idimy , ierror , iflg , Ipass , j , Kprint , Lun , m , mp , 
     &        n , np
!*** End of declarations inserted by SPAG
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
      DIMENSION y(75,105) , am(75) , bm(75) , cm(75) , an(105) , bn(105) , 
     &          cn(105) , w(1952) , s(75) , t(105)
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
      DO i = 1 , m
        s(i) = i*deltas
      ENDDO
      deltat = 1.0E0/(n+1)
      DO j = 1 , n
        t(j) = j*deltat
      ENDDO
!
!     COMPUTE THE COEFFICIENTS AM, BM AND CM CORRESPONDING TO THE S
!     DIRECTION.
!
      hds = deltas/2.
      tds = deltas + deltas
      DO i = 1 , m
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
      DO j = 1 , n
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
      DO j = 1 , n
        DO i = 1 , m
          y(i,j) = 3.75*s(i)*t(j)*(s(i)**4.+t(j)**4.)
        ENDDO
      ENDDO
!
!     INCLUDE NONHOMOGENEOUS BOUNDARY INTO RIGHT SIDE. NOTE THAT THE
!     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
!
      DO j = 1 , n
        y(m,j) = y(m,j) - cm(m)*t(j)**5.
      ENDDO
      DO i = 1 , m
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
          DO j = 1 , n
            DO i = 1 , m
              z = ABS(y(i,j)-(s(i)*t(j))**5.)
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
99001       FORMAT ('1',20X,'SUBROUTINE BLKTRI EXAMPLE'///10X,
     &              'THE OUTPUT FROM THE NCAR CONTROL DATA 7600 WAS'//32X,
     &              'IERROR = 0'/18X,'DISCRETIZATION ERROR = 1.6478E-05'/12X,
     &              'REQUIRED LENGTH OF W ARRAY = 823'//10X,
     &              'THE OUTPUT FROM YOUR COMPUTER IS'//32X,'IERROR =',I2/18X,
     &              'DISCRETIZATION ERROR =',1PE12.5/12X,
     &              'REQUIRED LENGTH OF W ARRAY =',I4)
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
