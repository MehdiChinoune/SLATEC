!*==PCHQK5.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK PCHQK5
      SUBROUTINE PCHQK5(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--PCHQK55
!***BEGIN PROLOGUE  PCHQK5
!***PURPOSE  Test the PCH to B-spline conversion routine PCHBS.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHQK5-S, DPCHQ5-D)
!***KEYWORDS  PCHIP CONVERSION ROUTINE QUICK CHECK
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!              PCHIP QUICK CHECK NUMBER 5
!
!     TESTS THE CONVERSION ROUTINE:  PCHBS.
! *Usage:
!
!        INTEGER  LUN, KPRINT, IPASS
!
!        CALL PCHQK5 (LUN, KPRINT, IPASS)
!
! *Arguments:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
! *Description:
!
!   This routine tests a constructed data set with four different
!   KNOTYP settings.  It computes the function and derivatives of the
!   resulting B-representation via BVALU and compares with PCH data.
!
! *Caution:
!   This routine assumes BVALU has already been successfully tested.
!
!***ROUTINES CALLED  BVALU, PCHBS, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   900411  DATE WRITTEN
!   900412  Corrected minor errors in initial implementation.
!   900430  Corrected errors in prologue.
!   900501  Corrected declarations.
!   930317  Improved output formats.  (FNF)
!***END PROLOGUE  PCHQK5
!
!*Internal Notes:
!  TOL  is the tolerance to use for quantities that should only
!       theoretically be equal.
!  TOLZ is the tolerance to use for quantities that should be exactly
!       equal.
!
!**End
!
!  Declare arguments.
!
      INTEGER Lun , Kprint , Ipass
!
!  Declare externals.
!
      REAL BVALU , R1MACH
      EXTERNAL BVALU , PCHBS , R1MACH
!
!  Declare variables.
!
      INTEGER i , ierr , ifail , inbv , j , knotyp , k , N , ndim , nknots
      PARAMETER (N=9)
      REAL bcoef(2*N) , d(N) , dcalc , derr , dermax , f(N) , fcalc , ferr , 
     &     fermax , t(2*N+4) , terr , termax , tol , tolz , tsave(2*N+4) , 
     &     work(16*N) , x(N) , ZERO
      PARAMETER (ZERO=0.0E0)
      LOGICAL fail
!
!  Define relative error function.
!
      REAL ans , err , RELERR
      RELERR(err,ans) = ABS(err)/MAX(1.0E-5,ABS(ans))
!
!  Define test data.
!
      DATA x/ - 2.2E0 , -1.2E0 , -1.0E0 , -0.5E0 , -0.01E0 , 0.5E0 , 1.0E0 , 
     &     2.0E0 , 2.2E0/
      DATA f/0.0079E0 , 0.2369E0 , 0.3679E0 , 0.7788E0 , 0.9999E0 , 0.7788E0 , 
     &     0.3679E0 , 0.1083E0 , 0.0079E0/
      DATA d/0.0000E0 , 0.3800E0 , 0.7173E0 , 0.5820E0 , 0.0177E0 , -0.5696E0 , 
     &     -0.5135E0 , -0.0778E0 , -0.0025E0/
!
!  Initialize.
!
!***FIRST EXECUTABLE STATEMENT  PCHQK5
      ifail = 0
      tol = 100*R1MACH(4)
      tolz = ZERO
!
      IF ( Kprint>=3 ) WRITE (Lun,99001)
!
!  FORMATS.
!
99001 FORMAT ('1'//10X,'TEST PCH TO B-SPLINE CONVERTER')
      IF ( Kprint>=2 ) WRITE (Lun,99002)
99002 FORMAT (//10X,'PCHQK5 RESULTS'/10X,'--------------')
!
!  Loop over a series of values of KNOTYP.
!
      IF ( Kprint>=3 ) WRITE (Lun,99003)
99003 FORMAT (/4X,'(Results should be the same for all KNOTYP values.)')
      DO knotyp = 2 , -1 , -1
!        ------------
        CALL PCHBS(N,x,f,d,1,knotyp,nknots,t,bcoef,ndim,k,ierr)
!        ------------
        IF ( Kprint>=3 ) WRITE (Lun,99004) knotyp , nknots , ndim , k , ierr
99004   FORMAT (/4X,'KNOTYP =',I2,':  NKNOTS =',I3,',  NDIM =',I3,',  K =',I2,
     &          ',  IERR =',I3)
        IF ( ierr/=0 ) THEN
          ifail = ifail + 1
          IF ( Kprint>=3 ) WRITE (Lun,99005)
99005     FORMAT (' *** Failed -- bad IERR value.')
        ELSE
!             Compare evaluated results with inputs to PCHBS.
          inbv = 1
          fermax = ZERO
          dermax = ZERO
          IF ( Kprint>=3 ) THEN
            WRITE (Lun,99006)
99006       FORMAT (/15X,'X',9X,'KNOTS',10X,'F',7X,'FERR',8X,'D',7X,'DERR')
            WRITE (Lun,99013) t(1) , t(2)
            j = 1
          ENDIF
          DO i = 1 , N
            fcalc = BVALU(t,bcoef,ndim,k,0,x(i),inbv,work)
            ferr = f(i) - fcalc
            fermax = MAX(fermax,RELERR(ferr,f(i)))
            dcalc = BVALU(t,bcoef,ndim,k,1,x(i),inbv,work)
            derr = d(i) - dcalc
            dermax = MAX(dermax,RELERR(derr,d(i)))
            IF ( Kprint>=3 ) THEN
              j = j + 2
              WRITE (Lun,99007) x(i) , t(j) , t(j+1) , f(i) , ferr , d(i) , derr
99007         FORMAT (10X,3F8.2,F10.4,1P,E10.2,0P,F10.4,1P,E10.2)
            ENDIF
          ENDDO
          IF ( Kprint>=3 ) THEN
            j = j + 2
            WRITE (Lun,99013) t(j) , t(j+1)
          ENDIF
          fail = (fermax>tol) .OR. (dermax>tol)
          IF ( fail ) ifail = ifail + 1
          IF ( (Kprint>=3).OR.(Kprint>=2).AND.fail ) WRITE (Lun,99008) fermax , 
     &         dermax , tol
99008     FORMAT (/5X,'Maximum relative errors:'/15X,'F-error =',1P,E13.5,5X,
     &            'D-error =',E13.5/5X,'Both should be less than  TOL =',E13.5)
        ENDIF
!
!          Special check for KNOTYP=-1.
        IF ( knotyp==0 ) THEN
!             Save knot vector for next test.
          DO i = 1 , nknots
            tsave(i) = t(i)
          ENDDO
        ELSEIF ( knotyp==-1 ) THEN
!             Check that knot vector is unchanged.
          termax = ZERO
          DO i = 1 , nknots
            terr = ABS(t(i)-tsave(i))
            termax = MAX(termax,terr)
          ENDDO
          IF ( termax>tolz ) THEN
            ifail = ifail + 1
            IF ( Kprint>=2 ) WRITE (Lun,99009) termax , tolz
99009       FORMAT (/' *** T-ARRAY MAXIMUM CHANGE =',1P,E13.5,
     &              ';  SHOULD NOT EXCEED TOLZ =',E13.5)
          ENDIF
        ENDIF
      ENDDO
!
!  PRINT SUMMARY AND TERMINATE.
!
      IF ( (Kprint>=2).AND.(ifail/=0) ) WRITE (Lun,99010) ifail
99010 FORMAT (/' *** TROUBLE ***',I5,' CONVERSION TESTS FAILED.')
!
      IF ( ifail==0 ) THEN
        Ipass = 1
        IF ( Kprint>=2 ) WRITE (Lun,99011)
99011   FORMAT (/' ------------  PCHIP PASSED  ALL CONVERSION TESTS',
     &          ' ------------')
      ELSE
        Ipass = 0
        IF ( Kprint>=1 ) WRITE (Lun,99012)
99012   FORMAT (/' ************  PCHIP FAILED SOME CONVERSION TESTS',
     &          ' ************')
      ENDIF
!
      RETURN
99013 FORMAT (18X,2F8.2)
!------------- LAST LINE OF PCHQK5 FOLLOWS -----------------------------
      END SUBROUTINE PCHQK5
