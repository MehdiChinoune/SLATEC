!*==OUTERR.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK OUTERR
      SUBROUTINE OUTERR(Method,Ierr,Iout,Nfail,Istdo,Iter,Err)
      IMPLICIT NONE
!*--OUTERR5
!***BEGIN PROLOGUE  OUTERR
!***SUBSIDIARY
!***PURPOSE  Output error messages for the SLAP Quick Check.
!***LIBRARY   SLATEC (SLAP)
!***TYPE      SINGLE PRECISION (OUTERR-S, DUTERR-D)
!***AUTHOR  Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-300
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   881010  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890920  Converted prologue to SLATEC 4.0 format.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   921021  Added 1P's to output formats.  (FNF)
!***END PROLOGUE  OUTERR
!     .. Scalar Arguments ..
      REAL Err
      INTEGER Ierr , Iout , Istdo , Iter , Nfail
      CHARACTER Method*6
!***FIRST EXECUTABLE STATEMENT  OUTERR
      IF ( Ierr/=0 ) Nfail = Nfail + 1
      IF ( Iout==1.AND.Ierr/=0 ) THEN
        WRITE (Istdo,99001) Method
99001   FORMAT (1X,A6,' : **** FAILURE ****')
      ENDIF
      IF ( Iout==2 ) THEN
        IF ( Ierr==0 ) THEN
          WRITE (Istdo,99002) Method
99002     FORMAT (1X,A6,' : **** PASSED  ****')
        ELSE
          WRITE (Istdo,99004) Method , Ierr , Iter , Err
        ENDIF
      ENDIF
      IF ( Iout>=3 ) THEN
        IF ( Ierr==0 ) THEN
          WRITE (Istdo,99003) Method , Ierr , Iter , Err
99003     FORMAT (' ***************** PASSED ***********************'/' **** ',
     &            A6,' Quick Test PASSED: IERR = ',I5,
     &            ' ****'/' ***************** PASSED ***********************'/
     &            ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
        ELSE
          WRITE (Istdo,99004) Method , Ierr , Iter , Err
        ENDIF
      ENDIF
      RETURN
99004 FORMAT (' **************** WARNING ***********************'/' **** ',A6,
     &        ' Quick Test FAILED: IERR = ',I5,
     &        ' ****'/' **************** WARNING ***********************'/
     &        ' Iteration Count = ',I3,' Stop Test = ',1P,E12.6)
!------------- LAST LINE OF OUTERR FOLLOWS ----------------------------
      END SUBROUTINE OUTERR
