!*==PCHQK1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK PCHQK1
      SUBROUTINE PCHQK1(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--PCHQK15
!***BEGIN PROLOGUE  PCHQK1
!***PURPOSE  Test the PCHIP evaluators CHFDV, CHFEV, PCHFD and PCHFE.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHQK1-S, DPCHQ1-D)
!***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!              PCHIP QUICK CHECK NUMBER 1
!
!     TESTS THE EVALUATORS:  CHFDV, CHFEV, PCHFD, PCHFE.
! *Usage:
!
!        INTEGER  LUN, KPRINT, IPASS
!
!        CALL PCHQK1 (LUN, KPRINT, IPASS)
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
!   This routine carries out three tests of the PCH evaluators:
!     EVCHCK tests the single-cubic evaluators.
!     EVPCCK tests the full PCH evaluators.
!     EVERCK exercises the error returns in all evaluators.
!
!***ROUTINES CALLED  EVCHCK, EVERCK, EVPCCK
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
!   890618  REVISION DATE from Version 3.2
!   890706  Cosmetic changes to prologue.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900309  Added EVERCK to list of routines called.  (FNF)
!   900314  Improved some output formats.
!   900315  Revised prologue and improved some output formats.  (FNF)
!   900316  Additional minor cosmetic changes.  (FNF)
!   900321  Removed IFAIL from call sequence for SLATEC standards and
!           made miscellaneous cosmetic changes.  (FNF)
!   930317  Improved output formats.  (FNF)
!***END PROLOGUE  PCHQK1
!
!  Declare arguments.
!
      INTEGER Lun , Kprint , Ipass
!
!  DECLARE LOCAL VARIABLES.
!
      INTEGER i1 , i2 , i3 , i4 , i5 , i6 , i7 , i8 , i9 , ifail , npts
      REAL work(4000)
      LOGICAL fail
!
!***FIRST EXECUTABLE STATEMENT  PCHQK1
      IF ( Kprint>=2 ) WRITE (Lun,99001) Kprint
!
!  FORMATS.
!
99001 FORMAT ('1'/' ------------  PCHIP QUICK CHECK OUTPUT',' ------------'//
     &        20X,'( KPRINT =',I2,' )')
!
!  TEST CHFDV AND CHFEV.
!
      ifail = 0
      npts = 1000
      i1 = 1 + npts
      i2 = i1 + npts
      i3 = i2 + npts
      CALL EVCHCK(Lun,Kprint,npts,work(1),work(i1),work(i2),work(i3),fail)
      IF ( fail ) ifail = ifail + 1
!
!  TEST PCHFD AND PCHFE.
!
      i1 = 1 + 10
      i2 = i1 + 10
      i3 = i2 + 100
      i4 = i3 + 100
      i5 = i4 + 100
      i6 = i5 + 51
      i7 = i6 + 51
      i8 = i7 + 51
      i9 = i8 + 51
      CALL EVPCCK(Lun,Kprint,work(1),work(i1),work(i2),work(i3),work(i4),
     &            work(i5),work(i6),work(i7),work(i8),work(i9),fail)
      IF ( fail ) ifail = ifail + 2
!
!  TEST ERROR RETURNS.
!
      CALL EVERCK(Lun,Kprint,fail)
      IF ( fail ) ifail = ifail + 4
!
!  PRINT SUMMARY AND TERMINATE.
!     At this point, IFAIL has the following value:
!        IFAIL = 0  IF ALL TESTS PASSED.
!        IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
!           IFAIL=1  IF SINGLE CUBIC TEST FAILED. (SEE EVCHCK OUTPUT.)
!           IFAIL=2  IF PCHFD/PCHFE  TEST FAILED. (SEE EVPCCK OUTPUT.)
!           IFAIL=4  IF ERROR RETURN TEST FAILED. (SEE EVERCK OUTPUT.)
!
      IF ( (Kprint>=2).AND.(ifail/=0) ) WRITE (Lun,99002) ifail
99002 FORMAT (/' *** TROUBLE ***',I5,' EVALUATION TESTS FAILED.')
!
      IF ( ifail==0 ) THEN
        Ipass = 1
        IF ( Kprint>=2 ) WRITE (Lun,99003)
99003   FORMAT (/' ------------  PCHIP PASSED  ALL EVALUATION TESTS',
     &          ' ------------')
      ELSE
        Ipass = 0
        IF ( Kprint>=1 ) WRITE (Lun,99004)
99004   FORMAT (/' ************  PCHIP FAILED SOME EVALUATION TESTS',
     &          ' ************')
      ENDIF
!
      RETURN
!------------- LAST LINE OF PCHQK1 FOLLOWS -----------------------------
      END SUBROUTINE PCHQK1
