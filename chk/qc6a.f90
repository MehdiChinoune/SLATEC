!*==QC6A.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QC6A
SUBROUTINE QC6A(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QC6A5
  !***BEGIN PROLOGUE  QC6A
  !***PURPOSE  Test subroutine AAAAAA.
  !***LIBRARY   SLATEC
  !***TYPE      ALL (QC6A-A)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***DESCRIPTION
  !
  ! *Usage:
  !
  !        INTEGER  LUN, KPRINT, IPASS
  !
  !        CALL  QC6A (LUN, KPRINT, IPASS)
  !
  ! *Arguments:
  !
  !     LUN   :IN  is the unit number to which output is to be written.
  !
  !     KPRINT:IN  controls the amount of output, as specified in the
  !                SLATEC Guidelines.
  !
  !     IPASS:OUT  indicates whether the test passed or failed.
  !                A value of one is good, indicating no failures.
  !
  ! *Description:
  !
  !   This routine tests the SLATEC routine AAAAAA to see if the version
  !   number in the SLATEC library source is the same as the quick check
  !   version number.
  !
  !***ROUTINES CALLED  AAAAAA
  !***REVISION HISTORY  (YYMMDD)
  !   890713  DATE WRITTEN
  !   921215  Updated for Version 4.0.  (WRB)
  !   930701  Updated for Version 4.1.  (WRB)
  !***END PROLOGUE  QC6A
  !
  !*Internal Notes:
  !
  !     Data set-up is done via a PARAMETER statement.
  !
  !**End
  !
  !  Declare arguments.
  !
  INTEGER Lun , Kprint , Ipass
  !
  !  DECLARE VARIABLES.
  !
  CHARACTER(16) :: ver , VERSN
  PARAMETER (VERSN=' 4.1')
  !
  !***FIRST EXECUTABLE STATEMENT  QC6A
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  !
  !     Formats.
  !
  99001 FORMAT ('1'//' CODE TO TEST SLATEC ROUTINE AAAAAA')
  CALL AAAAAA(ver)
  IF ( ver==VERSN ) THEN
    Ipass = 1
    IF ( Kprint>=3 ) THEN
      WRITE (Lun,99006)
      WRITE (Lun,99002) ver
      99002     FORMAT (' *** Passed -- version number = ',A16)
    ENDIF
  ELSE
    Ipass = 0
    IF ( Kprint>=3 ) WRITE (Lun,99006)
    IF ( Kprint>=2 ) WRITE (Lun,99003) ver , VERSN
    99003   FORMAT (' *** Failed -- version number from AAAAAA = ',A16,&
      ' but expected version number = ',A16)
  ENDIF
  !
  !     Terminate.
  !
  IF ( Kprint>=2.AND.Ipass==1 ) WRITE (Lun,99004)
  99004 FORMAT (/' ************QC6A   PASSED  ALL TESTS ****************')
  IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99005)
  99005 FORMAT (/' ************QC6A   FAILED SOME TESTS ****************')
  RETURN
  99006 FORMAT (/' QC6A RESULTS')
  !------------- LAST LINE OF QC6A FOLLOWS -----------------------------
END SUBROUTINE QC6A
