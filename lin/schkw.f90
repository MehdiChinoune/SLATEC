!*==SCHKW.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SCHKW
SUBROUTINE SCHKW(Name,Lociw,Leniw,Locw,Lenw,Ierr,Iter,Err)
  IMPLICIT NONE
  !*--SCHKW5
  !***BEGIN PROLOGUE  SCHKW
  !***SUBSIDIARY
  !***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
  !            This routine checks the work array lengths and interfaces
  !            to the SLATEC error handler if a problem is found.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  R2
  !***TYPE      SINGLE PRECISION (SCHKW-S, DCHKW-D)
  !***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING
  !***AUTHOR  Seager, Mark K., (LLNL)
  !             Lawrence Livermore National Laboratory
  !             PO BOX 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !             seager@llnl.gov
  !***DESCRIPTION
  !
  ! *Usage:
  !     CHARACTER*(*) NAME
  !     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
  !     REAL    ERR
  !
  !     CALL SCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
  !
  ! *Arguments:
  ! NAME   :IN       Character*(*).
  !         Name of the calling routine.  This is used in the output
  !         message, if an error is detected.
  ! LOCIW  :IN       Integer.
  !         Location of the first free element in the integer workspace
  !         array.
  ! LENIW  :IN       Integer.
  !         Length of the integer workspace array.
  ! LOCW   :IN       Integer.
  !         Location of the first free element in the real workspace
  !         array.
  ! LENRW  :IN       Integer.
  !         Length of the real workspace array.
  ! IERR   :OUT      Integer.
  !         Return error flag.
  !               IERR = 0 => All went well.
  !               IERR = 1 => Insufficient storage allocated for
  !                           WORK or IWORK.
  ! ITER   :OUT      Integer.
  !         Set to zero on return.
  ! ERR    :OUT      Real.
  !         Set to the smallest positive magnitude if all went well.
  !         Set to a very large number if an error is detected.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   880225  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI
  !           X3.9-1978.  (FNF)
  !   910506  Made subsidiary.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)
  !***END PROLOGUE  SCHKW
  !     .. Scalar Arguments ..
  REAL Err
  INTEGER Ierr, Iter, Leniw, Lenw, Lociw, Locw
  CHARACTER Name*(*)
  !     .. Local Scalars ..
  CHARACTER xern1*8, xern2*8, xernam*8
  !     .. External Functions ..
  REAL R1MACH
  EXTERNAL R1MACH
  !     .. External Subroutines ..
  EXTERNAL XERMSG
  !***FIRST EXECUTABLE STATEMENT  SCHKW
  !
  !         Check the Integer workspace situation.
  !
  Ierr = 0
  Iter = 0
  Err = R1MACH(1)
  IF ( Lociw>Leniw ) THEN
    Ierr = 1
    Err = R1MACH(2)
    xernam = Name
    WRITE (xern1,'(I8)') Lociw
    WRITE (xern2,'(I8)') Leniw
    CALL XERMSG('SLATEC','SCHKW',&
      'In '//xernam//', INTEGER work array too short.  '//&
      'IWORK needs '//xern1//'; have allocated '//xern2,1,1)
  ENDIF
  !
  !         Check the Real workspace situation.
  IF ( Locw>Lenw ) THEN
    Ierr = 1
    Err = R1MACH(2)
    xernam = Name
    WRITE (xern1,'(I8)') Locw
    WRITE (xern2,'(I8)') Lenw
    CALL XERMSG('SLATEC','SCHKW',&
      'In '//xernam//', REAL work array too short.  '//&
      'RWORK needs '//xern1//'; have allocated '//xern2,1,1)
  ENDIF
  !------------- LAST LINE OF SCHKW FOLLOWS ----------------------------
END SUBROUTINE SCHKW
