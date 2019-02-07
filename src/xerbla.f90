!*==XERBLA.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK XERBLA
SUBROUTINE XERBLA(Srname,Info)
  IMPLICIT NONE
  !*--XERBLA5
  !***BEGIN PROLOGUE  XERBLA
  !***SUBSIDIARY
  !***PURPOSE  Error handler for the Level 2 and Level 3 BLAS Routines.
  !***LIBRARY   SLATEC
  !***CATEGORY  R3
  !***TYPE      ALL (XERBLA-A)
  !***KEYWORDS  ERROR MESSAGE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !***DESCRIPTION
  !
  !  Purpose
  !  =======
  !
  !  It is called by Level 2 and 3 BLAS routines if an input parameter
  !  is invalid.
  !
  !  Parameters
  !  ==========
  !
  !  SRNAME - CHARACTER*6.
  !           On entry, SRNAME specifies the name of the routine which
  !           called XERBLA.
  !
  !  INFO   - INTEGER.
  !           On entry, INFO specifies the position of the invalid
  !           parameter in the parameter-list of the calling routine.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   860720  DATE WRITTEN
  !   910610  Routine rewritten to serve as an interface between the
  !           Level 2 and Level 3 BLAS routines and the SLATEC error
  !           handler XERMSG.  (BKS)
  !***END PROLOGUE  XERBLA
  !
  !     ..    Scalar Arguments ..
  INTEGER Info
  CHARACTER(6) :: Srname
  CHARACTER(2) :: xern1
  !
  !***FIRST EXECUTABLE STATEMENT  XERBLA
  !
  WRITE (xern1,'(I2)') Info
  CALL XERMSG('SLATEC',Srname,'On entry to '//Srname//' parameter number '//&
    xern1//' had an illegal value',Info,1)
  !
  !
  !     End of XERBLA.
  !
END SUBROUTINE XERBLA
