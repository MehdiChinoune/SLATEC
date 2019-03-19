!** XERBLA
SUBROUTINE XERBLA(Srname,Info)
  IMPLICIT NONE
  !>
  !***
  !  Error handler for the Level 2 and Level 3 BLAS Routines.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  R3
  !***
  ! **Type:**      ALL (XERBLA-A)
  !***
  ! **Keywords:**  ERROR MESSAGE
  !***
  ! **Author:**  Dongarra, J. J., (ANL)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   860720  DATE WRITTEN
  !   910610  Routine rewritten to serve as an interface between the
  !           Level 2 and Level 3 BLAS routines and the SLATEC error
  !           handler XERMSG.  (BKS)
  
  !
  !     ..    Scalar Arguments ..
  INTEGER Info
  CHARACTER(6) :: Srname
  CHARACTER(2) :: xern1
  !
  !* FIRST EXECUTABLE STATEMENT  XERBLA
  !
  WRITE (xern1,'(I2)') Info
  CALL XERMSG('SLATEC',Srname,'On entry to '//Srname//' parameter number '//&
    xern1//' had an illegal value',Info,1)
  !
  !
  !     End of XERBLA.
  !
END SUBROUTINE XERBLA
