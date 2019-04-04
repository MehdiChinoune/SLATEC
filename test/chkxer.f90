!** CHKXER
SUBROUTINE CHKXER(Srnamt,Infot,Nout,Fatal,Kprint)
  IMPLICIT NONE
  !>
  !***
  !  Test whether an error has been detected.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Author:**  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !  Tests whether an error has been detected.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  NUMXER

  !* REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)

  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Infot, Kprint, Nout
  CHARACTER(6) :: Srnamt
  !     .. Local Scalars ..
  INTEGER nerr
  !     .. External Functions ..
  INTEGER, EXTERNAL :: NUMXER
  !* FIRST EXECUTABLE STATEMENT  CHKXER
  IF ( NUMXER(nerr)/=Infot ) THEN
    Fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Nout,FMT=99001) Infot, Srnamt
  END IF
  RETURN
  !
  99001 FORMAT (' ***** ILLEGAL VALUE OF PARAMETER NUMBER ',I2,' NOT D',&
    'ETECTED BY ',A6,' *****')
  !
  !     End of CHKXER.
  !
END SUBROUTINE CHKXER
