!*==CHKXER.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CHKXER
SUBROUTINE CHKXER(Srnamt,Infot,Nout,Fatal,Kprint)
  IMPLICIT NONE
  !*--CHKXER5
  !***BEGIN PROLOGUE  CHKXER
  !***SUBSIDIARY
  !***PURPOSE  Test whether an error has been detected.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Tests whether an error has been detected.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CHKXER
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  INTEGER Infot , Kprint , Nout
  CHARACTER(6) :: Srnamt
  !     .. Local Scalars ..
  INTEGER nerr
  !     .. External Functions ..
  INTEGER NUMXER
  EXTERNAL NUMXER
  !***FIRST EXECUTABLE STATEMENT  CHKXER
  IF ( NUMXER(nerr)/=Infot ) THEN
    Fatal = .TRUE.
    IF ( Kprint>=3 ) WRITE (Nout,FMT=99001) Infot , Srnamt
  ENDIF
  RETURN
  !
  99001 FORMAT (' ***** ILLEGAL VALUE OF PARAMETER NUMBER ',I2,' NOT D',&
    'ETECTED BY ',A6,' *****')
  !
  !     End of CHKXER.
  !
END SUBROUTINE CHKXER
