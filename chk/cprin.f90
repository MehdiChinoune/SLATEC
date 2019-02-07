!*==CPRIN.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CPRIN
SUBROUTINE CPRIN(Lun,Num1,Kprint,Ip,Exact,Result,Abserr,Neval,Ierv,Lierv)
  IMPLICIT NONE
  !*--CPRIN5
  !***BEGIN PROLOGUE  CPRIN
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CQAG, CQAG, CQAGI, CQAGP, CQAGS, CQAWC,
  !            CQAWF, CQAWO, CQAWS, and CQNG.
  !***LIBRARY   SLATEC
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !   This program is called by the (single precision) Quadpack quick
  !   check routines for printing out their messages.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810401  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   910627  Code completely rewritten.  (WRB)
  !***END PROLOGUE  CPRIN
  !     .. Scalar Arguments ..
  REAL Abserr, Exact, Result
  INTEGER Ip, Kprint, Lierv, Lun, Neval, Num1
  !     .. Array Arguments ..
  INTEGER Ierv(*)
  !     .. Local Scalars ..
  REAL error
  INTEGER ier, k
  !     .. Intrinsic Functions ..
  INTRINSIC ABS
  !***FIRST EXECUTABLE STATEMENT  CPRIN
  ier = Ierv(1)
  error = ABS(Exact-Result)
  !
  IF ( Kprint>=2 ) THEN
    IF ( Ip/=1 ) THEN
      !
      !         Write failure messages.
      !
      WRITE (UNIT=Lun,FMT=99002) Num1
      IF ( Num1==0 ) WRITE (UNIT=Lun,FMT=99003)
      IF ( Num1>0 ) WRITE (UNIT=Lun,FMT=99004) Num1
      IF ( Lierv>1 ) WRITE (UNIT=Lun,FMT=99005) (Ierv(k),k=2,Lierv)
      IF ( Num1==6 ) WRITE (UNIT=Lun,FMT=99006)
      WRITE (UNIT=Lun,FMT=99007)
      WRITE (UNIT=Lun,FMT=99008)
      IF ( Num1/=5 ) THEN
        WRITE (UNIT=Lun,FMT=99009) Exact, Result, error, Abserr, ier, &
          Neval
      ELSE
        WRITE (Lun,FMT=99010) Result, Abserr, ier, Neval
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      !
      !           Write PASS message.
      !
      WRITE (UNIT=Lun,FMT=99001) Num1
    ENDIF
  ENDIF
  !
  RETURN
  !
  99001 FORMAT (' TEST ON IER = ',I2,' PASSED')
  99002 FORMAT (' TEST ON IER = ',I1,' FAILED.')
  99003 FORMAT (' WE MUST HAVE IER = 0, ERROR.LE.ABSERR AND ABSERR.LE',&
    '.MAX(EPSABS,EPSREL*ABS(EXACT))')
  99004 FORMAT (' WE MUST HAVE IER = ',I1)
  99005 FORMAT (' OR IER =     ',8(I1,2X))
  99006 FORMAT (' RESULT, ABSERR, NEVAL AND EVENTUALLY LAST SHOULD BE',' ZERO')
  99007 FORMAT (' WE HAVE   ')
  99008 FORMAT (7X,'EXACT',11X,'RESULT',6X,'ERROR',4X,'ABSERR',4X,'IER     NEVAL',&
    /,' ',42X,'(EST.ERR.)(FLAG)(NO F-EVAL)')
  99009 FORMAT (' ',2(E15.7,1X),2(E9.2,1X),I4,4X,I6)
  99010 FORMAT (5X,'INFINITY',4X,E15.7,11X,E9.2,I5,4X,I6)
END SUBROUTINE CPRIN
