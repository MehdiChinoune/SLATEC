!DECK QCPSI
SUBROUTINE QCPSI(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER Ipass, Kprint
  !***BEGIN PROLOGUE  QCPSI
  !***PURPOSE  Quick check for PSIFN.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     ABSTRACT
  !     QCPSI IS A QUICK CHECK ROUTINE WHICH EXERCISES THE MAJOR
  !     LOOPS IN SUBROUTINE PSIFN(X,N,KODE,M,ANS,NZ,IERR) FOR DERIVATIVES
  !     OF THE PSI FUNCTION.  FOR N=0, THE PSI FUNCTIONS ARE CALCULATED
  !     EXPLICITLY AND CHECKED AGAINST EVALUATIONS FROM PSIFN. FOR
  !     N.GT.0, CONSISTENCY CHECKS ARE MADE BY COMPARING A SEQUENCE
  !     AGAINST SINGLE EVALUATIONS OF PSIFN, ONE AT A TIME.
  !     IF THE RELATIVE ERROR IS LESS THAN 1000 TIMES UNIT ROUNDOFF,
  !     THEN THE TEST IS PASSED--IF NOT,
  !     THEN X, THE VALUES TO BE COMPARED, THE RELATIVE ERROR AND
  !     PARAMETERS KODE AND N ARE WRITTEN ON LOGICAL UNIT 6 WHERE N IS
  !     THE ORDER OF THE DERIVATIVE AND KODE IS A SELECTION PARAMETER
  !     DEFINED IN THE PROLOGUE TO PSIFN.
  !
  !     FUNCTIONS I1MACH AND R1MACH MUST BE INITIALIZED ACCORDING TO THE
  !     PROLOGUE IN EACH FUNCTION FOR THE MACHINE ENVIRONMENT BEFORE
  !     QCPSI OR PSIFN CAN BE EXECUTED.
  !
  !***ROUTINES CALLED  PSIFN, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  QCPSI
  INTEGER i, ierr, iflg, ix, kode, Lun, m, n, nm, nn, nz
  REAL er, euler, psi1, psi2, r1m4, s, tol, x
  REAL R1MACH
  DIMENSION psi1(3), psi2(20)
  DATA euler/0.5772156649015328606E0/
  !***FIRST EXECUTABLE STATEMENT  QCPSI
  r1m4 = R1MACH(4)
  tol = 1000.0E0*MAX(r1m4,1.0E-18)
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  99001 FORMAT ('1'//' QUICK CHECK DIAGNOSTICS FOR PSIFN'//)
  !-----------------------------------------------------------------------
  !     CHECK PSI(I) AND PSI(I-0.5), I=1,2,...
  !-----------------------------------------------------------------------
  iflg = 0
  n = 0
  DO kode = 1, 2
    DO m = 1, 2
      s = -euler + (m-1)*(-2.0E0*LOG(2.0E0))
      x = 1.0E0 - (m-1)*0.5E0
      DO i = 1, 20
        CALL PSIFN(x,n,kode,1,psi2,nz,ierr)
        psi1(1) = -s + (kode-1)*LOG(x)
        er = ABS((psi1(1)-psi2(1))/psi1(1))
        IF ( er>tol ) THEN
          IF ( iflg==0 ) THEN
            IF ( Kprint>=2 ) WRITE (Lun,99004)
          ENDIF
          iflg = iflg + 1
          IF ( Kprint>=2 ) WRITE (Lun,99005) x, psi1(1), psi2(i), er, &
            kode, n
          IF ( iflg>200 ) GOTO 100
        ENDIF
        s = s + 1.0E0/x
        x = x + 1.0E0
      ENDDO
    ENDDO
  ENDDO
  !-----------------------------------------------------------------------
  !     CHECK SMALL X.LT.UNIT ROUNDOFF
  !-----------------------------------------------------------------------
  kode = 1
  x = tol/10000.0E0
  n = 1
  CALL PSIFN(x,n,kode,1,psi2,nz,ierr)
  psi1(1) = x**(-n-1)
  er = ABS((psi1(1)-psi2(1))/psi1(1))
  IF ( er>tol ) THEN
    IF ( iflg==0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99004)
    ENDIF
    iflg = iflg + 1
    IF ( Kprint>=2 ) WRITE (Lun,99005) x, psi1(1), psi2(1), er, kode, n
  ENDIF
  !-----------------------------------------------------------------------
  !     CONSISTENCY TESTS FOR N.GE.0
  !-----------------------------------------------------------------------
  DO kode = 1, 2
    DO m = 1, 5
      DO n = 1, 16, 5
        nn = n - 1
        x = 0.1E0
        DO ix = 1, 25, 2
          x = x + 1.0E0
          CALL PSIFN(x,nn,kode,m,psi2,nz,ierr)
          DO i = 1, m
            nm = nn + i - 1
            CALL PSIFN(x,nm,kode,1,psi1,nz,ierr)
            er = ABS((psi2(i)-psi1(1))/psi1(1))
            IF ( er>=tol ) THEN
              IF ( iflg==0 ) THEN
                IF ( Kprint>=2 ) WRITE (Lun,99004)
              ENDIF
              iflg = iflg + 1
              IF ( Kprint>=2 ) WRITE (Lun,99005) x, psi1(1), psi2(i), &
                er, kode, nm
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  IF ( iflg==0.AND.Kprint>=3 ) THEN
    WRITE (Lun,99002)
    99002   FORMAT (//' QUICK CHECKS OK'//)
  ENDIF
  Ipass = 0
  IF ( iflg==0 ) Ipass = 1
  RETURN
  100 CONTINUE
  IF ( Kprint>=2 ) WRITE (Lun,99003)
  99003 FORMAT (//' PROCESSING OF MAIN LOOPS TERMINATED BECAUSE THE NUM',&
    'BER OF DIAGNOSTIC PRINTS EXCEEDS 200'//)
  Ipass = 0
  IF ( iflg==0 ) Ipass = 1
  99004 FORMAT (8X,'X',13X,'PSI1',11X,'PSI2',9X,'REL ERR',5X,'KODE',3X,'N')
  99005 FORMAT (4E15.6,2I5)
END SUBROUTINE QCPSI
