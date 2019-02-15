!DECK DQCKIN
SUBROUTINE DQCKIN(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER Ipass, Kprint, nz
  !***BEGIN PROLOGUE  DQCKIN
  !***PURPOSE  Quick check for DBSKIN.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     ABSTRACT     * A DOUBLE PRECISION ROUTINE *
  !     DQCKIN IS A QUICK CHECK ROUTINE WHICH EXERCISES THE MAJOR
  !     LOOPS IN SUBROUTINE DBSKIN (X,N,KODE,M,Y,NZ,IERR) FOR BICKLEY
  !     FUNCTIONS KI(J,X).  MORE PRECISELY, DQCKIN DOES CONSISTENCY CHECKS
  !     ON THE OUTPUT FROM DBSKIN BY COMPARING SINGLE EVALUATIONS (M=1)
  !     AGAINST SELECTED MEMBERS OF SEQUENCES WHICH ARE GENERATED BY
  !     RECURSION.  IF THE RELATIVE ERROR IS LESS THAN 1000 TIMES UNIT
  !     ROUND OFF, THEN THE TEST IS PASSED - IF NOT, THEN X, THE VALUES
  !     TO BE COMPARED, THE RELATIVE ERROR AND PARAMETERS KODE, N, M AND K
  !     ARE WRITTEN ON LOGICAL UNIT 6 WHERE K IS THE MEMBER OF THE
  !     SEQUENCE OF LENGTH M WHICH FAILED THE TEST.  THAT IS, THE INDEX
  !     OF THE FUNCTION WHICH FAILED THE TEST IS J=N+K-1.  UNDERFLOW
  !     TESTS ARE MADE AND ERROR CONDITIONS ARE TRIGGERED.
  !
  !     FUNCTIONS I1MACH AND D1MACH MUST BE INITIALIZED ACCORDING TO THE
  !     PROLOGUE IN EACH FUNCTION FOR THE MACHINE ENVIRONMENT BEFORE
  !     DQCKIN OR DBSKIN CAN BE EXECUTED.  FIFTEEN MACHINE ENVIRONMENTS
  !     CAN BE DEFINED IN I1MACH AND D1MACH.
  !
  !***ROUTINES CALLED  D1MACH, DBSKIN, I1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DQCKIN
  INTEGER i, ierr, iflg, ix, i1m12, j, k, kode, Lun, m, mdel, &
    mm, n, ndel, nn
  INTEGER I1MACH
  REAL(8) :: aix, er, tol, v, x, xinc, y
  REAL(8) :: D1MACH
  DIMENSION v(1), y(10)
  !***FIRST EXECUTABLE STATEMENT  DQCKIN
  tol = 1000.0D0*MAX(D1MACH(4),1.0D-18)
  iflg = 0
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  99001 FORMAT ('1'//' QUICK CHECK DIAGNOSTICS FOR DBSKIN'//)
  DO kode = 1, 2
    n = 0
    DO nn = 1, 7
      m = 1
      DO mm = 1, 4
        x = 0.0D0
        DO ix = 1, 6
          IF ( n/=0.OR.ix/=1 ) THEN
            CALL DBSKIN(x,n,kode,m,y,nz,ierr)
            DO k = 1, m, 2
              j = n + k - 1
              CALL DBSKIN(x,j,kode,1,v,nz,ierr)
              er = ABS((v(1)-y(k))/v(1))
              IF ( er>tol ) THEN
                IF ( iflg==0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Lun,99002)
                  99002                 FORMAT (8X,'X',13X,'V(1)',11X,'Y(K)',9X,'REL ER','R',5X,&
                    'KODE',3X,'N',4X,'M',4X,'K')
                ENDIF
                iflg = iflg + 1
                IF ( Kprint>=2 ) WRITE (Lun,99003) x, v(1), y(k), er, &
                  kode, n, m, k
                99003               FORMAT (4E15.6,4I5)
                IF ( iflg>200 ) GOTO 300
              ENDIF
            ENDDO
          ENDIF
          aix = 2*ix - 3
          xinc = MAX(1.0D0,aix)
          x = x + xinc
        ENDDO
        mdel = MAX(1,mm-1)
        m = m + mdel
      ENDDO
      ndel = MAX(1,2*n-2)
      n = n + ndel
    ENDDO
  ENDDO
  !-----------------------------------------------------------------------
  !     TEST UNDERFLOW
  !-----------------------------------------------------------------------
  kode = 1
  m = 10
  n = 10
  i1m12 = I1MACH(15)
  x = -2.302D0*D1MACH(5)*i1m12
  CALL DBSKIN(x,n,kode,m,y,nz,ierr)
  IF ( nz==m ) THEN
    DO i = 1, m
      IF ( y(i)/=0.0D0 ) GOTO 100
    ENDDO
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,99004)
    99004   FORMAT (//' NZ IN UNDERFLOW TEST IS NOT 1'//)
    iflg = iflg + 1
  ENDIF
  GOTO 200
  100  iflg = iflg + 1
  IF ( Kprint>=2 ) WRITE (Lun,99005)
  99005 FORMAT (//' SOME Y VALUE IN UNDERFLOW TEST IS NOT ZERO'//)
  200 CONTINUE
  IF ( iflg==0.AND.Kprint>=3 ) THEN
    WRITE (Lun,99006)
    99006   FORMAT (//' QUICK CHECKS OK'//)
  ENDIF
  Ipass = 0
  IF ( iflg==0 ) Ipass = 1
  RETURN
  300 CONTINUE
  IF ( Kprint>=2 ) WRITE (Lun,99007)
  99007 FORMAT (//' PROCESSING OF MAIN LOOPS TERMINATED BECAUSE THE NUM',&
    'BER OF DIAGNOSTIC PRINTS EXCEEDS 200'//)
  Ipass = 0
  IF ( iflg==0 ) Ipass = 1
END SUBROUTINE DQCKIN
