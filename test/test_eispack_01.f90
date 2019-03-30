MODULE TEST24_MOD
  IMPLICIT NONE

CONTAINS
  !** EISQX1
  SUBROUTINE EISQX1(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for SGEEV and CGEEV.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Description:**
    !
    !     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
    !     SGEEV AND CGEEV.  THE EIGENVALUES OF INPUT MATRIX A(.,.)
    !     ARE STORED IN EK(.).  RELERR IS THE RELATIVE ACCURACY
    !     REQUIRED FOR THEM TO PASS.
    !
    !***
    ! **Routines called:**  CGEEV, R1MACH, SGEEV

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900405  CALL to XERROR replaced by message to LUN.  (WRB)

    INTEGER info
    REAL R1MACH
    INTEGER Kprint, Ipass, Lun
    INTEGER job, i, j, id
    REAL w(9)
    REAL err, erri, relerr, recj
    COMPLEX ac(3,3), ec(3), vc(3,3)
    INTEGER, PARAMETER :: lda = 3, n = 3, ldv = 3
    REAL :: a(3,3) = RESHAPE( [ 1., -2., 6., -1., 0., -3., 2., 5., 6. ], [3,3] )
    REAL, PARAMETER :: ek(3) = [ -1., 3., 5. ]
    !* FIRST EXECUTABLE STATEMENT  EISQX1
    Ipass = 1
    relerr = SQRT(R1MACH(4))
    DO j = 1, n
      DO i = 1, n
        ac(i,j) = CMPLX(a(i,j),0.)
      ENDDO
    ENDDO
    job = 1
    CALL CGEEV(ac,lda,n,ec,vc,ldv,w,job,info)
    IF ( info/=0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99003) 'CGEEV', info
      Ipass = 0
    ENDIF
    DO j = 1, n
      err = ABS(AIMAG(ec(j)))
      IF ( err>=relerr ) Ipass = 0
      recj = REAL(ec(j))
      err = ABS(recj-ek(1))
      id = 1
      DO i = 2, n
        erri = ABS(recj-ek(i))
        IF ( erri<err ) id = i
        err = MIN(erri,err)
      ENDDO
      IF ( ABS(recj-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
    ENDDO
    job = 0
    CALL SGEEV(a,lda,n,ec,vc,ldv,w,job,info)
    IF ( info/=0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99003) 'SGEEV', info
      Ipass = 0
    ENDIF
    DO j = 1, n
      err = ABS(AIMAG(ec(j)))
      IF ( err>=relerr ) Ipass = 0
      recj = REAL(ec(j))
      err = ABS(recj-ek(1))
      id = 1
      DO i = 2, n
        erri = ABS(recj-ek(i))
        IF ( erri<err ) id = i
        err = MIN(erri,err)
      ENDDO
      IF ( ABS(recj-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
    ENDDO
    IF ( Kprint>=2.AND.Ipass/=0 ) WRITE (Lun,99001)
    99001 FORMAT (' EISQX1 PASSES ALL TESTS.')
    IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99002)
    99002 FORMAT (' EISQX1 FAILS SOME TESTS.')
    99003 FORMAT (1X,'Eigenvalue iteration failed to converge in ',A5,', INFO = ',&
      I4)
  END SUBROUTINE EISQX1
  !** EISQX2
  SUBROUTINE EISQX2(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for SSIEV, CHIEV and SSPEV.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Kahaner, D. K., (NBS)
    !***
    ! **Description:**
    !
    !     THIS QUICK CHECK ROUTINE IS WRITTEN FOR EISPACK DRIVERS
    !     SSIEV, CHIEV AND SSPEV.  THE EIGENVALUES OF INPUT MATRIX
    !     A(.,.) ARE STORED IN EK(.).  RELERR IS THE RELATIVE
    !     ACCURACY REQUIRED FOR THEM TO PASS.
    !
    !***
    ! **Routines called:**  CHIEV, R1MACH, SSIEV, SSPEV

    !* REVISION HISTORY  (YYMMDD)
    !   800808  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   900405  CALL to XERROR replaced by message to LUN.  (WRB)

    INTEGER info
    REAL R1MACH
    INTEGER Kprint, Ipass, Lun
    INTEGER job, i, j, id
    REAL a1(4,4), a2(10), e(4), v(4,4), w(16)
    REAL err, erri, relerr
    COMPLEX ac(4,4), vc(4,4)
    INTEGER, PARAMETER :: lda = 4, n = 4, ldv = 4
    REAL, PARAMETER :: ap(10) = [ 5., 4., 5., 1., 1., 4., 1., 1., 2., 4. ]
    REAL, PARAMETER :: ek(4) = [ 1., 2., 5., 10. ]
    !* FIRST EXECUTABLE STATEMENT  EISQX2
    Ipass = 1
    relerr = SQRT(R1MACH(4))
    id = 0
    DO j = 1, n
      DO i = 1, j
        id = id + 1
        a1(i,j) = ap(id)
        a2(id) = ap(id)
        ac(i,j) = CMPLX(ap(id),0.)
      ENDDO
    ENDDO
    job = 1
    CALL CHIEV(ac,lda,n,e,vc,ldv,w,job,info)
    IF ( info/=0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99003) 'CHIEV', info
      Ipass = 0
    ENDIF
    DO j = 1, n
      err = ABS(e(j)-ek(1))
      id = 1
      DO i = 2, n
        erri = ABS(e(j)-ek(i))
        IF ( erri<err ) id = i
        err = MIN(erri,err)
      ENDDO
      IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
    ENDDO
    CALL SSIEV(a1,lda,n,e,w,job,info)
    IF ( info/=0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99003) 'SSIEV', info
      Ipass = 0
    ENDIF
    DO j = 1, n
      err = ABS(e(j)-ek(1))
      id = 1
      DO i = 2, n
        erri = ABS(e(j)-ek(i))
        IF ( erri<err ) id = i
        err = MIN(erri,err)
      ENDDO
      IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
    ENDDO
    job = 0
    CALL SSPEV(a2,n,e,v,ldv,w,job,info)
    IF ( info/=0 ) THEN
      IF ( Kprint>=2 ) WRITE (Lun,99003) 'SSPEV', info
      Ipass = 0
    ENDIF
    DO j = 1, n
      err = ABS(e(j)-ek(1))
      id = 1
      DO i = 2, n
        erri = ABS(e(j)-ek(i))
        IF ( erri<err ) id = i
        err = MIN(erri,err)
      ENDDO
      IF ( ABS(e(j)-ek(id))/ABS(ek(id))>=relerr ) Ipass = 0
    ENDDO
    IF ( Kprint>=2.AND.Ipass/=0 ) WRITE (Lun,99001)
    99001 FORMAT (' EISQX2 PASSES ALL TESTS.')
    IF ( Kprint>=1.AND.Ipass==0 ) WRITE (Lun,99002)
    99002 FORMAT (' EISQX2 FAILS SOME TESTS.')
    99003 FORMAT (1X,'Eigenvalue iteration failed to converge in ',A5,', INFO = ',&
      I4)
  END SUBROUTINE EISQX2
END MODULE TEST24_MOD
!** TEST24
PROGRAM TEST24
  USE TEST24_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2
  !***
  ! **Type:**      ALL (TEST24-A)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        SGEEV    CGEEV
  !        SSIEV    CHIEV
  !        SSPEV
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  EISQX1, EISQX2, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST24
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test SGEEV and CGEEV
  !
  CALL EISQX1(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test SSIEV, CHIEV and SSPEV
  !
  CALL EISQX2(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST24 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST24 *************')
  ENDIF
  STOP
END PROGRAM TEST24
