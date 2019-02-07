!*==CCHK43.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CCHK43
SUBROUTINE CCHK43(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
  IMPLICIT NONE
  !*--CCHK436
  !***BEGIN PROLOGUE  CCHK43
  !***SUBSIDIARY
  !***PURPOSE  Quick check for CHERK and CSYRK.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Quick check for CHERK and CSYRK.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CHERK, CMAKE3, CMMCH, CSYRK, LCE, LCERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CCHK43
  !     .. Parameters ..
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0,0.0))
  REAL RZERO, RONE
  PARAMETER (RZERO=0.0,RONE=1.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps, Thresh
  INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax) ,&
    B(Nmax,Nmax), Bb(Nmax*Nmax), Bs(Nmax*Nmax), Bet(Nbet) ,&
    C(Nmax,Nmax), Cc(Nmax*Nmax), Ct(Nmax), Cs(Nmax*Nmax)
  REAL G(Nmax)
  INTEGER Idim(Nidim)
  !     .. Local Scalars ..
  COMPLEX alpha, als, beta, bets
  REAL err, errmax, ralpha, rals, rbeta, rbets
  INTEGER i, ia, ib, ict, icu, ik, in, j, jc, jj, k, ks, laa ,&
    lcc, lda, ldas, ldc, ldcs, lj, ma, n, na, nargs, nc ,&
    nerr, ns
  LOGICAL conj, ftl, null, reset, tran, upper
  CHARACTER :: trans, transs, uplo, transt, uplos
  CHARACTER(2) :: ichu, icht
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LCE, LCERES
  EXTERNAL LCE, LCERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL CHERK, CSYRK, CMAKE3, CMMCH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, MIN
  !     .. Data statements ..
  DATA ichu/'UL'/, icht/'NC'/
  !***FIRST EXECUTABLE STATEMENT  CCHK43
  conj = Sname(2:3)=='HE'
  !
  nargs = 10
  nc = 0
  reset = .TRUE.
  errmax = RZERO
  !
  DO in = 1, Nidim
    n = Idim(in)
    !        Set LDC to 1 more than minimum value if room.
    ldc = n
    IF ( ldc<Nmax ) ldc = ldc + 1
    !        Skip tests if not enough room.
    IF ( ldc<=Nmax ) THEN
      lcc = ldc*n
      !
      DO ik = 1, Nidim
        k = Idim(ik)
        !
        DO ict = 1, 2
          trans = icht(ict:ict)
          tran = trans=='C'
          IF ( tran.AND..NOT.conj ) trans = 'T'
          IF ( tran ) THEN
            ma = k
            na = n
          ELSE
            ma = n
            na = k
          ENDIF
          !              Set LDA to 1 more than minimum value if room.
          lda = ma
          IF ( lda<Nmax ) lda = lda + 1
          !              Skip tests if not enough room.
          IF ( lda<=Nmax ) THEN
            laa = lda*na
            !
            !              Generate the matrix A.
            !
            CALL CMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
            !
            DO icu = 1, 2
              uplo = ichu(icu:icu)
              upper = uplo=='U'
              !
              DO ia = 1, Nalf
                alpha = Alf(ia)
                IF ( conj ) THEN
                  ralpha = REAL(alpha)
                  alpha = CMPLX(ralpha,RZERO)
                ENDIF
                !
                DO ib = 1, Nbet
                  beta = Bet(ib)
                  IF ( conj ) THEN
                    rbeta = REAL(beta)
                    beta = CMPLX(rbeta,RZERO)
                  ENDIF
                  null = n<=0
                  IF ( conj ) null = null .OR.&
                    ((k<=0.OR.ralpha==RZERO).AND.rbeta==RONE)
                  !
                  !                       Generate the matrix C.
                  !
                  CALL CMAKE3(Sname(2:3),uplo,' ',n,n,C,Nmax,Cc,ldc,reset,&
                    ZERO)
                  !
                  nc = nc + 1
                  !
                  !                       Save every datum before calling the subroutine.
                  !
                  uplos = uplo
                  transs = trans
                  ns = n
                  ks = k
                  IF ( conj ) THEN
                    rals = ralpha
                  ELSE
                    als = alpha
                  ENDIF
                  DO i = 1, laa
                    As(i) = Aa(i)
                  ENDDO
                  ldas = lda
                  IF ( conj ) THEN
                    rbets = rbeta
                  ELSE
                    bets = beta
                  ENDIF
                  DO i = 1, lcc
                    Cs(i) = Cc(i)
                  ENDDO
                  ldcs = ldc
                  !
                  !                       Call the subroutine.
                  !
                  IF ( conj ) THEN
                    CALL CHERK(uplo,trans,n,k,ralpha,Aa,lda,rbeta,Cc,ldc)
                  ELSE
                    CALL CSYRK(uplo,trans,n,k,alpha,Aa,lda,beta,Cc,ldc)
                  ENDIF
                  !
                  !                       Check if error-exit was taken incorrectly.
                  !
                  IF ( NUMXER(nerr)/=0 ) THEN
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                    Fatal = .TRUE.
                  ENDIF
                  !
                  !                       See what data changed inside subroutines.
                  !
                  isame(1) = uplos==uplo
                  isame(2) = transs==trans
                  isame(3) = ns==n
                  isame(4) = ks==k
                  IF ( conj ) THEN
                    isame(5) = rals==ralpha
                  ELSE
                    isame(5) = als==alpha
                  ENDIF
                  isame(6) = LCE(As,Aa,laa)
                  isame(7) = ldas==lda
                  IF ( conj ) THEN
                    isame(8) = rbets==rbeta
                  ELSE
                    isame(8) = bets==beta
                  ENDIF
                  IF ( null ) THEN
                    isame(9) = LCE(Cs,Cc,lcc)
                  ELSE
                    isame(9) = LCERES(Sname(2:3),uplo,n,n,Cs,Cc,ldc)
                  ENDIF
                  isame(10) = ldcs==ldc
                  !
                  !                       If data was incorrectly changed, report and
                  !                       return.
                  !
                  DO i = 1, nargs
                    IF ( .NOT.isame(i) ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                    ENDIF
                  ENDDO
                  !
                  IF ( .NOT.null ) THEN
                    !
                    !                          Check the result column by column.
                    !
                    IF ( conj ) THEN
                      transt = 'C'
                    ELSE
                      transt = 'T'
                    ENDIF
                    jc = 1
                    DO j = 1, n
                      IF ( upper ) THEN
                        jj = 1
                        lj = j
                      ELSE
                        jj = j
                        lj = n - j + 1
                      ENDIF
                      IF ( tran ) THEN
                        ftl = .FALSE.
                        CALL CMMCH(transt,'N',lj,1,k,alpha,A(1,jj),Nmax,&
                          A(1,j),Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc)&
                          ,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                      ELSE
                        ftl = .FALSE.
                        CALL CMMCH('N',transt,lj,1,k,alpha,A(jj,1),Nmax,&
                          A(j,1),Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc)&
                          ,ldc,Eps,err,ftl,Nout,.TRUE.,Kprint)
                      ENDIF
                      IF ( upper ) THEN
                        jc = jc + ldc
                      ELSE
                        jc = jc + ldc + 1
                      ENDIF
                      errmax = MAX(errmax,err)
                      IF ( ftl ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=3 ) THEN
                          WRITE (Nout,FMT=99004) Sname
                          IF ( conj ) THEN
                            WRITE (Nout,FMT=99005) nc, Sname, uplo ,&
                              trans, n, k, ralpha, lda, rbeta ,&
                              ldc
                          ELSE
                            WRITE (Nout,FMT=99006) nc, Sname, uplo ,&
                              trans, n, k, alpha, lda, beta, ldc
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF
                  !
                ENDDO
                !
              ENDDO
              !
            ENDDO
          ENDIF
          !
        ENDDO
        !
      ENDDO
    ENDIF
    !
  ENDDO
  !
  !     Report result.
  !
  IF ( .NOT.Fatal ) THEN
    IF ( Kprint>=3 ) THEN
      IF ( errmax<Thresh ) THEN
        WRITE (Nout,FMT=99001) Sname, nc
      ELSE
        WRITE (Nout,FMT=99003) Sname, nc, errmax
      ENDIF
    ENDIF
  ENDIF
  RETURN
  !
  99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
  99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
    'ANGED INCORRECTLY *******')
  99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
    /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
  99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
  99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,',',&
    F4.1,', C,',I3,')               ','          .')
  99006 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
    '), A,',I3,',(',F4.1,',',F4.1,'), C,',I3,')          .')
  99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of CCHK43.
  !
END SUBROUTINE CCHK43
