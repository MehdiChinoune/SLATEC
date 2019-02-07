!*==CCHK23.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CCHK23
SUBROUTINE CCHK23(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
  IMPLICIT NONE
  !*--CCHK236
  !***BEGIN PROLOGUE  CCHK23
  !***SUBSIDIARY
  !***PURPOSE  Quick check for CHEMM and CSYMM.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Quick check for CHEMM and CSYMM.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CHEMM, CMAKE3, CMMCH, CSYMM, LCE, LCERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CCHK23
  !     .. Parameters ..
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0,0.0))
  REAL RZERO
  PARAMETER (RZERO=0.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps, Thresh
  INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax) ,&
    B(Nmax,Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax) ,&
    C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax)
  REAL G(Nmax)
  INTEGER Idim(Nidim)
  !     .. Local Scalars ..
  COMPLEX alpha, als, beta, bls
  REAL err, errmax
  INTEGER i, ia, ib, ics, icu, im, in, laa, lbb, lcc, lda, ldas ,&
    ldb, ldbs, ldc, ldcs, m, ms, n, na, nargs, nc, nerr, ns
  LOGICAL conj, ftl, left, null, reset
  CHARACTER :: side, sides, uplo, uplos
  CHARACTER(2) :: ichu, ichs
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LCE, LCERES
  EXTERNAL LCE, LCERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL CHEMM, CSYMM, CMAKE3, CMMCH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, MIN
  !     .. Data statements ..
  DATA ichs/'LR'/, ichu/'UL'/
  !***FIRST EXECUTABLE STATEMENT  CCHK23
  conj = Sname(2:3)=='HE'
  !
  nargs = 12
  nc = 0
  reset = .TRUE.
  errmax = RZERO
  !
  DO im = 1, Nidim
    m = Idim(im)
    !
    DO in = 1, Nidim
      n = Idim(in)
      !           Set LDC to 1 more than minimum value if room.
      ldc = m
      IF ( ldc<Nmax ) ldc = ldc + 1
      !           Skip tests if not enough room.
      IF ( ldc<=Nmax ) THEN
        lcc = ldc*n
        null = n<=0 .OR. m<=0
        !           Set LDB to 1 more than minimum value if room.
        ldb = m
        IF ( ldb<Nmax ) ldb = ldb + 1
        !           Skip tests if not enough room.
        IF ( ldb<=Nmax ) THEN
          lbb = ldb*n
          !
          !           Generate the matrix B.
          !
          CALL CMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
          !
          DO ics = 1, 2
            side = ichs(ics:ics)
            left = side=='L'
            !
            IF ( left ) THEN
              na = m
            ELSE
              na = n
            ENDIF
            !              Set LDA to 1 more than minimum value if room.
            lda = na
            IF ( lda<Nmax ) lda = lda + 1
            !              Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
              laa = lda*na
              !
              DO icu = 1, 2
                uplo = ichu(icu:icu)
                !
                !                 Generate the hermitian or symmetric matrix A.
                !
                CALL CMAKE3(Sname(2:3),uplo,' ',na,na,A,Nmax,Aa,lda,reset,&
                  ZERO)
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL CMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
                    !
                    nc = nc + 1
                    !
                    !                       Save every datum before calling the
                    !                       subroutine.
                    !
                    sides = side
                    uplos = uplo
                    ms = m
                    ns = n
                    als = alpha
                    DO i = 1, laa
                      As(i) = Aa(i)
                    ENDDO
                    ldas = lda
                    DO i = 1, lbb
                      Bs(i) = Bb(i)
                    ENDDO
                    ldbs = ldb
                    bls = beta
                    DO i = 1, lcc
                      Cs(i) = Cc(i)
                    ENDDO
                    ldcs = ldc
                    !
                    !                       Call the subroutine.
                    !
                    IF ( conj ) THEN
                      CALL CHEMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,&
                        ldc)
                    ELSE
                      CALL CSYMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,&
                        ldc)
                    ENDIF
                    !
                    !                       Check if error-exit was taken incorrectly.
                    !
                    IF ( NUMXER(nerr)/=0 ) THEN
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                      Fatal = .TRUE.
                    ENDIF
                    !
                    !                       See what data changed inside subroutines.
                    !
                    isame(1) = sides==side
                    isame(2) = uplos==uplo
                    isame(3) = ms==m
                    isame(4) = ns==n
                    isame(5) = als==alpha
                    isame(6) = LCE(As,Aa,laa)
                    isame(7) = ldas==lda
                    isame(8) = LCE(Bs,Bb,lbb)
                    isame(9) = ldbs==ldb
                    isame(10) = bls==beta
                    IF ( null ) THEN
                      isame(11) = LCE(Cs,Cc,lcc)
                    ELSE
                      isame(11) = LCERES('GE',' ',m,n,Cs,Cc,ldc)
                    ENDIF
                    isame(12) = ldcs==ldc
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
                    ftl = .FALSE.
                    IF ( .NOT.null ) THEN
                      !
                      !                          Check the result.
                      !
                      IF ( left ) THEN
                        CALL CMMCH('N','N',m,n,m,alpha,A,Nmax,B,Nmax,beta,C,&
                          Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,&
                          Kprint)
                      ELSE
                        CALL CMMCH('N','N',m,n,n,alpha,B,Nmax,A,Nmax,beta,C,&
                          Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,&
                          Kprint)
                      ENDIF
                      errmax = MAX(errmax,err)
                    ENDIF
                    !
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, side, uplo ,&
                          m, n, alpha, lda, ldb, beta, ldc
                      ENDIF
                    ENDIF
                  ENDDO
                  !
                ENDDO
                !
              ENDDO
            ENDIF
            !
          ENDDO
        ENDIF
      ENDIF
      !
    ENDDO
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
  99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',F4.1,&
    '), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,')    .')
  99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of CCHK23.
  !
END SUBROUTINE CCHK23
