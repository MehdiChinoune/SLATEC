!DECK DCHK23
SUBROUTINE DCHK23(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCHK23
  !***SUBSIDIARY
  !***PURPOSE  Test DSYMM.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Quick check for DSYMM.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DMAKE3, DMMCH, DSYMM, LDE, LDERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
  !***END PROLOGUE  DCHK23
  !     .. Parameters ..
  REAL(8) :: ZERO
  PARAMETER (ZERO=0.0D0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL(8) :: Eps, Thresh
  INTEGER Kprint, Nalf, Nbet, Nidim, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax)&
    , G(Nmax), Bb(Nmax*Nmax), Bet(Nbet), Bs(Nmax*Nmax) ,&
    C(Nmax,Nmax), Cc(Nmax*Nmax), Cs(Nmax*Nmax), Ct(Nmax) ,&
    B(Nmax,Nmax)
  INTEGER Idim(Nidim)
  !     .. Local Scalars ..
  REAL(8) :: alpha, als, beta, bls, err, errmax
  INTEGER i, ia, ib, ics, icu, im, in, laa, lbb, lcc, lda, ldas ,&
    ldb, ldbs, ldc, ldcs, m, ms, n, na, nargs, nc, nerr, ns
  LOGICAL ftl, left, null, reset
  CHARACTER :: side, sides, uplo, uplos
  CHARACTER(2) :: ichs, ichu
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LDE, LDERES
  EXTERNAL LDE, LDERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL DSYMM, DMAKE3, DMMCH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, MIN
  !     .. Data statements ..
  DATA ichs/'LR'/, ichu/'UL'/
  !***FIRST EXECUTABLE STATEMENT  DCHK23
  nargs = 12
  nc = 0
  reset = .TRUE.
  errmax = ZERO
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
        !
        !           Set LDB to 1 more than minimum value if room.
        ldb = m
        IF ( ldb<Nmax ) ldb = ldb + 1
        !           Skip tests if not enough room.
        IF ( ldb<=Nmax ) THEN
          lbb = ldb*n
          !
          !           Generate the matrix B.
          !
          CALL DMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
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
                !                 Generate the symmetric matrix A.
                !
                CALL DMAKE3('SY',uplo,' ',na,na,A,Nmax,Aa,lda,reset,ZERO)
                !
                DO ia = 1, Nalf
                  alpha = Alf(ia)
                  !
                  DO ib = 1, Nbet
                    beta = Bet(ib)
                    !
                    !                       Generate the matrix C.
                    !
                    CALL DMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
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
                    CALL DSYMM(side,uplo,m,n,alpha,Aa,lda,Bb,ldb,beta,Cc,&
                      ldc)
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
                    isame(6) = LDE(As,Aa,laa)
                    isame(7) = ldas==lda
                    isame(8) = LDE(Bs,Bb,lbb)
                    isame(9) = ldbs==ldb
                    isame(10) = bls==beta
                    IF ( null ) THEN
                      isame(11) = LDE(Cs,Cc,lcc)
                    ELSE
                      isame(11) = LDERES('GE',' ',m,n,Cs,Cc,ldc)
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
                        CALL DMMCH('N','N',m,n,m,alpha,A,Nmax,B,Nmax,beta,C,&
                          Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,&
                          Kprint)
                      ELSE
                        CALL DMMCH('N','N',m,n,n,alpha,B,Nmax,A,Nmax,beta,C,&
                          Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,.TRUE.,&
                          Kprint)
                      ENDIF
                      errmax = MAX(errmax,err)
                    ENDIF
                    IF ( ftl ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=3 ) THEN
                        WRITE (Nout,FMT=99004) Sname
                        WRITE (Nout,FMT=99005) nc, Sname, side, uplo ,&
                          m, n, alpha, lda, ldb, beta, ldc
                      ENDIF
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
  99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,&
    ', B,',I3,',',F4.1,', C,',I3,')   ',' .')
  99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of DCHK23.
  !
END SUBROUTINE DCHK23
