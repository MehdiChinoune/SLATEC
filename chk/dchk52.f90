!DECK DCHK52
SUBROUTINE DCHK52(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCHK52
  !***SUBSIDIARY
  !***PURPOSE  Quick check for DSYR and DSPR.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Du Croz, J. (NAG)
  !           Hanson, R. J. (SNLA)
  !***DESCRIPTION
  !
  !  Quick check for DSYR and DSPR.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DMAKE2, DMVCH, DSPR, DSYR, LDE, LDERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
  !***END PROLOGUE  DCHK52
  !     .. Parameters ..
  REAL(8) :: ZERO, HALF, ONE
  PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL(8) :: Eps, Thresh
  INTEGER Incmax, Kprint, Nalf, Nidim, Ninc, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL(8) :: A(Nmax,Nmax), Aa(Nmax*Nmax), Alf(Nalf), As(Nmax*Nmax)&
    , G(Nmax), X(Nmax), Xs(Nmax*Incmax), Xx(Nmax*Incmax) ,&
    Y(Nmax), Ys(Nmax*Incmax), Yt(Nmax), Yy(Nmax*Incmax) ,&
    Z(Nmax)
  INTEGER Idim(Nidim), Inc(Ninc)
  !     .. Local Scalars ..
  REAL(8) :: alpha, als, err, errmax, transl
  INTEGER i, ia, ic, in, incx, incxs, ix, j, ja, jj, laa, lda ,&
    ldas, lj, lx, n, nargs, nc, ns, nerr
  LOGICAL ftl, full, null, packed, reset, upper
  CHARACTER :: uplo, uplos
  CHARACTER(2) :: ich
  !     .. Local Arrays ..
  REAL(8) :: w(1)
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LDE, LDERES
  EXTERNAL LDE, LDERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL DMAKE2, DMVCH, DSPR, DSYR
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX
  !     .. Data statements ..
  DATA ich/'UL'/
  !***FIRST EXECUTABLE STATEMENT  DCHK52
  full = Sname(3:3)=='Y'
  packed = Sname(3:3)=='P'
  !     Define the number of arguments.
  IF ( full ) THEN
    nargs = 7
  ELSEIF ( packed ) THEN
    nargs = 6
  ENDIF
  !
  nc = 0
  reset = .TRUE.
  errmax = ZERO
  !
  DO in = 1, Nidim
    n = Idim(in)
    !        Set LDA to 1 more than minimum value if room.
    lda = n
    IF ( lda<Nmax ) lda = lda + 1
    !        Skip tests if not enough room.
    IF ( lda<=Nmax ) THEN
      IF ( packed ) THEN
        laa = (n*(n+1))/2
      ELSE
        laa = lda*n
      ENDIF
      !
      DO ic = 1, 2
        uplo = ich(ic:ic)
        upper = uplo=='U'
        !
        DO ix = 1, Ninc
          incx = Inc(ix)
          lx = ABS(incx)*n
          !
          !              Generate the vector X.
          !
          transl = HALF
          CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,transl)
          IF ( n>1 ) THEN
            X(n/2) = ZERO
            Xx(1+ABS(incx)*(n/2-1)) = ZERO
          ENDIF
          !
          DO ia = 1, Nalf
            alpha = Alf(ia)
            null = n<=0 .OR. alpha==ZERO
            !
            !                 Generate the matrix A.
            !
            transl = ZERO
            CALL DMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,n-1,n-1,reset,&
              transl)
            !
            nc = nc + 1
            !
            !                 Save every datum before calling the subroutine.
            !
            uplos = uplo
            ns = n
            als = alpha
            DO i = 1, laa
              As(i) = Aa(i)
            ENDDO
            ldas = lda
            DO i = 1, lx
              Xs(i) = Xx(i)
            ENDDO
            incxs = incx
            !
            !                 Call the subroutine.
            !
            IF ( full ) THEN
              CALL DSYR(uplo,n,alpha,Xx,incx,Aa,lda)
            ELSEIF ( packed ) THEN
              CALL DSPR(uplo,n,alpha,Xx,incx,Aa)
            ENDIF
            !
            !                 Check if error-exit was taken incorrectly.
            !
            IF ( NUMXER(nerr)/=0 ) THEN
              IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
              Fatal = .TRUE.
            ENDIF
            !
            !                 See what data changed inside subroutines.
            !
            isame(1) = uplo==uplos
            isame(2) = ns==n
            isame(3) = als==alpha
            isame(4) = LDE(Xs,Xx,lx)
            isame(5) = incxs==incx
            IF ( null ) THEN
              isame(6) = LDE(As,Aa,laa)
            ELSE
              isame(6) = LDERES(Sname(2:3),uplo,n,n,As,Aa,lda)
            ENDIF
            IF ( .NOT.packed ) isame(7) = ldas==lda
            !
            !                 If data was incorrectly changed, report and return.
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
              !                    Check the result column by column.
              !
              IF ( incx>0 ) THEN
                DO i = 1, n
                  Z(i) = X(i)
                ENDDO
              ELSE
                DO i = 1, n
                  Z(i) = X(n-i+1)
                ENDDO
              ENDIF
              ja = 1
              DO j = 1, n
                w(1) = Z(j)
                IF ( upper ) THEN
                  jj = 1
                  lj = j
                ELSE
                  jj = j
                  lj = n - j + 1
                ENDIF
                ftl = .FALSE.
                CALL DMVCH('N',lj,1,alpha,Z(jj),lj,w,1,ONE,A(jj,j),1,Yt,G,&
                  Aa(ja),Eps,err,ftl,Nout,.TRUE.,Kprint)
                IF ( .NOT.(full) ) THEN
                  ja = ja + lj
                ELSEIF ( upper ) THEN
                  ja = ja + lda
                ELSE
                  ja = ja + lda + 1
                ENDIF
                errmax = MAX(errmax,err)
              ENDDO
            ENDIF
            IF ( ftl ) THEN
              Fatal = .TRUE.
              IF ( Kprint>=3 ) THEN
                WRITE (Nout,FMT=99005) j
                WRITE (Nout,FMT=99004) Sname
                IF ( full ) THEN
                  WRITE (Nout,FMT=99007) nc, Sname, uplo, n, alpha ,&
                    incx, lda
                ELSEIF ( packed ) THEN
                  WRITE (Nout,FMT=99006) nc, Sname, uplo, n, alpha ,&
                    incx
                ENDIF
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
  99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
  99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,&
    ', AP)                           .')
  99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', A,',I3,&
    ')                        .')
  99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of DCHK52.
  !
END SUBROUTINE DCHK52
