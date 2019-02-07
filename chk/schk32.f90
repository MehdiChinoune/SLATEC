!*==SCHK32.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHK32
SUBROUTINE SCHK32(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,&
    Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Xt,G,Z)
  IMPLICIT NONE
  !*--SCHK326
  !***BEGIN PROLOGUE  SCHK32
  !***SUBSIDIARY
  !***PURPOSE  Quick check for STRMV, STBMV, STPMV, STRSV, STBSV and
  !            STPSV.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Du Croz, J. (NAG)
  !           Hanson, R. J. (SNLA)
  !***DESCRIPTION
  !
  !  Quick check for STRMV, STBMV, STPMV, STRSV, STBSV and STPSV.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE2, SMVCH, STBMV, STBSV,
  !                    STPMV, STPSV, STRMV, STRSV
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  SCHK32
  !     .. Parameters ..
  REAL ZERO, HALF, ONE
  PARAMETER (ZERO=0.0,HALF=0.5,ONE=1.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps, Thresh
  INTEGER Incmax, Kprint, Nidim, Ninc, Nkb, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL A(Nmax,Nmax), Aa(Nmax*Nmax), As(Nmax*Nmax), G(Nmax), X(Nmax) ,&
    Xs(Nmax*Incmax), Xt(Nmax), Xx(Nmax*Incmax), Z(Nmax)
  INTEGER Idim(Nidim), Inc(Ninc), Kb(Nkb)
  !     .. Local Scalars ..
  REAL err, errmax, transl
  INTEGER i, icd, ict, icu, ik, in, incx, incxs, ix, k, ks, laa ,&
    lda, ldas, lx, n, nargs, nc, nk, ns, nerr
  LOGICAL banded, ftl, full, null, packed, reset
  CHARACTER :: diag, diags, trans, transs, uplo, uplos
  CHARACTER(2) :: ichd, ichu
  CHARACTER(3) :: icht
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LSE, LSERES
  EXTERNAL LSE, LSERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL SMAKE2, SMVCH, STBMV, STBSV, STPMV, STPSV, STRMV, STRSV
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX
  !     .. Data statements ..
  DATA ichu/'UL'/, icht/'NTC'/, ichd/'UN'/
  !***FIRST EXECUTABLE STATEMENT  SCHK32
  full = Sname(3:3)=='R'
  banded = Sname(3:3)=='B'
  packed = Sname(3:3)=='P'
  !     Define the number of arguments.
  IF ( full ) THEN
    nargs = 8
  ELSEIF ( banded ) THEN
    nargs = 9
  ELSEIF ( packed ) THEN
    nargs = 7
  ENDIF
  !
  nc = 0
  reset = .TRUE.
  errmax = ZERO
  !     Set up zero vector for SMVCH.
  DO i = 1, Nmax
    Z(i) = ZERO
  ENDDO
  !
  DO in = 1, Nidim
    n = Idim(in)
    !
    IF ( banded ) THEN
      nk = Nkb
    ELSE
      nk = 1
    ENDIF
    DO ik = 1, nk
      IF ( banded ) THEN
        k = Kb(ik)
      ELSE
        k = n - 1
      ENDIF
      !           Set LDA to 1 more than minimum value if room.
      IF ( banded ) THEN
        lda = k + 1
      ELSE
        lda = n
      ENDIF
      IF ( lda<Nmax ) lda = lda + 1
      !           Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        IF ( packed ) THEN
          laa = (n*(n+1))/2
        ELSE
          laa = lda*n
        ENDIF
        null = n<=0
        !
        DO icu = 1, 2
          uplo = ichu(icu:icu)
          !
          DO ict = 1, 3
            trans = icht(ict:ict)
            !
            DO icd = 1, 2
              diag = ichd(icd:icd)
              !
              !                    Generate the matrix A.
              !
              transl = ZERO
              CALL SMAKE2(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,lda,k,k,reset,&
                transl)
              !
              DO ix = 1, Ninc
                incx = Inc(ix)
                lx = ABS(incx)*n
                !
                !                       Generate the vector X.
                !
                transl = HALF
                CALL SMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,&
                  transl)
                IF ( n>1 ) THEN
                  X(n/2) = ZERO
                  Xx(1+ABS(incx)*(n/2-1)) = ZERO
                ENDIF
                !
                nc = nc + 1
                !
                !                       Save every datum before calling the subroutine.
                !
                uplos = uplo
                transs = trans
                diags = diag
                ns = n
                ks = k
                DO i = 1, laa
                  As(i) = Aa(i)
                ENDDO
                ldas = lda
                DO i = 1, lx
                  Xs(i) = Xx(i)
                ENDDO
                incxs = incx
                !
                !                       Call the subroutine.
                !
                IF ( Sname(4:5)=='MV' ) THEN
                  IF ( full ) THEN
                    CALL STRMV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                  ELSEIF ( banded ) THEN
                    CALL STBMV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                  ELSEIF ( packed ) THEN
                    CALL STPMV(uplo,trans,diag,n,Aa,Xx,incx)
                  ENDIF
                ELSEIF ( Sname(4:5)=='SV' ) THEN
                  IF ( full ) THEN
                    CALL STRSV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                  ELSEIF ( banded ) THEN
                    CALL STBSV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                  ELSEIF ( packed ) THEN
                    CALL STPSV(uplo,trans,diag,n,Aa,Xx,incx)
                  ENDIF
                ENDIF
                !
                !                       Check if error-exit was taken incorrectly.
                !
                IF ( NUMXER(nerr)/=0 ) THEN
                  IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                  Fatal = .TRUE.
                ENDIF
                !
                !                       See what data changed inside subroutines.
                !
                isame(1) = uplo==uplos
                isame(2) = trans==transs
                isame(3) = diag==diags
                isame(4) = ns==n
                IF ( full ) THEN
                  isame(5) = LSE(As,Aa,laa)
                  isame(6) = ldas==lda
                  IF ( null ) THEN
                    isame(7) = LSE(Xs,Xx,lx)
                  ELSE
                    isame(7) = LSERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                  ENDIF
                  isame(8) = incxs==incx
                ELSEIF ( banded ) THEN
                  isame(5) = ks==k
                  isame(6) = LSE(As,Aa,laa)
                  isame(7) = ldas==lda
                  IF ( null ) THEN
                    isame(8) = LSE(Xs,Xx,lx)
                  ELSE
                    isame(8) = LSERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                  ENDIF
                  isame(9) = incxs==incx
                ELSEIF ( packed ) THEN
                  isame(5) = LSE(As,Aa,laa)
                  IF ( null ) THEN
                    isame(6) = LSE(Xs,Xx,lx)
                  ELSE
                    isame(6) = LSERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                  ENDIF
                  isame(7) = incxs==incx
                ENDIF
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
                  IF ( Sname(4:5)=='MV' ) THEN
                    !
                    !                             Check the result.
                    !
                    CALL SMVCH(trans,n,n,ONE,A,Nmax,X,incx,ZERO,Z,incx,Xt,G,&
                      Xx,Eps,err,ftl,Nout,.TRUE.,Kprint)
                  ELSEIF ( Sname(4:5)=='SV' ) THEN
                    !
                    !                             Compute approximation to original vector.
                    !
                    DO i = 1, n
                      Z(i) = Xx(1+(i-1)*ABS(incx))
                      Xx(1+(i-1)*ABS(incx)) = X(i)
                    ENDDO
                    CALL SMVCH(trans,n,n,ONE,A,Nmax,Z,incx,ZERO,X,incx,Xt,G,&
                      Xx,Eps,err,ftl,Nout,.FALSE.,Kprint)
                  ENDIF
                  errmax = MAX(errmax,err)
                ENDIF
                IF ( ftl ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=3 ) THEN
                    WRITE (Nout,FMT=99004) Sname
                    IF ( full ) THEN
                      WRITE (Nout,FMT=99007) nc, Sname, uplo, trans ,&
                        diag, n, lda, incx
                    ELSEIF ( banded ) THEN
                      WRITE (Nout,FMT=99006) nc, Sname, uplo, trans ,&
                        diag, n, k, lda, incx
                    ELSEIF ( packed ) THEN
                      WRITE (Nout,FMT=99005) nc, Sname, uplo, trans ,&
                        diag, n, incx
                    ENDIF
                  ENDIF
                ENDIF
                !
              ENDDO
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
  !
  !     Report result.
  !
  IF ( .NOT.(Fatal) ) THEN
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
  99005 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', AP, ','X,',I2,&
    ')                        .')
  99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,', X,',I2,&
    ')                 .')
  99007 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', A,',I3,', X,',I2,&
    ')                     .')
  99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of SCHK32.
  !
END SUBROUTINE SCHK32
