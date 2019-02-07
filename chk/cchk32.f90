!*==CCHK32.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CCHK32
SUBROUTINE CCHK32(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,&
    Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Xt,G,Z)
  IMPLICIT NONE
  !*--CCHK326
  !***BEGIN PROLOGUE  CCHK32
  !***SUBSIDIARY
  !***PURPOSE  Quick check for CTRMV, CTBMV, CTPMV, CTRSV, CTBSV and
  !            CTPSV.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Du Croz, J. (NAG)
  !           Hanson, R. J. (SNLA)
  !***DESCRIPTION
  !
  !  Quick check for CTRMV, CTBMV, CTPMV, CTRSV, CTBSV and CTPSV.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CMAKE2, CMVCH, CTBMV, CTBSV, CTPMV, CTPSV, CTRMV,
  !                    CTRSV, LCE, LCERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CCHK32
  !     .. Parameters ..
  COMPLEX ZERO, HALF, ONE
  PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0),ONE=(1.0,0.0))
  REAL RZERO
  PARAMETER (RZERO=0.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps, Thresh
  INTEGER Incmax, Kprint, Nidim, Ninc, Nkb, Nmax, Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  COMPLEX A(Nmax,Nmax), Aa(Nmax*Nmax), As(Nmax*Nmax), X(Nmax) ,&
    Xs(Nmax*Incmax), Xt(Nmax), Xx(Nmax*Incmax), Z(Nmax)
  REAL G(Nmax)
  INTEGER Idim(Nidim), Inc(Ninc), Kb(Nkb)
  !     .. Local Scalars ..
  COMPLEX transl
  REAL err, errmax
  INTEGER i, icd, ict, icu, ik, in, incx, incxs, ix, k, ks, laa ,&
    lda, ldas, lx, n, nargs, nc, nerr, nk, ns
  LOGICAL banded, ftl, full, null, packed, reset
  CHARACTER :: diag, diags, trans, transs, uplo, uplos
  CHARACTER(2) :: ichd, ichu
  CHARACTER(3) :: icht
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LCE, LCERES
  EXTERNAL LCE, LCERES, NUMXER
  !     .. External Subroutines ..
  EXTERNAL CMAKE2, CMVCH, CTBMV, CTBSV, CTPMV, CTPSV, CTRMV, CTRSV
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX
  !     .. Data statements ..
  DATA ichu/'UL'/, icht/'NTC'/, ichd/'UN'/
  !***FIRST EXECUTABLE STATEMENT  CCHK32
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
  errmax = RZERO
  !     Set up zero vector for CMVCH.
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
              CALL CMAKE2(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,lda,k,k,reset,&
                transl)
              !
              DO ix = 1, Ninc
                incx = Inc(ix)
                lx = ABS(incx)*n
                !
                !                       Generate the vector X.
                !
                transl = HALF
                CALL CMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,&
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
                    CALL CTRMV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                  ELSEIF ( banded ) THEN
                    CALL CTBMV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                  ELSEIF ( packed ) THEN
                    CALL CTPMV(uplo,trans,diag,n,Aa,Xx,incx)
                  ENDIF
                ELSEIF ( Sname(4:5)=='SV' ) THEN
                  IF ( full ) THEN
                    CALL CTRSV(uplo,trans,diag,n,Aa,lda,Xx,incx)
                  ELSEIF ( banded ) THEN
                    CALL CTBSV(uplo,trans,diag,n,k,Aa,lda,Xx,incx)
                  ELSEIF ( packed ) THEN
                    CALL CTPSV(uplo,trans,diag,n,Aa,Xx,incx)
                  ENDIF
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
                isame(1) = uplo==uplos
                isame(2) = trans==transs
                isame(3) = diag==diags
                isame(4) = ns==n
                IF ( full ) THEN
                  isame(5) = LCE(As,Aa,laa)
                  isame(6) = ldas==lda
                  IF ( null ) THEN
                    isame(7) = LCE(Xs,Xx,lx)
                  ELSE
                    isame(7) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                  ENDIF
                  isame(8) = incxs==incx
                ELSEIF ( banded ) THEN
                  isame(5) = ks==k
                  isame(6) = LCE(As,Aa,laa)
                  isame(7) = ldas==lda
                  IF ( null ) THEN
                    isame(8) = LCE(Xs,Xx,lx)
                  ELSE
                    isame(8) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
                  ENDIF
                  isame(9) = incxs==incx
                ELSEIF ( packed ) THEN
                  isame(5) = LCE(As,Aa,laa)
                  IF ( null ) THEN
                    isame(6) = LCE(Xs,Xx,lx)
                  ELSE
                    isame(6) = LCERES('GE',' ',1,n,Xs,Xx,ABS(incx))
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
                    CALL CMVCH(trans,n,n,ONE,A,Nmax,X,incx,ZERO,Z,incx,Xt,G,&
                      Xx,Eps,err,ftl,Nout,.TRUE.,Kprint)
                  ELSEIF ( Sname(4:5)=='SV' ) THEN
                    !
                    !                             Compute approximation to original vector.
                    !
                    DO i = 1, n
                      Z(i) = Xx(1+(i-1)*ABS(incx))
                      Xx(1+(i-1)*ABS(incx)) = X(i)
                    ENDDO
                    CALL CMVCH(trans,n,n,ONE,A,Nmax,Z,incx,ZERO,X,incx,Xt,G,&
                      Xx,Eps,err,ftl,Nout,.FALSE.,Kprint)
                  ENDIF
                  errmax = MAX(errmax,err)
                ENDIF
                IF ( ftl ) THEN
                  Fatal = .TRUE.
                  IF ( Kprint>=3 ) THEN
                    WRITE (Nout,FMT=99004) Sname
                    IF ( full ) THEN
                      WRITE (Nout,FMT=99006) nc, Sname, uplo, trans ,&
                        diag, n, lda, incx
                    ELSEIF ( banded ) THEN
                      WRITE (Nout,FMT=99005) nc, Sname, uplo, trans ,&
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
  99005 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', AP, ','X,',I2,&
    ')                                      .')
  99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,', X,',I2,&
    ')                               .')
  99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of CCHK32.
  !
END SUBROUTINE CCHK32
