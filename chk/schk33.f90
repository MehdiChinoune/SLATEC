!*==SCHK33.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHK33
SUBROUTINE SCHK33(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Nmax,A,Aa,As,B,Bb,Bs,Ct,G,C)
  IMPLICIT NONE
  !*--SCHK336
  !***BEGIN PROLOGUE  SCHK33
  !***SUBSIDIARY
  !***PURPOSE  Quick check for STRMM and STRSM.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Quick check for STRMM and STRSM.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE3, SMMCH, STRMM, STRSM
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  SCHK33
  !     .. Parameters ..
  REAL ZERO , ONE
  PARAMETER (ZERO=0.0,ONE=1.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps , Thresh
  INTEGER Kprint , Nalf , Nidim , Nmax , Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) , G(Nmax) ,&
    B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bs(Nmax*Nmax) , C(Nmax,Nmax) ,&
    Ct(Nmax)
  INTEGER Idim(Nidim)
  !     .. Local Scalars ..
  REAL alpha , als , err , errmax
  INTEGER i , ia , icd , ics , ict , icu , im , in , j , laa , lbb , lda ,&
    ldas , ldb , ldbs , m , ms , n , na , nargs , nc , nerr , ns
  LOGICAL ftl , null , reset , left
  CHARACTER :: side , sides , uplo , uplos , tranas , transa , diag , diags
  CHARACTER(2) :: ichs , ichu , ichd
  CHARACTER(3) :: icht
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LSE , LSERES
  EXTERNAL LSE , LSERES , NUMXER
  !     .. External Subroutines ..
  EXTERNAL STRMM , STRSM , SMAKE3 , SMMCH
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !     .. Data statements ..
  DATA icht/'NTC'/ , ichu/'UL'/ , ichd/'UN'/ , ichs/'LR'/
  !***FIRST EXECUTABLE STATEMENT  SCHK33
  nargs = 11
  !
  nc = 0
  reset = .TRUE.
  errmax = ZERO
  !
  !     Set up zero matrix for SMMCH.
  DO j = 1 , Nmax
    DO i = 1 , Nmax
      C(i,j) = ZERO
    ENDDO
  ENDDO
  !
  DO im = 1 , Nidim
    m = Idim(im)
    !
    DO in = 1 , Nidim
      n = Idim(in)
      !           Set LDB to 1 more than minimum value if room.
      ldb = m
      IF ( ldb<Nmax ) ldb = ldb + 1
      !           Skip tests if not enough room.
      IF ( ldb<=Nmax ) THEN
        lbb = ldb*n
        null = m<=0 .OR. n<=0
        !
        DO ics = 1 , 2
          side = ichs(ics:ics)
          left = side=='L'
          IF ( left ) THEN
            na = m
          ELSE
            na = n
          ENDIF
          !              Set LDA to 1 more than minimum value if room.
          lda = na
          IF ( lda<Nmax ) lda = lda + 1
          !              Skip tests if not enough room.
          IF ( lda>Nmax ) EXIT
          laa = lda*na
          !
          DO icu = 1 , 2
            uplo = ichu(icu:icu)
            !
            DO ict = 1 , 3
              transa = icht(ict:ict)
              !
              DO icd = 1 , 2
                diag = ichd(icd:icd)
                !
                DO ia = 1 , Nalf
                  alpha = Alf(ia)
                  !
                  !                          Generate the matrix A.
                  !
                  CALL SMAKE3('TR',uplo,diag,na,na,A,Nmax,Aa,lda,reset,ZERO)
                  !
                  !                          Generate the matrix B.
                  !
                  CALL SMAKE3('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
                  !
                  nc = nc + 1
                  !
                  !                          Save every datum before calling the
                  !                          subroutine.
                  !
                  sides = side
                  uplos = uplo
                  tranas = transa
                  diags = diag
                  ms = m
                  ns = n
                  als = alpha
                  DO i = 1 , laa
                    As(i) = Aa(i)
                  ENDDO
                  ldas = lda
                  DO i = 1 , lbb
                    Bs(i) = Bb(i)
                  ENDDO
                  ldbs = ldb
                  !
                  !                          Call the subroutine.
                  !
                  IF ( Sname(4:5)=='MM' ) THEN
                    CALL STRMM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,&
                      ldb)
                  ELSEIF ( Sname(4:5)=='SM' ) THEN
                    CALL STRSM(side,uplo,transa,diag,m,n,alpha,Aa,lda,Bb,&
                      ldb)
                  ENDIF
                  !
                  !                          Check if error-exit was taken incorrectly.
                  !
                  IF ( NUMXER(nerr)/=0 ) THEN
                    IF ( Kprint>=2 ) WRITE (Nout,FMT=99006)
                    Fatal = .TRUE.
                  ENDIF
                  !
                  !                          See what data changed inside subroutines.
                  !
                  isame(1) = sides==side
                  isame(2) = uplos==uplo
                  isame(3) = tranas==transa
                  isame(4) = diags==diag
                  isame(5) = ms==m
                  isame(6) = ns==n
                  isame(7) = als==alpha
                  isame(8) = LSE(As,Aa,laa)
                  isame(9) = ldas==lda
                  IF ( null ) THEN
                    isame(10) = LSE(Bs,Bb,lbb)
                  ELSE
                    isame(10) = LSERES('GE',' ',m,n,Bs,Bb,ldb)
                  ENDIF
                  isame(11) = ldbs==ldb
                  !
                  !                          If data was incorrectly changed, report and
                  !                          return.
                  !
                  DO i = 1 , nargs
                    IF ( .NOT.isame(i) ) THEN
                      Fatal = .TRUE.
                      IF ( Kprint>=2 ) WRITE (Nout,FMT=99002) i
                    ENDIF
                  ENDDO
                  !
                  ftl = .FALSE.
                  IF ( .NOT.null ) THEN
                    IF ( Sname(4:5)=='MM' ) THEN
                      !
                      !                                Check the result.
                      !
                      IF ( left ) THEN
                        CALL SMMCH(transa,'N',m,n,m,alpha,A,Nmax,B,Nmax,&
                          ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,&
                          .TRUE.,Kprint)
                      ELSE
                        CALL SMMCH('N',transa,m,n,n,alpha,B,Nmax,A,Nmax,&
                          ZERO,C,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,&
                          .TRUE.,Kprint)
                      ENDIF
                    ELSEIF ( Sname(4:5)=='SM' ) THEN
                      !
                      !                                Compute approximation to original
                      !                                matrix.
                      !
                      DO j = 1 , n
                        DO i = 1 , m
                          C(i,j) = Bb(i+(j-1)*ldb)
                          Bb(i+(j-1)*ldb) = alpha*B(i,j)
                        ENDDO
                      ENDDO
                      !
                      IF ( left ) THEN
                        CALL SMMCH(transa,'N',m,n,m,ONE,A,Nmax,C,Nmax,ZERO,&
                          B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,&
                          .FALSE.,Kprint)
                      ELSE
                        CALL SMMCH('N',transa,m,n,n,ONE,C,Nmax,A,Nmax,ZERO,&
                          B,Nmax,Ct,G,Bb,ldb,Eps,err,ftl,Nout,&
                          .FALSE.,Kprint)
                      ENDIF
                    ENDIF
                    errmax = MAX(errmax,err)
                  ENDIF
                  IF ( ftl ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=3 ) THEN
                      WRITE (Nout,FMT=99004) Sname
                      WRITE (Nout,FMT=99005) nc , Sname , side , uplo ,&
                        transa , diag , m , n , alpha , lda , ldb
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
        WRITE (Nout,FMT=99001) Sname , nc
      ELSE
        WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
    ENDIF
  ENDIF
  !
  RETURN
  !
  99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
  99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
    'ANGED INCORRECTLY *******')
  99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
    /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
  99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
  99005 FORMAT (1X,I6,': ',A6,'(',4('''',A1,''','),2(I3,','),F4.1,', A,',I3,&
    ', B,',I3,')        .')
  99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of SCHK33.
  !
END SUBROUTINE SCHK33
