!*==SCHK43.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHK43
SUBROUTINE SCHK43(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
  IMPLICIT NONE
  !*--SCHK436
  !***BEGIN PROLOGUE  SCHK43
  !***SUBSIDIARY
  !***PURPOSE  Quick check for SSYRK.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !***DESCRIPTION
  !
  !  Quick check for SSYRK.
  !
  !  Auxiliary routine for test program for Level 3 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE3, SMMCH, SSYRK
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  SCHK43
  !     .. Parameters ..
  REAL ZERO , ONE
  PARAMETER (ZERO=0.0,ONE=1.0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL Eps , Thresh
  INTEGER Kprint , Nalf , Nbet , Nidim , Nmax , Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) , Bet(Nbet)&
    , G(Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bs(Nmax*Nmax) ,&
    C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
  INTEGER Idim(Nidim)
  !     .. Local Scalars ..
  REAL alpha , als , beta , bets , err , errmax
  INTEGER i , ia , ib , ict , icu , ik , in , j , jc , jj , k , laa , lcc ,&
    lda , ldas , ldc , ldcs , n , na , nargs , nc , nerr , ns , ks ,&
    lj , ma
  LOGICAL ftl , null , reset , tran , upper
  CHARACTER :: uplo , uplos , trans , transs
  CHARACTER(2) :: ichu
  CHARACTER(3) :: icht
  !     .. Local Arrays ..
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LSE , LSERES
  EXTERNAL LSE , LSERES , NUMXER
  !     .. External Subroutines ..
  EXTERNAL SSYRK , SMAKE3 , SMMCH
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !     .. Data statements ..
  DATA icht/'NTC'/ , ichu/'UL'/
  !***FIRST EXECUTABLE STATEMENT  SCHK43
  nargs = 10
  !
  nc = 0
  reset = .TRUE.
  errmax = ZERO
  !
  !
  DO in = 1 , Nidim
    n = Idim(in)
    !        Set LDC to 1 more than minimum value if room.
    ldc = n
    IF ( ldc<Nmax ) ldc = ldc + 1
    !        Skip tests if not enough room.
    IF ( ldc<=Nmax ) THEN
      lcc = ldc*n
      null = n<=0
      !
      DO ik = 1 , Nidim
        k = Idim(ik)
        !
        DO ict = 1 , 3
          trans = icht(ict:ict)
          tran = trans=='T' .OR. trans=='C'
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
            CALL SMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
            !
            DO icu = 1 , 2
              uplo = ichu(icu:icu)
              upper = uplo=='U'
              !
              DO ia = 1 , Nalf
                alpha = Alf(ia)
                !
                DO ib = 1 , Nbet
                  beta = Bet(ib)
                  !
                  !                       Generate the matrix C.
                  !
                  CALL SMAKE3('SY',uplo,' ',n,n,C,Nmax,Cc,ldc,reset,ZERO)
                  !
                  nc = nc + 1
                  !
                  !                       Save every datum before calling the subroutine.
                  !
                  uplos = uplo
                  transs = trans
                  ns = n
                  ks = k
                  als = alpha
                  DO i = 1 , laa
                    As(i) = Aa(i)
                  ENDDO
                  ldas = lda
                  bets = beta
                  DO i = 1 , lcc
                    Cs(i) = Cc(i)
                  ENDDO
                  ldcs = ldc
                  !
                  !                       Call the subroutine.
                  !
                  CALL SSYRK(uplo,trans,n,k,alpha,Aa,lda,beta,Cc,ldc)
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
                  isame(1) = uplos==uplo
                  isame(2) = transs==trans
                  isame(3) = ns==n
                  isame(4) = ks==k
                  isame(5) = als==alpha
                  isame(6) = LSE(As,Aa,laa)
                  isame(7) = ldas==lda
                  isame(8) = bets==beta
                  IF ( null ) THEN
                    isame(9) = LSE(Cs,Cc,lcc)
                  ELSE
                    isame(9) = LSERES('SY',uplo,n,n,Cs,Cc,ldc)
                  ENDIF
                  isame(10) = ldcs==ldc
                  !
                  !                       If data was incorrectly changed, report and
                  !                       return.
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
                    !
                    !                          Check the result column by column.
                    !
                    jc = 1
                    DO j = 1 , n
                      IF ( upper ) THEN
                        jj = 1
                        lj = j
                      ELSE
                        jj = j
                        lj = n - j + 1
                      ENDIF
                      IF ( tran ) THEN
                        CALL SMMCH('T','N',lj,1,k,alpha,A(1,jj),Nmax,A(1,j),&
                          Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                          Eps,err,ftl,Nout,.TRUE.,Kprint)
                      ELSE
                        CALL SMMCH('N','T',lj,1,k,alpha,A(jj,1),Nmax,A(j,1),&
                          Nmax,beta,C(jj,j),Nmax,Ct,G,Cc(jc),ldc,&
                          Eps,err,ftl,Nout,.TRUE.,Kprint)
                      ENDIF
                      IF ( upper ) THEN
                        jc = jc + ldc
                      ELSE
                        jc = jc + ldc + 1
                      ENDIF
                      errmax = MAX(errmax,err)
                    ENDDO
                  ENDIF
                  IF ( ftl ) THEN
                    Fatal = .TRUE.
                    IF ( Kprint>=3 ) THEN
                      WRITE (Nout,FMT=99004) Sname
                      WRITE (Nout,FMT=99005) nc , Sname , uplo , trans , n ,&
                        k , alpha , lda , beta , ldc
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
      ENDDO
    ENDIF
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
  RETURN
  !
  99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL','S)')
  99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
    'ANGED INCORRECTLY *******')
  99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
    /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
  99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
  99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,',I3,',',&
    F4.1,', C,',I3,')           .')
  99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of SCHK43.
  !
END SUBROUTINE SCHK43
