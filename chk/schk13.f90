!*==SCHK13.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SCHK13
      SUBROUTINE SCHK13(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,
     &                  Nbet,Bet,Nmax,A,Aa,As,B,Bb,Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--SCHK136
!***BEGIN PROLOGUE  SCHK13
!***SUBSIDIARY
!***PURPOSE  Quick check for SGEMM.
!***LIBRARY   SLATEC (BLAS)
!***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
!***AUTHOR  Dongarra, J. J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!***DESCRIPTION
!
!  Quick check for SGEMM.
!
!  Auxiliary routine for test program for Level 3 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  LSE, LSERES, NUMXER, SGEMM, SMAKE3, SMMCH
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  SCHK13
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0)
!     .. Scalar Arguments ..
      LOGICAL Fatal
      REAL Eps , Thresh
      INTEGER Kprint , Nalf , Nbet , Nidim , Nmax , Nout
      CHARACTER*6 Sname
!     .. Array Arguments ..
      REAL A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) , Bet(Nbet)
     &     , G(Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bs(Nmax*Nmax) , 
     &     C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      REAL alpha , als , beta , bls , err , errmax
      INTEGER i , ia , ib , ica , icb , ik , im , in , k , ks , laa , lbb , 
     &        lcc , lda , ldas , ldb , ldbs , ldc , ldcs , m , ma , mb , ms , 
     &        n , na , nargs , nb , nc , nerr , ns
      LOGICAL ftl , null , reset , trana , tranb
      CHARACTER*1 tranas , tranbs , transa , transb
      CHARACTER*3 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      INTEGER NUMXER
      LOGICAL LSE , LSERES
      EXTERNAL LSE , LSERES , NUMXER
!     .. External Subroutines ..
      EXTERNAL SGEMM , SMAKE3 , SMMCH
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Data statements ..
      DATA ich/'NTC'/
!***FIRST EXECUTABLE STATEMENT  SCHK13
      nargs = 13
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
!
      DO im = 1 , Nidim
        m = Idim(im)
!
        DO in = 1 , Nidim
          n = Idim(in)
!           Set LDC to 1 more than minimum value if room.
          ldc = m
          IF ( ldc<Nmax ) ldc = ldc + 1
!           Skip tests if not enough room.
          IF ( ldc<=Nmax ) THEN
            lcc = ldc*n
            null = n<=0 .OR. m<=0
!
            DO ik = 1 , Nidim
              k = Idim(ik)
!
              DO ica = 1 , 3
                transa = ich(ica:ica)
                trana = transa=='T' .OR. transa=='C'
!
                IF ( trana ) THEN
                  ma = k
                  na = m
                ELSE
                  ma = m
                  na = k
                ENDIF
!                 Set LDA to 1 more than minimum value if room.
                lda = ma
                IF ( lda<Nmax ) lda = lda + 1
!                 Skip tests if not enough room.
                IF ( lda<=Nmax ) THEN
                  laa = lda*na
!
!                 Generate the matrix A.
!
                  CALL SMAKE3('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset,ZERO)
!
                  DO icb = 1 , 3
                    transb = ich(icb:icb)
                    tranb = transb=='T' .OR. transb=='C'
!
                    IF ( tranb ) THEN
                      mb = n
                      nb = k
                    ELSE
                      mb = k
                      nb = n
                    ENDIF
!                    Set LDB to 1 more than minimum value if room.
                    ldb = mb
                    IF ( ldb<Nmax ) ldb = ldb + 1
!                    Skip tests if not enough room.
                    IF ( ldb<=Nmax ) THEN
                      lbb = ldb*nb
!
!                    Generate the matrix B.
!
                      CALL SMAKE3('GE',' ',' ',mb,nb,B,Nmax,Bb,ldb,reset,ZERO)
!
                      DO ia = 1 , Nalf
                        alpha = Alf(ia)
!
                        DO ib = 1 , Nbet
                          beta = Bet(ib)
!
!                          Generate the matrix C.
!
                          CALL SMAKE3('GE',' ',' ',m,n,C,Nmax,Cc,ldc,reset,ZERO)
!
                          nc = nc + 1
!
!                          Save every datum before calling the
!                          subroutine.
!
                          tranas = transa
                          tranbs = transb
                          ms = m
                          ns = n
                          ks = k
                          als = alpha
                          DO i = 1 , laa
                            As(i) = Aa(i)
                          ENDDO
                          ldas = lda
                          DO i = 1 , lbb
                            Bs(i) = Bb(i)
                          ENDDO
                          ldbs = ldb
                          bls = beta
                          DO i = 1 , lcc
                            Cs(i) = Cc(i)
                          ENDDO
                          ldcs = ldc
!
!                          Call the subroutine.
!
                          CALL SGEMM(transa,transb,m,n,k,alpha,Aa,lda,Bb,ldb,
     &                               beta,Cc,ldc)
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
                          isame(1) = transa==tranas
                          isame(2) = transb==tranbs
                          isame(3) = ms==m
                          isame(4) = ns==n
                          isame(5) = ks==k
                          isame(6) = als==alpha
                          isame(7) = LSE(As,Aa,laa)
                          isame(8) = ldas==lda
                          isame(9) = LSE(Bs,Bb,lbb)
                          isame(10) = ldbs==ldb
                          isame(11) = bls==beta
                          IF ( null ) THEN
                            isame(12) = LSE(Cs,Cc,lcc)
                          ELSE
                            isame(12) = LSERES('GE',' ',m,n,Cs,Cc,ldc)
                          ENDIF
                          isame(13) = ldcs==ldc
!
!                          If data was incorrectly changed, report
!                          and return.
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
!                             Check the result.
!
                            CALL SMMCH(transa,transb,m,n,k,alpha,A,Nmax,B,Nmax,
     &                                 beta,C,Nmax,Ct,G,Cc,ldc,Eps,err,ftl,Nout,
     &                                 .TRUE.,Kprint)
                            errmax = MAX(errmax,err)
                          ENDIF
                          IF ( ftl ) THEN
                            Fatal = .TRUE.
                            IF ( Kprint>=3 ) THEN
                              WRITE (Nout,FMT=99004) Sname
                              WRITE (Nout,FMT=99005) nc , Sname , transa , 
     &                               transb , m , n , k , alpha , lda , ldb , 
     &                               beta , ldc
                            ENDIF
                          ENDIF
!
                        ENDDO
!
                      ENDDO
                    ENDIF
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
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',
     &        /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',''',A1,''',',3(I3,','),F4.1,', A,',I3,
     &        ', B,',I3,',',F4.1,', ','C,',I3,').')
99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     &        '******')
!
!     End of SCHK13.
!
      END SUBROUTINE SCHK13
