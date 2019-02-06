!*==DCHK12.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DCHK12
      SUBROUTINE DCHK12(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,
     &                  Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,
     &                  Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--DCHK127
!***BEGIN PROLOGUE  DCHK12
!***SUBSIDIARY
!***PURPOSE  Test DGEMV and DGBMV.
!***LIBRARY   SLATEC (BLAS)
!***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
!***AUTHOR  Du Croz, J. (NAG)
!           Hanson, R. J. (SNLA)
!***DESCRIPTION
!
!  Quick check for DGEMV and DGBMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DGBMV, DGEMV, DMAKE2, DMVCH, LDE, LDERES, NUMXER
!***REVISION HISTORY  (YYMMDD)
!   870810  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards. (BKS)
!***END PROLOGUE  DCHK12
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF
      PARAMETER (ZERO=0.0D0,HALF=0.5D0)
!     .. Scalar Arguments ..
      LOGICAL Fatal
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Kprint , Nalf , Nbet , Nidim , Ninc , Nkb , Nmax , Nout
      CHARACTER(6) :: Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax)
     &                 , Bet(Nbet) , G(Nmax) , X(Nmax) , Xs(Nmax*Incmax) ,
     &                 Xx(Nmax*Incmax) , Y(Nmax) , Ys(Nmax*Incmax) , Yt(Nmax) ,
     &                 Yy(Nmax*Incmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , beta , bls , err , errmax , transl
      INTEGER i , ia , ib , ic , iku , im , in , incx , incxs , incy , incys ,
     &        ix , iy , kl , kls , ku , kus , laa , lda , ldas , lx , ly , m ,
     &        ml , ms , n , nargs , nc , nd , nk , nl , ns , nerr
      LOGICAL banded , ftl , full , null , reset , tran
      CHARACTER :: trans , transs
      CHARACTER(3) :: ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      INTEGER NUMXER
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES , NUMXER
!     .. External Subroutines ..
      EXTERNAL DGBMV , DGEMV , DMAKE2 , DMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Data statements ..
      DATA ich/'NTC'/
!***FIRST EXECUTABLE STATEMENT  DCHK12
      full = Sname(3:3)=='E'
      banded = Sname(3:3)=='B'
!     Define the number of arguments.
      IF ( full ) THEN
        nargs = 11
      ELSEIF ( banded ) THEN
        nargs = 13
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
        n = Idim(in)
        nd = n/2 + 1
!
        DO im = 1 , 2
          IF ( im==1 ) m = MAX(n-nd,0)
          IF ( im==2 ) m = MIN(n+nd,Nmax)
!
          IF ( banded ) THEN
            nk = Nkb
          ELSE
            nk = 1
          ENDIF
          DO iku = 1 , nk
            IF ( banded ) THEN
              ku = Kb(iku)
              kl = MAX(ku-1,0)
            ELSE
              ku = n - 1
              kl = m - 1
            ENDIF
!              Set LDA to 1 more than minimum value if room.
            IF ( banded ) THEN
              lda = kl + ku + 1
            ELSE
              lda = m
            ENDIF
            IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
              laa = lda*n
              null = n<=0 .OR. m<=0
!
!              Generate the matrix A.
!
              transl = ZERO
              CALL DMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,kl,ku,reset,
     &                    transl)
!
              DO ic = 1 , 3
                trans = ich(ic:ic)
                tran = trans=='T' .OR. trans=='C'
!
                IF ( tran ) THEN
                  ml = n
                  nl = m
                ELSE
                  ml = m
                  nl = n
                ENDIF
!
                DO ix = 1 , Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*nl
!
!                    Generate the vector X.
!
                  transl = HALF
                  CALL DMAKE2('GE',' ',' ',1,nl,X,1,Xx,ABS(incx),0,nl-1,reset,
     &                        transl)
                  IF ( nl>1 ) THEN
                    X(nl/2) = ZERO
                    Xx(1+ABS(incx)*(nl/2-1)) = ZERO
                  ENDIF
!
                  DO iy = 1 , Ninc
                    incy = Inc(iy)
                    ly = ABS(incy)*ml
!
                    DO ia = 1 , Nalf
                      alpha = Alf(ia)
!
                      DO ib = 1 , Nbet
                        beta = Bet(ib)
!
!                             Generate the vector Y.
!
                        transl = ZERO
                        CALL DMAKE2('GE',' ',' ',1,ml,Y,1,Yy,ABS(incy),0,ml-1,
     &                              reset,transl)
!
                        nc = nc + 1
!
!                             Save every datum before calling the
!                             subroutine.
!
                        transs = trans
                        ms = m
                        ns = n
                        kls = kl
                        kus = ku
                        als = alpha
                        DO i = 1 , laa
                          As(i) = Aa(i)
                        ENDDO
                        ldas = lda
                        DO i = 1 , lx
                          Xs(i) = Xx(i)
                        ENDDO
                        incxs = incx
                        bls = beta
                        DO i = 1 , ly
                          Ys(i) = Yy(i)
                        ENDDO
                        incys = incy
!
!                             Call the subroutine.
!
                        IF ( full ) THEN
                          CALL DGEMV(trans,m,n,alpha,Aa,lda,Xx,incx,beta,Yy,
     &                               incy)
                        ELSEIF ( banded ) THEN
                          CALL DGBMV(trans,m,n,kl,ku,alpha,Aa,lda,Xx,incx,beta,
     &                               Yy,incy)
                        ENDIF
!
!                             Check if error-exit was taken incorrectly.
!
                        IF ( NUMXER(nerr)/=0 ) THEN
                          IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                          Fatal = .TRUE.
                        ENDIF
!
!                             See what data changed inside subroutines.
!
                        isame(1) = trans==transs
                        isame(2) = ms==m
                        isame(3) = ns==n
                        IF ( full ) THEN
                          isame(4) = als==alpha
                          isame(5) = LDE(As,Aa,laa)
                          isame(6) = ldas==lda
                          isame(7) = LDE(Xs,Xx,lx)
                          isame(8) = incxs==incx
                          isame(9) = bls==beta
                          IF ( null ) THEN
                            isame(10) = LDE(Ys,Yy,ly)
                          ELSE
                            isame(10) = LDERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                          ENDIF
                          isame(11) = incys==incy
                        ELSEIF ( banded ) THEN
                          isame(4) = kls==kl
                          isame(5) = kus==ku
                          isame(6) = als==alpha
                          isame(7) = LDE(As,Aa,laa)
                          isame(8) = ldas==lda
                          isame(9) = LDE(Xs,Xx,lx)
                          isame(10) = incxs==incx
                          isame(11) = bls==beta
                          IF ( null ) THEN
                            isame(12) = LDE(Ys,Yy,ly)
                          ELSE
                            isame(12) = LDERES('GE',' ',1,ml,Ys,Yy,ABS(incy))
                          ENDIF
                          isame(13) = incys==incy
                        ENDIF
!
!                             If data was incorrectly changed, report
!                             and return.
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
!                                Check the result.
!
                          CALL DMVCH(trans,m,n,alpha,A,Nmax,X,incx,beta,Y,incy,
     &                               Yt,G,Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                          errmax = MAX(errmax,err)
                        ENDIF
                        IF ( ftl ) THEN
                          Fatal = .TRUE.
                          IF ( Kprint>=3 ) THEN
                            WRITE (Nout,FMT=99004) Sname
                            IF ( full ) THEN
                              WRITE (Nout,FMT=99006) nc , Sname , trans , m ,
     &                               n , alpha , lda , incx , beta , incy
                            ELSEIF ( banded ) THEN
                              WRITE (Nout,FMT=99005) nc , Sname , trans , m ,
     &                               n , kl , ku , alpha , lda , incx , beta ,
     &                               incy
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
!
              ENDDO
            ENDIF
!
          ENDDO
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
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',4(I3,','),F4.1,', A,',I3,', X,',I2,
     &        ',',F4.1,', Y,',I2,') .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,', X,',I2,
     &        ',',F4.1,', Y,',I2,')         .')
99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     &        '******')
!
!     End of DCHK12.
!
      END SUBROUTINE DCHK12
