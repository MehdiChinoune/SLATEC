!*==DCHK22.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DCHK22
      SUBROUTINE DCHK22(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nkb,Kb,
     &                  Nalf,Alf,Nbet,Bet,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,
     &                  Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--DCHK227
!***BEGIN PROLOGUE  DCHK22
!***SUBSIDIARY
!***PURPOSE  Test DSYMV, DSBMV and DSPMV.
!***LIBRARY   SLATEC (BLAS)
!***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
!***AUTHOR  Du Croz, J. (NAG)
!           Hanson, R. J. (SNLA)
!***DESCRIPTION
!
!  Quick check for DSYMV, DSBMV and DSPMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DMAKE2, DMVCH, DSBMV, DSPMV, DSYMV, LDE, LDERES,
!                    NUMXER
!***REVISION HISTORY  (YYMMDD)
!   870810  DATE WRITTEN
!   910619  Modified to meet SLATEC code and prologue standards. (BKS)
!***END PROLOGUE  DCHK22
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
      INTEGER i , ia , ib , ic , ik , in , incx , incxs , incy , incys , ix ,
     &        iy , k , ks , laa , lda , ldas , lx , ly , n , nargs , nc , nk ,
     &        ns , nerr
      LOGICAL banded , ftl , full , null , packed , reset
      CHARACTER :: uplo , uplos
      CHARACTER(2) :: ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      INTEGER NUMXER
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES , NUMXER
!     .. External Subroutines ..
      EXTERNAL DMAKE2 , DMVCH , DSBMV , DSPMV , DSYMV
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     .. Data statements ..
      DATA ich/'UL'/
!***FIRST EXECUTABLE STATEMENT  DCHK22
      full = Sname(3:3)=='Y'
      banded = Sname(3:3)=='B'
      packed = Sname(3:3)=='P'
!     Define the number of arguments.
      IF ( full ) THEN
        nargs = 10
      ELSEIF ( banded ) THEN
        nargs = 11
      ELSEIF ( packed ) THEN
        nargs = 9
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
        n = Idim(in)
!
        IF ( banded ) THEN
          nk = Nkb
        ELSE
          nk = 1
        ENDIF
        DO ik = 1 , nk
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
            DO ic = 1 , 2
              uplo = ich(ic:ic)
!
!              Generate the matrix A.
!
              transl = ZERO
              CALL DMAKE2(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,k,k,reset,
     &                    transl)
!
              DO ix = 1 , Ninc
                incx = Inc(ix)
                lx = ABS(incx)*n
!
!                 Generate the vector X.
!
                transl = HALF
                CALL DMAKE2('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,reset,
     &                      transl)
                IF ( n>1 ) THEN
                  X(n/2) = ZERO
                  Xx(1+ABS(incx)*(n/2-1)) = ZERO
                ENDIF
!
                DO iy = 1 , Ninc
                  incy = Inc(iy)
                  ly = ABS(incy)*n
!
                  DO ia = 1 , Nalf
                    alpha = Alf(ia)
!
                    DO ib = 1 , Nbet
                      beta = Bet(ib)
!
!                          Generate the vector Y.
!
                      transl = ZERO
                      CALL DMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,
     &                            transl)
!
                      nc = nc + 1
!
!                          Save every datum before calling the
!                          subroutine.
!
                      uplos = uplo
                      ns = n
                      ks = k
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
!                          Call the subroutine.
!
                      IF ( full ) THEN
                        CALL DSYMV(uplo,n,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                      ELSEIF ( banded ) THEN
                        CALL DSBMV(uplo,n,k,alpha,Aa,lda,Xx,incx,beta,Yy,incy)
                      ELSEIF ( packed ) THEN
                        CALL DSPMV(uplo,n,alpha,Aa,Xx,incx,beta,Yy,incy)
                      ENDIF
!
!                          Check if error-exit was taken incorrectly.
!
                      IF ( NUMXER(nerr)/=0 ) THEN
                        IF ( Kprint>=2 ) WRITE (Nout,FMT=99008)
                        Fatal = .TRUE.
                      ENDIF
!
!                          See what data changed inside subroutines.
!
                      isame(1) = uplo==uplos
                      isame(2) = ns==n
                      IF ( full ) THEN
                        isame(3) = als==alpha
                        isame(4) = LDE(As,Aa,laa)
                        isame(5) = ldas==lda
                        isame(6) = LDE(Xs,Xx,lx)
                        isame(7) = incxs==incx
                        isame(8) = bls==beta
                        IF ( null ) THEN
                          isame(9) = LDE(Ys,Yy,ly)
                        ELSE
                          isame(9) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                        ENDIF
                        isame(10) = incys==incy
                      ELSEIF ( banded ) THEN
                        isame(3) = ks==k
                        isame(4) = als==alpha
                        isame(5) = LDE(As,Aa,laa)
                        isame(6) = ldas==lda
                        isame(7) = LDE(Xs,Xx,lx)
                        isame(8) = incxs==incx
                        isame(9) = bls==beta
                        IF ( null ) THEN
                          isame(10) = LDE(Ys,Yy,ly)
                        ELSE
                          isame(10) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                        ENDIF
                        isame(11) = incys==incy
                      ELSEIF ( packed ) THEN
                        isame(3) = als==alpha
                        isame(4) = LDE(As,Aa,laa)
                        isame(5) = LDE(Xs,Xx,lx)
                        isame(6) = incxs==incx
                        isame(7) = bls==beta
                        IF ( null ) THEN
                          isame(8) = LDE(Ys,Yy,ly)
                        ELSE
                          isame(8) = LDERES('GE',' ',1,n,Ys,Yy,ABS(incy))
                        ENDIF
                        isame(9) = incys==incy
                      ENDIF
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
!
!                             Check the result.
!
                        CALL DMVCH('N',n,n,alpha,A,Nmax,X,incx,beta,Y,incy,Yt,G,
     &                             Yy,Eps,err,ftl,Nout,.TRUE.,Kprint)
                        errmax = MAX(errmax,err)
                      ENDIF
                      IF ( ftl ) THEN
                        Fatal = .TRUE.
                        IF ( Kprint>=3 ) THEN
                          WRITE (Nout,FMT=99004) Sname
                          IF ( full ) THEN
                            WRITE (Nout,FMT=99007) nc , Sname , uplo , n ,
     &                             alpha , lda , incx , beta , incy
                          ELSEIF ( banded ) THEN
                            WRITE (Nout,FMT=99006) nc , Sname , uplo , n ,
     &                             alpha , incx , beta , incy
                          ELSEIF ( packed ) THEN
                            WRITE (Nout,FMT=99005) nc , Sname , uplo , n ,
     &                             alpha , incx , beta , incy
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
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', AP',', X,',I2,',',
     &        F4.1,', Y,',I2,')                .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,', X,',I2,
     &        ',',F4.1,', Y,',I2,')         .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', A,',I3,', X,',I2,',',
     &        F4.1,', Y,',I2,')             .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     &        '******')
!
!     End of DCHK22.
!
      END SUBROUTINE DCHK22
