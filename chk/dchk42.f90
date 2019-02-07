!*==DCHK42.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DCHK42
SUBROUTINE DCHK42(Sname,Eps,Thresh,Nout,Kprint,Fatal,Nidim,Idim,Nalf,Alf,&
    Ninc,Inc,Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
  IMPLICIT NONE
  !*--DCHK426
  !***BEGIN PROLOGUE  DCHK42
  !***SUBSIDIARY
  !***PURPOSE  Test DGER.
  !***LIBRARY   SLATEC (BLAS)
  !***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
  !***AUTHOR  Du Croz, J. (NAG)
  !           Hanson, R. J. (SNLA)
  !***DESCRIPTION
  !
  !  Quick check for DGER.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DGER, DMAKE2, DMVCH, LDE, LDERES, NUMXER
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards. (BKS)
  !***END PROLOGUE  DCHK42
  !     .. Parameters ..
  REAL(8) :: ZERO , HALF , ONE
  PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
  !     .. Scalar Arguments ..
  LOGICAL Fatal
  REAL(8) :: Eps , Thresh
  INTEGER Incmax , Kprint , Nalf , Nidim , Ninc , Nmax , Nout
  CHARACTER(6) :: Sname
  !     .. Array Arguments ..
  REAL(8) :: A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax)&
    , G(Nmax) , X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) ,&
    Y(Nmax) , Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) ,&
    Z(Nmax)
  INTEGER Idim(Nidim) , Inc(Ninc)
  !     .. Local Scalars ..
  REAL(8) :: alpha , als , err , errmax , transl
  INTEGER i , ia , im , in , incx , incxs , incy , incys , ix , iy , j ,&
    laa , lda , ldas , lx , ly , m , ms , n , nargs , nc , nd , ns ,&
    nerr
  LOGICAL ftl , null , reset
  !     .. Local Arrays ..
  REAL(8) :: w(1)
  LOGICAL isame(13)
  !     .. External Functions ..
  INTEGER NUMXER
  LOGICAL LDE , LDERES
  EXTERNAL LDE , LDERES , NUMXER
  !     .. External Subroutines ..
  EXTERNAL DGER , DMAKE2 , DMVCH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , MAX , MIN
  !***FIRST EXECUTABLE STATEMENT  DCHK42
  !     Define the number of arguments.
  nargs = 9
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
      !           Set LDA to 1 more than minimum value if room.
      lda = m
      IF ( lda<Nmax ) lda = lda + 1
      !           Skip tests if not enough room.
      IF ( lda<=Nmax ) THEN
        laa = lda*n
        null = n<=0 .OR. m<=0
        !
        DO ix = 1 , Ninc
          incx = Inc(ix)
          lx = ABS(incx)*m
          !
          !              Generate the vector X.
          !
          transl = HALF
          CALL DMAKE2('GE',' ',' ',1,m,X,1,Xx,ABS(incx),0,m-1,reset,transl)
          IF ( m>1 ) THEN
            X(m/2) = ZERO
            Xx(1+ABS(incx)*(m/2-1)) = ZERO
          ENDIF
          !
          DO iy = 1 , Ninc
            incy = Inc(iy)
            ly = ABS(incy)*n
            !
            !                 Generate the vector Y.
            !
            transl = ZERO
            CALL DMAKE2('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,reset,&
              transl)
            IF ( n>1 ) THEN
              Y(n/2) = ZERO
              Yy(1+ABS(incy)*(n/2-1)) = ZERO
            ENDIF
            !
            DO ia = 1 , Nalf
              alpha = Alf(ia)
              !
              !                    Generate the matrix A.
              !
              transl = ZERO
              CALL DMAKE2(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,m-1,n-1,&
                reset,transl)
              !
              nc = nc + 1
              !
              !                    Save every datum before calling the subroutine.
              !
              ms = m
              ns = n
              als = alpha
              DO i = 1 , laa
                As(i) = Aa(i)
              ENDDO
              ldas = lda
              DO i = 1 , lx
                Xs(i) = Xx(i)
              ENDDO
              incxs = incx
              DO i = 1 , ly
                Ys(i) = Yy(i)
              ENDDO
              incys = incy
              !
              !                    Call the subroutine.
              !
              CALL DGER(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
              !
              !                    Check if error-exit was taken incorrectly.
              !
              IF ( NUMXER(nerr)/=0 ) THEN
                IF ( Kprint>=2 ) WRITE (Nout,FMT=99007)
                Fatal = .TRUE.
              ENDIF
              !
              !                    See what data changed inside subroutine.
              !
              isame(1) = ms==m
              isame(2) = ns==n
              isame(3) = als==alpha
              isame(4) = LDE(Xs,Xx,lx)
              isame(5) = incxs==incx
              isame(6) = LDE(Ys,Yy,ly)
              isame(7) = incys==incy
              IF ( null ) THEN
                isame(8) = LDE(As,Aa,laa)
              ELSE
                isame(8) = LDERES('GE',' ',m,n,As,Aa,lda)
              ENDIF
              isame(9) = ldas==lda
              !
              !                    If data was incorrectly changed, report and return.
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
                !                       Check the result column by column.
                !
                IF ( incx>0 ) THEN
                  DO i = 1 , m
                    Z(i) = X(i)
                  ENDDO
                ELSE
                  DO i = 1 , m
                    Z(i) = X(m-i+1)
                  ENDDO
                ENDIF
                DO j = 1 , n
                  IF ( incy>0 ) THEN
                    w(1) = Y(j)
                  ELSE
                    w(1) = Y(n-j+1)
                  ENDIF
                  CALL DMVCH('N',m,1,alpha,Z,Nmax,w,1,ONE,A(1,j),1,Yt,G,&
                    Aa(1+(j-1)*lda),Eps,err,ftl,Nout,.TRUE.,Kprint)
                  errmax = MAX(errmax,err)
                ENDDO
              ENDIF
              IF ( ftl ) THEN
                Fatal = .TRUE.
                IF ( Kprint>=3 ) THEN
                  WRITE (Nout,FMT=99005) j
                  WRITE (Nout,FMT=99004) Sname
                  WRITE (Nout,FMT=99006) nc , Sname , m , n , alpha , incx ,&
                    incy , lda
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
  99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',&
    'ANGED INCORRECTLY *******')
  99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C','ALLS)',&
    /' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,' - SUSPECT *******')
  99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
  99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
  99006 FORMAT (1X,I6,': ',A6,'(',2(I3,','),F4.1,', X,',I2,', Y,',I2,', A,',I3,&
    ')                  .')
  99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
    '******')
  !
  !     End of DCHK42.
  !
END SUBROUTINE DCHK42
