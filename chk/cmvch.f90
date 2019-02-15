!DECK CMVCH
SUBROUTINE CMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,Yy,Eps,&
    Err,Ftl,Nout,Mv,Kprint)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CMVCH
  !***SUBSIDIARY
  !***PURPOSE  Check the results of the computational tests.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. J., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  Checks the results of the computational tests.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CMVCH
  !     .. Parameters ..
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0,0.0))
  REAL RZERO, RONE
  PARAMETER (RZERO=0.0,RONE=1.0)
  !     .. Scalar Arguments ..
  COMPLEX Alpha, Beta
  REAL Eps, Err
  INTEGER Incx, Incy, Kprint, M, N, Nmax, Nout
  LOGICAL Mv, Ftl
  CHARACTER :: Trans
  !     .. Array Arguments ..
  COMPLEX A(Nmax,*), X(*), Y(*), Yt(*), Yy(*)
  REAL G(*)
  !     .. Local Scalars ..
  COMPLEX c
  REAL erri
  INTEGER i, incxl, incyl, iy, j, jx, k, kx, ky, ml, nl
  LOGICAL ctran, tran
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, AIMAG, CONJG, MAX, REAL, SQRT
  REAL CABS1
  !***FIRST EXECUTABLE STATEMENT  CMVCH
  tran = Trans=='T'
  ctran = Trans=='C'
  IF ( tran.OR.ctran ) THEN
    ml = N
    nl = M
  ELSE
    ml = M
    nl = N
  ENDIF
  IF ( Incx<0 ) THEN
    kx = nl
    incxl = -1
  ELSE
    kx = 1
    incxl = 1
  ENDIF
  IF ( Incy<0 ) THEN
    ky = ml
    incyl = -1
  ELSE
    ky = 1
    incyl = 1
  ENDIF
  !
  !     Compute expected result in YT using data in A, X and Y.
  !     Compute gauges in G.
  !
  iy = ky
  DO i = 1, ml
    Yt(iy) = ZERO
    G(iy) = RZERO
    jx = kx
    IF ( tran ) THEN
      DO j = 1, nl
        Yt(iy) = Yt(iy) + A(j,i)*X(jx)
        G(iy) = G(iy) + CABS1(A(j,i))*CABS1(X(jx))
        jx = jx + incxl
      ENDDO
    ELSEIF ( ctran ) THEN
      DO j = 1, nl
        Yt(iy) = Yt(iy) + CONJG(A(j,i))*X(jx)
        G(iy) = G(iy) + CABS1(A(j,i))*CABS1(X(jx))
        jx = jx + incxl
      ENDDO
    ELSE
      DO j = 1, nl
        Yt(iy) = Yt(iy) + A(i,j)*X(jx)
        G(iy) = G(iy) + CABS1(A(i,j))*CABS1(X(jx))
        jx = jx + incxl
      ENDDO
    ENDIF
    Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
    G(iy) = CABS1(Alpha)*G(iy) + CABS1(Beta)*CABS1(Y(iy))
    iy = iy + incyl
  ENDDO
  !
  !     Compute the error ratio for this result.
  !
  Err = ZERO
  DO i = 1, ml
    erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
    IF ( G(i)/=RZERO ) erri = erri/G(i)
    Err = MAX(Err,erri)
    IF ( Err*SQRT(Eps)>=RONE ) THEN
      Ftl = .TRUE.
      IF ( Kprint>=2 ) THEN
        WRITE (Nout,FMT=99001)
        DO k = 1, ml
          IF ( Mv ) THEN
            WRITE (Nout,FMT=99002) k, Yt(k), Yy(1+(k-1)*ABS(Incy))
          ELSE
            WRITE (Nout,FMT=99002) i, Yy(1+(k-1)*ABS(Incy)), Yt(k)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
  ENDDO
  RETURN
  !
  99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
    'F ACCURATE *******',/'                       EXPECTED RE',&
    'SULT                    COMPUTED RESULT')
  99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
  !
  !     End of CMVCH.
  !
END SUBROUTINE CMVCH
