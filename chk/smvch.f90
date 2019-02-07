!*==SMVCH.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SMVCH
SUBROUTINE SMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,Yy,Eps,&
    Err,Ftl,Nout,Mv,Kprint)
  IMPLICIT NONE
  !*--SMVCH6
  !***BEGIN PROLOGUE  SMVCH
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
  !***END PROLOGUE  SMVCH
  !     .. Parameters ..
  REAL ZERO, ONE
  PARAMETER (ZERO=0.0,ONE=1.0)
  !     .. Scalar Arguments ..
  REAL Alpha, Beta, Eps, Err
  INTEGER Incx, Incy, Kprint, M, N, Nmax, Nout
  LOGICAL Mv, Ftl
  CHARACTER :: Trans
  !     .. Array Arguments ..
  REAL A(Nmax,*), G(*), X(*), Y(*), Yt(*), Yy(*)
  !     .. Local Scalars ..
  REAL erri
  INTEGER i, incxl, incyl, iy, j, jx, k, kx, ky, ml, nl
  LOGICAL tran
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, SQRT
  !***FIRST EXECUTABLE STATEMENT  SMVCH
  tran = Trans=='T' .OR. Trans=='C'
  IF ( tran ) THEN
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
    G(iy) = ZERO
    jx = kx
    IF ( tran ) THEN
      DO j = 1, nl
        Yt(iy) = Yt(iy) + A(j,i)*X(jx)
        G(iy) = G(iy) + ABS(A(j,i)*X(jx))
        jx = jx + incxl
      ENDDO
    ELSE
      DO j = 1, nl
        Yt(iy) = Yt(iy) + A(i,j)*X(jx)
        G(iy) = G(iy) + ABS(A(i,j)*X(jx))
        jx = jx + incxl
      ENDDO
    ENDIF
    Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
    G(iy) = ABS(Alpha)*G(iy) + ABS(Beta*Y(iy))
    iy = iy + incyl
  ENDDO
  !
  !     Compute the error ratio for this result.
  !
  Err = ZERO
  DO i = 1, ml
    erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
    IF ( G(i)/=ZERO ) erri = erri/G(i)
    Err = MAX(Err,erri)
    IF ( Err*SQRT(Eps)>=ONE ) THEN
      Ftl = .TRUE.
      IF ( Kprint>=2 ) THEN
        WRITE (Nout,FMT=99001)
        DO k = 1, ml
          IF ( Mv ) THEN
            WRITE (Nout,FMT=99002) k, Yt(k), Yy(1+(k-1)*ABS(Incy))
          ELSE
            WRITE (Nout,FMT=99002) k, Yy(1+(k-1)*ABS(Incy)), Yt(k)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  RETURN
  !
  99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
    'F ACCURATE *******',/'           EXPECTED RESULT   COMPU',&
    'TED RESULT')
  99002 FORMAT (1X,I7,2G18.6)
  !
  !     End of SMVCH.
  !
END SUBROUTINE SMVCH
