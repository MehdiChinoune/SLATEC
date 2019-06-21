!** U11LS
SUBROUTINE U11LS(A,Mda,M,N,Ub,Db,Mode,Np,Krank,Ksure,H,W,Eb,Ic,Ir)
  !> Subsidiary to LLSIA
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (U11LS-S, DU11LS-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !       This routine performs a QR factorization of A
  !       using Householder transformations. Row and
  !       column pivots are chosen to reduce the growth
  !       of round-off and to help detect possible rank
  !       deficiency.
  !
  !***
  ! **See also:**  LLSIA
  !***
  ! **Routines called:**  ISAMAX, ISWAP, SAXPY, SDOT, SNRM2, SSCAL, SSWAP,
  !                    XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : XERMSG
  USE blas, ONLY : SAXPY, SSWAP
  INTEGER :: Mda, Mode, N, Np, Krank, Ksure, M
  INTEGER :: Ic(N), Ir(M)
  REAL(SP) :: A(Mda,N), Db(N), Eb(N), H(N), Ub(N), W(N)
  INTEGER :: mm, nmk, i, ii, im1, imin, is, j, jm1, jmax, jp1, kk, km1, kmi, &
    kp1, kz, l, lm1
  REAL(SP) :: bb, r2, rmin, summ, t, temp, tn, tt
  !
  !        INITIALIZATION
  !
  !* FIRST EXECUTABLE STATEMENT  U11LS
  j = 0
  Krank = N
  DO i = 1, N
    Ic(i) = i
  END DO
  DO i = 1, M
    Ir(i) = i
  END DO
  !
  !        DETERMINE REL AND ABS ERROR VECTORS
  !
  !
  !
  !        CALCULATE COL LENGTH
  !
  DO i = 1, N
    H(i) = NORM2(A(1:M,i))
    W(i) = H(i)
  END DO
  !
  !         INITIALIZE ERROR BOUNDS
  !
  DO i = 1, N
    Eb(i) = MAX(Db(i),Ub(i)*H(i))
    Ub(i) = Eb(i)
    Db(i) = 0._SP
  END DO
  !
  !          DISCARD SELF DEPENDENT COLUMNS
  !
  i = 1
  DO
    IF( Eb(i)>=H(i) ) THEN
      !
      !          MATRIX REDUCTION
      !
      kk = Krank
      Krank = Krank - 1
      IF( Mode==0 ) RETURN
      IF( i>Np ) THEN
        IF( i>Krank ) EXIT
        CALL SSWAP(1,Eb(i),1,Eb(kk),1)
        CALL SSWAP(1,Ub(i),1,Ub(kk),1)
        CALL SSWAP(1,W(i),1,W(kk),1)
        CALL SSWAP(1,H(i),1,H(kk),1)
        CALL ISWAP(1,Ic(i),1,Ic(kk),1)
        CALL SSWAP(M,A(1,i),1,A(1,kk),1)
      ELSE
        CALL XERMSG('U11LS','FIRST NP COLUMNS ARE LINEARLY DEPENDENT',8,0)
        Krank = i - 1
        RETURN
      END IF
    ELSE
      IF( i==Krank ) EXIT
      i = i + 1
    END IF
  END DO
  !
  !           TEST FOR ZERO RANK
  !
  IF( Krank<=0 ) THEN
    Krank = 0
    Ksure = 0
    RETURN
  END IF
  !
  !        M A I N    L O O P
  !
  100  j = j + 1
  jp1 = j + 1
  jm1 = j - 1
  kz = Krank
  IF( j<=Np ) kz = j
  !
  !        EACH COL HAS MM=M-J+1 COMPONENTS
  !
  mm = M - j + 1
  !
  !         UB DETERMINES COLUMN PIVOT
  !
  200  imin = j
  IF( H(j)/=0. ) THEN
    rmin = Ub(j)/H(j)
    DO i = j, kz
      IF( Ub(i)<H(i)*rmin ) THEN
        rmin = Ub(i)/H(i)
        imin = i
      END IF
    END DO
    !
    !     TEST FOR RANK DEFICIENCY
    !
    IF( rmin<1._SP ) GOTO 400
    tt = (Eb(imin)+ABS(Db(imin)))/H(imin)
    IF( tt<1._SP ) THEN
      !     COMPUTE EXACT UB
      DO i = 1, jm1
        W(i) = A(i,imin)
      END DO
      l = jm1
      DO
        W(l) = W(l)/A(l,l)
        IF( l==1 ) THEN
          tt = Eb(imin)
          DO i = 1, jm1
            tt = tt + ABS(W(i))*Eb(i)
          END DO
          Ub(imin) = tt
          IF( Ub(imin)/H(imin)<1._SP ) GOTO 400
          EXIT
        ELSE
          lm1 = l - 1
          DO i = l, jm1
            W(lm1) = W(lm1) - A(lm1,i)*W(i)
          END DO
          l = lm1
        END IF
      END DO
    END IF
  END IF
  !
  !        MATRIX REDUCTION
  !
  300  kk = Krank
  Krank = Krank - 1
  kz = Krank
  IF( Mode==0 ) RETURN
  IF( j>Np ) THEN
    IF( imin<=Krank ) THEN
      CALL ISWAP(1,Ic(imin),1,Ic(kk),1)
      CALL SSWAP(M,A(1,imin),1,A(1,kk),1)
      CALL SSWAP(1,Eb(imin),1,Eb(kk),1)
      CALL SSWAP(1,Ub(imin),1,Ub(kk),1)
      CALL SSWAP(1,Db(imin),1,Db(kk),1)
      CALL SSWAP(1,W(imin),1,W(kk),1)
      CALL SSWAP(1,H(imin),1,H(kk),1)
    END IF
    IF( j<=Krank ) GOTO 200
    GOTO 500
  ELSE
    CALL XERMSG('U11LS','FIRST NP COLUMNS ARE LINEARLY DEPENDENT',8,0)
    Krank = j - 1
    RETURN
  END IF
  !
  !        COLUMN PIVOT
  !
  400 CONTINUE
  IF( imin/=j ) THEN
    CALL SSWAP(1,H(j),1,H(imin),1)
    CALL SSWAP(M,A(1,j),1,A(1,imin),1)
    CALL SSWAP(1,Eb(j),1,Eb(imin),1)
    CALL SSWAP(1,Ub(j),1,Ub(imin),1)
    CALL SSWAP(1,Db(j),1,Db(imin),1)
    CALL SSWAP(1,W(j),1,W(imin),1)
    CALL ISWAP(1,Ic(j),1,Ic(imin),1)
  END IF
  !
  !        ROW PIVOT
  !
  jmax = MAXLOC(A(j:M,j),1)
  jmax = jmax + j - 1
  IF( jmax/=j ) THEN
    CALL SSWAP(N,A(j,1),Mda,A(jmax,1),Mda)
    CALL ISWAP(1,Ir(j),1,Ir(jmax),1)
  END IF
  !
  !     APPLY HOUSEHOLDER TRANSFORMATION
  !
  tn = NORM2(A(j:M,j))
  IF( tn==0._SP ) GOTO 300
  IF( A(j,j)/=0._SP ) tn = SIGN(tn,A(j,j))
  A(j:M,j) = A(j:M,j)/tn
  A(j,j) = A(j,j) + 1._SP
  IF( j/=N ) THEN
    DO i = jp1, N
      bb = -DOT_PRODUCT(A(j:M,j),A(j:M,i))/A(j,j)
      CALL SAXPY(mm,bb,A(j:M,j),1,A(j:M,i),1)
      IF( i>Np ) THEN
        IF( H(i)/=0._SP ) THEN
          tt = 1._SP - (ABS(A(j,i))/H(i))**2
          tt = MAX(tt,0._SP)
          t = tt
          tt = 1._SP + 0.05_SP*tt*(H(i)/W(i))**2
          IF( tt==1._SP ) THEN
            H(i) = NORM2(A(j+1:M,i))
            W(i) = H(i)
          ELSE
            H(i) = H(i)*SQRT(t)
          END IF
        END IF
      END IF
    END DO
  END IF
  H(j) = A(j,j)
  A(j,j) = -tn
  !
  !
  !          UPDATE UB, DB
  !
  Ub(j) = Ub(j)/ABS(A(j,j))
  Db(j) = (SIGN(Eb(j),Db(j))+Db(j))/A(j,j)
  IF( j/=Krank ) THEN
    DO i = jp1, Krank
      Ub(i) = Ub(i) + ABS(A(j,i))*Ub(j)
      Db(i) = Db(i) - A(j,i)*Db(j)
    END DO
    GOTO 100
  END IF
  !
  !        E N D    M A I N    L O O P
  !
  !
  !        COMPUTE KSURE
  !
  500  km1 = Krank - 1
  DO i = 1, km1
    is = 0
    kmi = Krank - i
    DO ii = 1, kmi
      IF( Ub(ii)>Ub(ii+1) ) THEN
        is = 1
        temp = Ub(ii)
        Ub(ii) = Ub(ii+1)
        Ub(ii+1) = temp
      END IF
    END DO
    IF( is==0 ) EXIT
  END DO
  Ksure = 0
  summ = 0._SP
  DO i = 1, Krank
    r2 = Ub(i)*Ub(i)
    IF( r2+summ>=1._SP ) EXIT
    summ = summ + r2
    Ksure = Ksure + 1
  END DO
  !
  !     IF SYSTEM IS OF REDUCED RANK AND MODE = 2
  !     COMPLETE THE DECOMPOSITION FOR SHORTEST LEAST SQUARES SOLUTION
  !
  IF( Krank/=N .AND. Mode>=2 ) THEN
    nmk = N - Krank
    kp1 = Krank + 1
    i = Krank
    DO
      tn = NORM2(A(i,kp1:N))/A(i,i)
      tn = A(i,i)*SQRT(1._SP+tn*tn)
      A(i,kp1:N) = A(i,kp1:N)/tn
      W(i) = A(i,i)/tn + 1._SP
      A(i,i) = -tn
      IF( i==1 ) EXIT
      im1 = i - 1
      DO ii = 1, im1
        tt = -DOT_PRODUCT(A(ii,kp1:N),A(i,kp1:N))/W(i)
        tt = tt - A(ii,i)
        CALL SAXPY(nmk,tt,A(i,kp1:N),1,A(ii,kp1:N),1)
        A(ii,i) = A(ii,i) + tt*W(i)
      END DO
      i = i - 1
    END DO
  END IF
END SUBROUTINE U11LS
