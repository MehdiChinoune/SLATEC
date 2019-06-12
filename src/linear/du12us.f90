!** DU12US
SUBROUTINE DU12US(A,Mda,M,N,B,Mdb,Nb,Mode,Krank,Rnorm,H,W,Ir,Ic)
  !>
  !  Subsidiary to DULSIA
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (U12US-S, DU12US-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !        Given the Householder LQ factorization of A, this
  !        subroutine solves the system AX=B. If the system
  !        is of reduced rank, this routine returns a solution
  !        according to the selected mode.
  !
  !       Note - If MODE.NE.2, W is never accessed.
  !
  !***
  ! **See also:**  DULSIA
  !***
  ! **Routines called:**  DAXPY, DDOT, DNRM2, DSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  USE blas, ONLY : DAXPY, DSWAP
  INTEGER :: Krank, M, Mda, Mdb, Mode, N, Nb
  INTEGER :: Ic(N), Ir(M)
  REAL(DP) :: A(Mda,N), B(Mdb,Nb), H(M), Rnorm(Nb), W(4*M)
  INTEGER :: i, ij, ip1, j, jb, k, kp1, mmk
  REAL(DP) :: bb, tt
  !* FIRST EXECUTABLE STATEMENT  DU12US
  k = Krank
  kp1 = k + 1
  !
  !        RANK=0
  !
  IF ( k>0 ) THEN
    !
    !     REORDER B TO REFLECT ROW INTERCHANGES
    !
    i = 0
    DO
      i = i + 1
      IF ( i==M ) THEN
        DO i = 1, M
          Ir(i) = ABS(Ir(i))
        END DO
        !
        !     IF A IS OF REDUCED RANK AND MODE=2,
        !     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
        !
        IF ( Mode>=2.AND.k/=M ) THEN
          mmk = M - k
          DO jb = 1, Nb
            DO j = 1, k
              i = kp1 - j
              tt = -DOT_PRODUCT(A(kp1:M,i),B(kp1:M,jb))/W(i)
              tt = tt - B(i,jb)
              CALL DAXPY(mmk,tt,A(kp1,i),1,B(kp1,jb),1)
              B(i,jb) = B(i,jb) + tt*W(i)
            END DO
          END DO
        END IF
        !
        !     FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B)
        !
        DO jb = 1, Nb
          Rnorm(jb) = NORM2(B(kp1:M,jb))
        END DO
        !
        !     BACK SOLVE LOWER TRIANGULAR L
        !
        DO jb = 1, Nb
          DO i = 1, k
            B(i,jb) = B(i,jb)/A(i,i)
            IF ( i==k ) EXIT
            ip1 = i + 1
            CALL DAXPY(k-i,-B(i,jb),A(ip1,i),1,B(ip1,jb),1)
          END DO
        END DO
        !
        !
        !      TRUNCATED SOLUTION
        !
        IF ( k/=N ) THEN
          DO jb = 1, Nb
            DO i = kp1, N
              B(i,jb) = 0.0D0
            END DO
          END DO
        END IF
        !
        !     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
        !
        DO i = 1, k
          j = kp1 - i
          tt = A(j,j)
          A(j,j) = H(j)
          DO jb = 1, Nb
            bb = -DOT_PRODUCT(A(j,j:N),B(j:N,jb))/H(j)
            CALL DAXPY(N-j+1,bb,A(j,j),Mda,B(j,jb),1)
          END DO
          A(j,j) = tt
        END DO
        !
        !
        !     REORDER B TO REFLECT COLUMN INTERCHANGES
        !
        i = 0
        DO
          i = i + 1
          IF ( i==N ) THEN
            DO i = 1, N
              Ic(i) = ABS(Ic(i))
            END DO
            RETURN
          ELSE
            j = Ic(i)
            IF ( j/=i ) THEN
              IF ( j>=0 ) THEN
                Ic(i) = -Ic(i)
                DO
                  CALL DSWAP(Nb,B(j,1),Mdb,B(i,1),Mdb)
                  ij = Ic(j)
                  Ic(j) = -Ic(j)
                  j = ij
                  IF ( j==i ) EXIT
                END DO
              END IF
            END IF
          END IF
        END DO
      ELSE
        j = Ir(i)
        IF ( j/=i ) THEN
          IF ( j>=0 ) THEN
            Ir(i) = -Ir(i)
            DO jb = 1, Nb
              Rnorm(jb) = B(i,jb)
            END DO
            ij = i
            DO
              DO jb = 1, Nb
                B(ij,jb) = B(j,jb)
              END DO
              ij = j
              j = Ir(ij)
              Ir(ij) = -Ir(ij)
              IF ( j==i ) THEN
                DO jb = 1, Nb
                  B(ij,jb) = Rnorm(jb)
                END DO
                EXIT
              END IF
            END DO
          END IF
        END IF
      END IF
    END DO
  ELSE
    DO jb = 1, Nb
      Rnorm(jb) = NORM2(B(1:M,jb))
    END DO
    DO jb = 1, Nb
      DO i = 1, N
        B(i,jb) = 0.0D0
      END DO
    END DO
    RETURN
  END IF
  !
  !        SOLUTION VECTORS ARE IN FIRST N ROWS OF B(,)
  !
  RETURN
END SUBROUTINE DU12US
