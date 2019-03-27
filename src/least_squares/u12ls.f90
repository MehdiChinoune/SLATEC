!** U12LS
SUBROUTINE U12LS(A,Mda,M,N,B,Mdb,Nb,Mode,Krank,Rnorm,H,W,Ic,Ir)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to LLSIA
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (U12LS-S, DU12LS-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !        Given the Householder QR factorization of A, this
  !        subroutine solves the system AX=B. If the system
  !        is of reduced rank, this routine returns a solution
  !        according to the selected mode.
  !
  !       Note - If MODE.NE.2, W is never accessed.
  !
  !***
  ! **See also:**  LLSIA
  !***
  ! **Routines called:**  SAXPY, SDOT, SNRM2, SSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  
  INTEGER i, ij, im1, j, jb, k, kp1, Krank, M, Mda, Mdb, Mode, N, Nb, nmk
  REAL A(Mda,*), B(Mdb,*), bb, H(*), Rnorm(*), SDOT, SNRM2, tt, W(*)
  INTEGER Ic(*), Ir(*)
  !* FIRST EXECUTABLE STATEMENT  U12LS
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
        ENDDO
        !
        !     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
        !
        DO j = 1, k
          tt = A(j,j)
          A(j,j) = H(j)
          DO i = 1, Nb
            bb = -SDOT(M-j+1,A(j,j),1,B(j,i),1)/H(j)
            CALL SAXPY(M-j+1,bb,A(j,j),1,B(j,i),1)
          ENDDO
          A(j,j) = tt
        ENDDO
        !
        !        FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B)
        !
        DO jb = 1, Nb
          Rnorm(jb) = SNRM2((M-k),B(kp1,jb),1)
        ENDDO
        !
        !     BACK SOLVE UPPER TRIANGULAR R
        !
        i = k
        DO
          DO jb = 1, Nb
            B(i,jb) = B(i,jb)/A(i,i)
          ENDDO
          IF ( i==1 ) THEN
            !
            !     RANK LT N
            !
            !      TRUNCATED SOLUTION
            !
            IF ( k/=N ) THEN
              DO jb = 1, Nb
                DO i = kp1, N
                  B(i,jb) = 0.0
                ENDDO
              ENDDO
              IF ( Mode/=1 ) THEN
                !
                !      MINIMAL LENGTH SOLUTION
                !
                nmk = N - k
                DO jb = 1, Nb
                  DO i = 1, k
                    tt = -SDOT(nmk,A(i,kp1),Mda,B(kp1,jb),1)/W(i)
                    tt = tt - B(i,jb)
                    CALL SAXPY(nmk,tt,A(i,kp1),Mda,B(kp1,jb),1)
                    B(i,jb) = B(i,jb) + tt*W(i)
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
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
                ENDDO
                RETURN
              ELSE
                j = Ic(i)
                IF ( j/=i ) THEN
                  IF ( j>=0 ) THEN
                    Ic(i) = -Ic(i)
                    DO
                      CALL SSWAP(Nb,B(j,1),Mdb,B(i,1),Mdb)
                      ij = Ic(j)
                      Ic(j) = -Ic(j)
                      j = ij
                      IF ( j==i ) EXIT
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE
            im1 = i - 1
            DO jb = 1, Nb
              CALL SAXPY(im1,-B(i,jb),A(1,i),1,B(1,jb),1)
            ENDDO
            i = im1
          ENDIF
        ENDDO
      ELSE
        j = Ir(i)
        IF ( j/=i ) THEN
          IF ( j>=0 ) THEN
            Ir(i) = -Ir(i)
            DO jb = 1, Nb
              Rnorm(jb) = B(i,jb)
            ENDDO
            ij = i
            DO
              DO jb = 1, Nb
                B(ij,jb) = B(j,jb)
              ENDDO
              ij = j
              j = Ir(ij)
              Ir(ij) = -Ir(ij)
              IF ( j==i ) THEN
                DO jb = 1, Nb
                  B(ij,jb) = Rnorm(jb)
                ENDDO
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ELSE
    DO jb = 1, Nb
      Rnorm(jb) = SNRM2(M,B(1,jb),1)
    ENDDO
    DO jb = 1, Nb
      DO i = 1, N
        B(i,jb) = 0.0
      ENDDO
    ENDDO
    RETURN
  ENDIF
  !
  !        SOLUTION VECTORS ARE IN FIRST N ROWS OF B(,)
  !
  RETURN
END SUBROUTINE U12LS
