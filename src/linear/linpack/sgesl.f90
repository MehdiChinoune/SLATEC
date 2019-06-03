!** SGESL
SUBROUTINE SGESL(A,Lda,N,Ipvt,B,Job)
  !>
  !  Solve the real system A*X=B or TRANS(A)*X=B using the
  !            factors of SGECO or SGEFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A1
  !***
  ! **Type:**      SINGLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     SGESL solves the real system
  !     A * X = B  or  TRANS(A) * X = B
  !     using the factors computed by SGECO or SGEFA.
  !
  !     On Entry
  !
  !        A       REAL(LDA, N)
  !                the output from SGECO or SGEFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from SGECO or SGEFA.
  !
  !        B       REAL(N)
  !                the right hand side vector.
  !
  !        JOB     INTEGER
  !                = 0         to solve  A*X = B ,
  !                = nonzero   to solve  TRANS(A)*X = B  where
  !                            TRANS(A)  is the transpose.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains a
  !        zero on the diagonal.  Technically, this indicates singularity,
  !        but it is often caused by improper arguments or improper
  !        setting of LDA .  It will not occur if the subroutines are
  !        called correctly and if SGECO has set RCOND .GT. 0.0
  !        or SGEFA has set INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
  !           IF (RCOND is too small) GO TO ...
  !           DO 10 J = 1, P
  !              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
  !        10 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  SAXPY, SDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: Lda, N, Job, Ipvt(N)
  REAL(SP) :: A(Lda,N), B(N)
  !
  INTEGER :: k, kb, l, nm1
  REAL(SP) :: t
  !* FIRST EXECUTABLE STATEMENT  SGESL
  nm1 = N - 1
  IF ( Job/=0 ) THEN
    !
    !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
    !        FIRST SOLVE  TRANS(U)*Y = B
    !
    DO k = 1, N
      t = SDOT(k-1,A(1,k),1,B(1),1)
      B(k) = (B(k)-t)/A(k,k)
    END DO
    !
    !        NOW SOLVE TRANS(L)*X = Y
    !
    IF ( nm1>=1 ) THEN
      DO kb = 1, nm1
        k = N - kb
        B(k) = B(k) + SDOT(N-k,A(k+1,k),1,B(k+1),1)
        l = Ipvt(k)
        IF ( l/=k ) THEN
          t = B(l)
          B(l) = B(k)
          B(k) = t
        END IF
      END DO
    END IF
  ELSE
    !
    !        JOB = 0, SOLVE  A * X = B
    !        FIRST SOLVE  L*Y = B
    !
    IF ( nm1>=1 ) THEN
      DO k = 1, nm1
        l = Ipvt(k)
        t = B(l)
        IF ( l/=k ) THEN
          B(l) = B(k)
          B(k) = t
        END IF
        CALL SAXPY(N-k,t,A(k+1,k),1,B(k+1),1)
      END DO
    END IF
    !
    !        NOW SOLVE  U*X = Y
    !
    DO kb = 1, N
      k = N + 1 - kb
      B(k) = B(k)/A(k,k)
      t = -B(k)
      CALL SAXPY(k-1,t,A(1,k),1,B(1),1)
    END DO
  END IF
END SUBROUTINE SGESL
