!** DSPSL
SUBROUTINE DSPSL(Ap,N,Kpvt,B)
  !>
  !***
  !  Solve a real symmetric system using the factors obtained
  !            from DSPFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1A
  !***
  ! **Type:**      DOUBLE PRECISION (SSPSL-S, DSPSL-D, CHPSL-C, CSPSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, SOLVE, SYMMETRIC
  !***
  ! **Author:**  Bunch, J., (UCSD)
  !***
  ! **Description:**
  !
  !     DSISL solves the double precision symmetric system
  !     A * X = B
  !     using the factors computed by DSPFA.
  !
  !     On Entry
  !
  !        AP      DOUBLE PRECISION(N*(N+1)/2)
  !                the output from DSPFA.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        KPVT    INTEGER(N)
  !                the pivot vector from DSPFA.
  !
  !        B       DOUBLE PRECISION(N)
  !                the right hand side vector.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero may occur if  DSPCO  has set RCOND .EQ. 0.0
  !        or  DSPFA  has set INFO .NE. 0  .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DSPFA(AP,N,KPVT,INFO)
  !           IF (INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DSPSL(AP,N,KPVT,C(1,J))
  !        10 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Modified routine equivalence list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N, Kpvt(*)
  REAL(8) :: Ap(*), B(*)
  !
  REAL(8) :: ak, akm1, bk, bkm1, denom, temp
  INTEGER ik, ikm1, ikp1, k, kk, km1k, km1km1, kp
  !
  !     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
  !     D INVERSE TO B.
  !
  !* FIRST EXECUTABLE STATEMENT  DSPSL
  k = N
  ik = (N*(N-1))/2
  DO WHILE ( k/=0 )
    kk = ik + k
    IF ( Kpvt(k)<0 ) THEN
      !
      !           2 X 2 PIVOT BLOCK.
      !
      ikm1 = ik - (k-1)
      IF ( k/=2 ) THEN
        kp = ABS(Kpvt(k))
        IF ( kp/=k-1 ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k-1)
          B(k-1) = B(kp)
          B(kp) = temp
        END IF
        !
        !              APPLY THE TRANSFORMATION.
        !
        CALL DAXPY(k-2,B(k),Ap(ik+1),1,B(1),1)
        CALL DAXPY(k-2,B(k-1),Ap(ikm1+1),1,B(1),1)
      END IF
      !
      !           APPLY D INVERSE.
      !
      km1k = ik + k - 1
      kk = ik + k
      ak = Ap(kk)/Ap(km1k)
      km1km1 = ikm1 + k - 1
      akm1 = Ap(km1km1)/Ap(km1k)
      bk = B(k)/Ap(km1k)
      bkm1 = B(k-1)/Ap(km1k)
      denom = ak*akm1 - 1.0D0
      B(k) = (akm1*bk-bkm1)/denom
      B(k-1) = (ak*bkm1-bk)/denom
      k = k - 2
      ik = ik - (k+1) - k
    ELSE
      !
      !           1 X 1 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        kp = Kpvt(k)
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        END IF
        !
        !              APPLY THE TRANSFORMATION.
        !
        CALL DAXPY(k-1,B(k),Ap(ik+1),1,B(1),1)
      END IF
      !
      !           APPLY D INVERSE.
      !
      B(k) = B(k)/Ap(kk)
      k = k - 1
      ik = ik - k
    END IF
  END DO
  !
  !     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
  !
  k = 1
  ik = 0
  DO WHILE ( k<=N )
    IF ( Kpvt(k)<0 ) THEN
      !
      !           2 X 2 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        !
        !              APPLY THE TRANSFORMATION.
        !
        B(k) = B(k) + DDOT(k-1,Ap(ik+1),1,B(1),1)
        ikp1 = ik + k
        B(k+1) = B(k+1) + DDOT(k-1,Ap(ikp1+1),1,B(1),1)
        kp = ABS(Kpvt(k))
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        END IF
      END IF
      ik = ik + k + k + 1
      k = k + 2
    ELSE
      !
      !           1 X 1 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        !
        !              APPLY THE TRANSFORMATION.
        !
        B(k) = B(k) + DDOT(k-1,Ap(ik+1),1,B(1),1)
        kp = Kpvt(k)
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        END IF
      END IF
      ik = ik + k
      k = k + 1
    END IF
  END DO
END SUBROUTINE DSPSL
