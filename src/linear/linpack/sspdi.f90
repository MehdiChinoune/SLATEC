!** SSPDI
SUBROUTINE SSPDI(Ap,N,Kpvt,Det,Inert,Work,Job)
  IMPLICIT NONE
  !>
  !***
  !  Compute the determinant, inertia, inverse of a real
  !            symmetric matrix stored in packed form using the factors
  !            from SSPFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1A, D3B1A
  !***
  ! **Type:**      SINGLE PRECISION (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C)
  !***
  ! **Keywords:**  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             PACKED, SYMMETRIC
  !***
  ! **Author:**  Bunch, J., (UCSD)
  !***
  ! **Description:**
  !
  !     SSPDI computes the determinant, inertia and inverse
  !     of a real symmetric matrix using the factors from SSPFA,
  !     where the matrix is stored in packed form.
  !
  !     On Entry
  !
  !        AP      REAL (N*(N+1)/2)
  !                the output from SSPFA.
  !
  !        N       INTEGER
  !                the order of the matrix A.
  !
  !        KPVT    INTEGER(N)
  !                the pivot vector from SSPFA.
  !
  !        WORK    REAL(N)
  !                work vector.  Contents ignored.
  !
  !        JOB     INTEGER
  !                JOB has the decimal expansion  ABC  where
  !                   If  C .NE. 0, the inverse is computed,
  !                   If  B .NE. 0, the determinant is computed,
  !                   If  A .NE. 0, the inertia is computed.
  !
  !                For example, JOB = 111  gives all three.
  !
  !     On Return
  !
  !        Variables not requested by JOB are not used.
  !
  !        AP     contains the upper triangle of the inverse of
  !               the original matrix, stored in packed form.
  !               The columns of the upper triangle are stored
  !               sequentially in a one-dimensional array.
  !
  !        DET    REAL(2)
  !               determinant of original matrix.
  !               Determinant = DET(1) * 10.0**DET(2)
  !               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
  !               or DET(1) = 0.0.
  !
  !        INERT  INTEGER(3)
  !               the inertia of the original matrix.
  !               INERT(1)  =  number of positive eigenvalues.
  !               INERT(2)  =  number of negative eigenvalues.
  !               INERT(3)  =  number of zero eigenvalues.
  !
  !     Error Condition
  !
  !        A division by zero will occur if the inverse is requested
  !        and  SSPCO  has set RCOND .EQ. 0.0
  !        or  SSPFA  has set  INFO .NE. 0 .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  SAXPY, SCOPY, SDOT, SSWAP

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

  INTEGER N, Job
  REAL Ap(*), Work(*)
  REAL Det(2)
  INTEGER Kpvt(*), Inert(3)
  !
  REAL akkp1, SDOT, temp
  REAL ten, d, t, ak, akp1
  INTEGER ij, ik, ikp1, iks, j, jb, jk, jkp1
  INTEGER k, kk, kkp1, km1, ks, ksj, kskp1, kstep
  LOGICAL noinv, nodet, noert
  !* FIRST EXECUTABLE STATEMENT  SSPDI
  noinv = MOD(Job,10)==0
  nodet = MOD(Job,100)/10==0
  noert = MOD(Job,1000)/100==0
  !
  IF ( .NOT.(nodet.AND.noert) ) THEN
    IF ( .NOT.(noert) ) THEN
      Inert(1) = 0
      Inert(2) = 0
      Inert(3) = 0
    END IF
    IF ( .NOT.(nodet) ) THEN
      Det(1) = 1.0E0
      Det(2) = 0.0E0
      ten = 10.0E0
    END IF
    t = 0.0E0
    ik = 0
    DO k = 1, N
      kk = ik + k
      d = Ap(kk)
      !
      !           CHECK IF 1 BY 1
      !
      IF ( Kpvt(k)<=0 ) THEN
        !
        !              2 BY 2 BLOCK
        !              USE DET (D  S)  =  (D/T * C - T) * T ,  T = ABS(S)
        !                      (S  C)
        !              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
        !              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
        !
        IF ( t/=0.0E0 ) THEN
          d = t
          t = 0.0E0
        ELSE
          ikp1 = ik + k
          kkp1 = ikp1 + k
          t = ABS(Ap(kkp1))
          d = (d/t)*Ap(kkp1+1) - t
        END IF
      END IF
      !
      IF ( .NOT.(noert) ) THEN
        IF ( d>0.0E0 ) Inert(1) = Inert(1) + 1
        IF ( d<0.0E0 ) Inert(2) = Inert(2) + 1
        IF ( d==0.0E0 ) Inert(3) = Inert(3) + 1
      END IF
      !
      IF ( .NOT.(nodet) ) THEN
        Det(1) = d*Det(1)
        IF ( Det(1)/=0.0E0 ) THEN
          DO WHILE ( ABS(Det(1))<1.0E0 )
            Det(1) = ten*Det(1)
            Det(2) = Det(2) - 1.0E0
          END DO
          DO WHILE ( ABS(Det(1))>=ten )
            Det(1) = Det(1)/ten
            Det(2) = Det(2) + 1.0E0
          END DO
        END IF
      END IF
      ik = ik + k
    END DO
  END IF
  !
  !     COMPUTE INVERSE(A)
  !
  IF ( .NOT.(noinv) ) THEN
    k = 1
    ik = 0
    DO WHILE ( k<=N )
      km1 = k - 1
      kk = ik + k
      ikp1 = ik + k
      kkp1 = ikp1 + k
      IF ( Kpvt(k)<0 ) THEN
        !
        !              2 BY 2
        !
        t = ABS(Ap(kkp1))
        ak = Ap(kk)/t
        akp1 = Ap(kkp1+1)/t
        akkp1 = Ap(kkp1)/t
        d = t*(ak*akp1-1.0E0)
        Ap(kk) = akp1/d
        Ap(kkp1+1) = ak/d
        Ap(kkp1) = -akkp1/d
        IF ( km1>=1 ) THEN
          CALL SCOPY(km1,Ap(ikp1+1),1,Work,1)
          ij = 0
          DO j = 1, km1
            jkp1 = ikp1 + j
            Ap(jkp1) = SDOT(j,Ap(ij+1),1,Work,1)
            CALL SAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ikp1+1),1)
            ij = ij + j
          END DO
          Ap(kkp1+1) = Ap(kkp1+1) + SDOT(km1,Work,1,Ap(ikp1+1),1)
          Ap(kkp1) = Ap(kkp1) + SDOT(km1,Ap(ik+1),1,Ap(ikp1+1),1)
          CALL SCOPY(km1,Ap(ik+1),1,Work,1)
          ij = 0
          DO j = 1, km1
            jk = ik + j
            Ap(jk) = SDOT(j,Ap(ij+1),1,Work,1)
            CALL SAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ik+1),1)
            ij = ij + j
          END DO
          Ap(kk) = Ap(kk) + SDOT(km1,Work,1,Ap(ik+1),1)
        END IF
        kstep = 2
      ELSE
        !
        !              1 BY 1
        !
        Ap(kk) = 1.0E0/Ap(kk)
        IF ( km1>=1 ) THEN
          CALL SCOPY(km1,Ap(ik+1),1,Work,1)
          ij = 0
          DO j = 1, km1
            jk = ik + j
            Ap(jk) = SDOT(j,Ap(ij+1),1,Work,1)
            CALL SAXPY(j-1,Work(j),Ap(ij+1),1,Ap(ik+1),1)
            ij = ij + j
          END DO
          Ap(kk) = Ap(kk) + SDOT(km1,Work,1,Ap(ik+1),1)
        END IF
        kstep = 1
      END IF
      !
      !           SWAP
      !
      ks = ABS(Kpvt(k))
      IF ( ks/=k ) THEN
        iks = (ks*(ks-1))/2
        CALL SSWAP(ks,Ap(iks+1),1,Ap(ik+1),1)
        ksj = ik + ks
        DO jb = ks, k
          j = k + ks - jb
          jk = ik + j
          temp = Ap(jk)
          Ap(jk) = Ap(ksj)
          Ap(ksj) = temp
          ksj = ksj - (j-1)
        END DO
        IF ( kstep/=1 ) THEN
          kskp1 = ikp1 + ks
          temp = Ap(kskp1)
          Ap(kskp1) = Ap(kkp1)
          Ap(kkp1) = temp
        END IF
      END IF
      ik = ik + k
      IF ( kstep==2 ) ik = ik + k + 1
      k = k + kstep
    END DO
  END IF
END SUBROUTINE SSPDI
