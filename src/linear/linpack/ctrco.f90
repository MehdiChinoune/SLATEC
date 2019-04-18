!** CTRCO
SUBROUTINE CTRCO(T,Ldt,N,Rcond,Z,Job)
  !>
  !  Estimate the condition number of a triangular matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2C3
  !***
  ! **Type:**      COMPLEX (STRCO-S, DTRCO-D, CTRCO-C)
  !***
  ! **Keywords:**  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
  !             TRIANGULAR MATRIX
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     CTRCO estimates the condition of a complex triangular matrix.
  !
  !     On Entry
  !
  !        T       COMPLEX(LDT,N)
  !                T contains the triangular matrix.  The zero
  !                elements of the matrix are not referenced, and
  !                the corresponding elements of the array can be
  !                used to store other information.
  !
  !        LDT     INTEGER
  !                LDT is the leading dimension of the array T.
  !
  !        N       INTEGER
  !                N is the order of the system.
  !
  !        JOB     INTEGER
  !                = 0         T  is lower triangular.
  !                = nonzero   T  is upper triangular.
  !
  !     On Return
  !
  !        RCOND   REAL
  !                an estimate of the reciprocal condition of  T .
  !                For the system  T*X = B, relative perturbations
  !                in  T  and  B  of size  EPSILON  may cause
  !                relative perturbations in  X  of size  EPSILON/RCOND .
  !                If  RCOND  is so small that the logical expression
  !                           1.0 + RCOND .EQ. 1.0
  !                is true, then  T  may be singular to working
  !                precision.  In particular,  RCOND  is zero  if
  !                exact singularity is detected or the estimate
  !                underflows.
  !
  !        Z       COMPLEX(N)
  !                a work vector whose contents are usually unimportant.
  !                If  T  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CAXPY, CSSCAL, SCASUM

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Ldt, N, Job
  COMPLEX T(Ldt,*), Z(*)
  REAL Rcond
  !
  COMPLEX w, wk, wkm, ek
  REAL tnorm, ynorm, s, sm
  INTEGER i1, j, j1, j2, k, kk, l
  LOGICAL lower
  !
  !* FIRST EXECUTABLE STATEMENT  CTRCO
  lower = Job==0
  !
  !     COMPUTE 1-NORM OF T
  !
  tnorm = 0.0E0
  DO j = 1, N
    l = j
    IF ( lower ) l = N + 1 - j
    i1 = 1
    IF ( lower ) i1 = j
    tnorm = MAX(tnorm,SCASUM(l,T(i1,j),1))
  END DO
  !
  !     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  CTRANS(T)*Y = E .
  !     CTRANS(T)  IS THE CONJUGATE TRANSPOSE OF T .
  !     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
  !     GROWTH IN THE ELEMENTS OF Y .
  !     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  !
  !     SOLVE CTRANS(T)*Y = E
  !
  ek = (1.0E0,0.0E0)
  DO j = 1, N
    Z(j) = (0.0E0,0.0E0)
  END DO
  DO kk = 1, N
    k = kk
    IF ( lower ) k = N + 1 - kk
    IF ( CABS1(Z(k))/=0.0E0 ) ek = CSIGN1(ek,-Z(k))
    IF ( CABS1(ek-Z(k))>CABS1(T(k,k)) ) THEN
      s = CABS1(T(k,k))/CABS1(ek-Z(k))
      CALL CSSCAL(N,s,Z,1)
      ek = CMPLX(s,0.0E0)*ek
    END IF
    wk = ek - Z(k)
    wkm = -ek - Z(k)
    s = CABS1(wk)
    sm = CABS1(wkm)
    IF ( CABS1(T(k,k))==0.0E0 ) THEN
      wk = (1.0E0,0.0E0)
      wkm = (1.0E0,0.0E0)
    ELSE
      wk = wk/CONJG(T(k,k))
      wkm = wkm/CONJG(T(k,k))
    END IF
    IF ( kk/=N ) THEN
      j1 = k + 1
      IF ( lower ) j1 = 1
      j2 = N
      IF ( lower ) j2 = k - 1
      DO j = j1, j2
        sm = sm + CABS1(Z(j)+wkm*CONJG(T(k,j)))
        Z(j) = Z(j) + wk*CONJG(T(k,j))
        s = s + CABS1(Z(j))
      END DO
      IF ( s<sm ) THEN
        w = wkm - wk
        wk = wkm
        DO j = j1, j2
          Z(j) = Z(j) + w*CONJG(T(k,j))
        END DO
      END IF
    END IF
    Z(k) = wk
  END DO
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  !
  ynorm = 1.0E0
  !
  !     SOLVE T*Z = Y
  !
  DO kk = 1, N
    k = N + 1 - kk
    IF ( lower ) k = kk
    IF ( CABS1(Z(k))>CABS1(T(k,k)) ) THEN
      s = CABS1(T(k,k))/CABS1(Z(k))
      CALL CSSCAL(N,s,Z,1)
      ynorm = s*ynorm
    END IF
    IF ( CABS1(T(k,k))/=0.0E0 ) Z(k) = Z(k)/T(k,k)
    IF ( CABS1(T(k,k))==0.0E0 ) Z(k) = (1.0E0,0.0E0)
    i1 = 1
    IF ( lower ) i1 = k + 1
    IF ( kk<N ) THEN
      w = -Z(k)
      CALL CAXPY(N-kk,w,T(i1,k),1,Z(i1),1)
    END IF
  END DO
  !     MAKE ZNORM = 1.0
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  IF ( tnorm/=0.0E0 ) Rcond = ynorm/tnorm
  IF ( tnorm==0.0E0 ) Rcond = 0.0E0
END SUBROUTINE CTRCO
