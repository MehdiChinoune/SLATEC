!** DQRDC
SUBROUTINE DQRDC(X,Ldx,N,P,Qraux,Jpvt,Work,Job)
  IMPLICIT NONE
  !>
  !***
  !  Use Householder transformations to compute the QR
  !            factorization of an N by P matrix.  Column pivoting is a
  !            users option.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D5
  !***
  ! **Type:**      DOUBLE PRECISION (SQRDC-S, DQRDC-D, CQRDC-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
  !             QR DECOMPOSITION
  !***
  ! **Author:**  Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     DQRDC uses Householder transformations to compute the QR
  !     factorization of an N by P matrix X.  Column pivoting
  !     based on the 2-norms of the reduced columns may be
  !     performed at the user's option.
  !
  !     On Entry
  !
  !        X       DOUBLE PRECISION(LDX,P), where LDX .GE. N.
  !                X contains the matrix whose decomposition is to be
  !                computed.
  !
  !        LDX     INTEGER.
  !                LDX is the leading dimension of the array X.
  !
  !        N       INTEGER.
  !                N is the number of rows of the matrix X.
  !
  !        P       INTEGER.
  !                P is the number of columns of the matrix X.
  !
  !        JPVT    INTEGER(P).
  !                JPVT contains integers that control the selection
  !                of the pivot columns.  The K-th column X(K) of X
  !                is placed in one of three classes according to the
  !                value of JPVT(K).
  !
  !                   If JPVT(K) .GT. 0, then X(K) is an initial
  !                                      column.
  !
  !                   If JPVT(K) .EQ. 0, then X(K) is a free column.
  !
  !                   If JPVT(K) .LT. 0, then X(K) is a final column.
  !
  !                Before the decomposition is computed, initial columns
  !                are moved to the beginning of the array X and final
  !                columns to the end.  Both initial and final columns
  !                are frozen in place during the computation and only
  !                free columns are moved.  At the K-th stage of the
  !                reduction, if X(K) is occupied by a free column
  !                it is interchanged with the free column of largest
  !                reduced norm.  JPVT is not referenced if
  !                JOB .EQ. 0.
  !
  !        WORK    DOUBLE PRECISION(P).
  !                WORK is a work array.  WORK is not referenced if
  !                JOB .EQ. 0.
  !
  !        JOB     INTEGER.
  !                JOB is an integer that initiates column pivoting.
  !                If JOB .EQ. 0, no pivoting is done.
  !                If JOB .NE. 0, pivoting is done.
  !
  !     On Return
  !
  !        X       X contains in its upper triangle the upper
  !                triangular matrix R of the QR factorization.
  !                Below its diagonal X contains information from
  !                which the orthogonal part of the decomposition
  !                can be recovered.  Note that if pivoting has
  !                been requested, the decomposition is not that
  !                of the original matrix X but that of X
  !                with its columns permuted as described by JPVT.
  !
  !        QRAUX   DOUBLE PRECISION(P).
  !                QRAUX contains further information required to recover
  !                the orthogonal part of the decomposition.
  !
  !        JPVT    JPVT(K) contains the index of the column of the
  !                original matrix that has been interchanged into
  !                the K-th column, if pivoting was requested.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DDOT, DNRM2, DSCAL, DSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER jj
  INTEGER Ldx, N, P, Job
  INTEGER Jpvt(*)
  REAL(8) :: X(Ldx,*), Qraux(*), Work(*)
  !
  INTEGER j, jp, l, lp1, lup, maxj, pl, pu
  REAL(8) :: maxnrm, DNRM2, tt
  REAL(8) :: DDOT, nrmxl, t
  LOGICAL negj, swapj
  !
  !* FIRST EXECUTABLE STATEMENT  DQRDC
  pl = 1
  pu = 0
  IF ( Job/=0 ) THEN
    !
    !        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
    !        ACCORDING TO JPVT.
    !
    DO j = 1, P
      swapj = Jpvt(j)>0
      negj = Jpvt(j)<0
      Jpvt(j) = j
      IF ( negj ) Jpvt(j) = -j
      IF ( swapj ) THEN
        IF ( j/=pl ) CALL DSWAP(N,X(1,pl),1,X(1,j),1)
        Jpvt(j) = Jpvt(pl)
        Jpvt(pl) = j
        pl = pl + 1
      END IF
    END DO
    pu = P
    DO jj = 1, P
      j = P - jj + 1
      IF ( Jpvt(j)<0 ) THEN
        Jpvt(j) = -Jpvt(j)
        IF ( j/=pu ) THEN
          CALL DSWAP(N,X(1,pu),1,X(1,j),1)
          jp = Jpvt(pu)
          Jpvt(pu) = Jpvt(j)
          Jpvt(j) = jp
        END IF
        pu = pu - 1
      END IF
    END DO
  END IF
  !
  !     COMPUTE THE NORMS OF THE FREE COLUMNS.
  !
  IF ( pu>=pl ) THEN
    DO j = pl, pu
      Qraux(j) = DNRM2(N,X(1,j),1)
      Work(j) = Qraux(j)
    END DO
  END IF
  !
  !     PERFORM THE HOUSEHOLDER REDUCTION OF X.
  !
  lup = MIN(N,P)
  DO l = 1, lup
    IF ( l>=pl.AND.l<pu ) THEN
      !
      !           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
      !           INTO THE PIVOT POSITION.
      !
      maxnrm = 0.0D0
      maxj = l
      DO j = l, pu
        IF ( Qraux(j)>maxnrm ) THEN
          maxnrm = Qraux(j)
          maxj = j
        END IF
      END DO
      IF ( maxj/=l ) THEN
        CALL DSWAP(N,X(1,l),1,X(1,maxj),1)
        Qraux(maxj) = Qraux(l)
        Work(maxj) = Work(l)
        jp = Jpvt(maxj)
        Jpvt(maxj) = Jpvt(l)
        Jpvt(l) = jp
      END IF
    END IF
    Qraux(l) = 0.0D0
    IF ( l/=N ) THEN
      !
      !           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
      !
      nrmxl = DNRM2(N-l+1,X(l,l),1)
      IF ( nrmxl/=0.0D0 ) THEN
        IF ( X(l,l)/=0.0D0 ) nrmxl = SIGN(nrmxl,X(l,l))
        CALL DSCAL(N-l+1,1.0D0/nrmxl,X(l,l),1)
        X(l,l) = 1.0D0 + X(l,l)
        !
        !              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
        !              UPDATING THE NORMS.
        !
        lp1 = l + 1
        IF ( P>=lp1 ) THEN
          DO j = lp1, P
            t = -DDOT(N-l+1,X(l,l),1,X(l,j),1)/X(l,l)
            CALL DAXPY(N-l+1,t,X(l,l),1,X(l,j),1)
            IF ( j>=pl.AND.j<=pu ) THEN
              IF ( Qraux(j)/=0.0D0 ) THEN
                tt = 1.0D0 - (ABS(X(l,j))/Qraux(j))**2
                tt = MAX(tt,0.0D0)
                t = tt
                tt = 1.0D0 + 0.05D0*tt*(Qraux(j)/Work(j))**2
                IF ( tt==1.0D0 ) THEN
                  Qraux(j) = DNRM2(N-l,X(l+1,j),1)
                  Work(j) = Qraux(j)
                ELSE
                  Qraux(j) = Qraux(j)*SQRT(t)
                END IF
              END IF
            END IF
          END DO
        END IF
        !
        !              SAVE THE TRANSFORMATION.
        !
        Qraux(l) = X(l,l)
        X(l,l) = -nrmxl
      END IF
    END IF
  END DO
END SUBROUTINE DQRDC