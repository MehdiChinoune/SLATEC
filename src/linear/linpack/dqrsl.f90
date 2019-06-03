!** DQRSL
SUBROUTINE DQRSL(X,Ldx,N,K,Qraux,Y,Qy,Qty,B,Rsd,Xb,Job,Info)
  !>
  !  Apply the output of DQRDC to compute coordinate transfor-
  !            mations, projections, and least squares solutions.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D9, D2A1
  !***
  ! **Type:**      DOUBLE PRECISION (SQRSL-S, DQRSL-D, CQRSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
  !             SOLVE
  !***
  ! **Author:**  Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     DQRSL applies the output of DQRDC to compute coordinate
  !     transformations, projections, and least squares solutions.
  !     For K .LE. MIN(N,P), let XK be the matrix
  !
  !            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
  !
  !     formed from columns JPVT(1), ... ,JPVT(K) of the original
  !     N X P matrix X that was input to DQRDC (if no pivoting was
  !     done, XK consists of the first K columns of X in their
  !     original order).  DQRDC produces a factored orthogonal matrix Q
  !     and an upper triangular matrix R such that
  !
  !              XK = Q * (R)
  !                       (0)
  !
  !     This information is contained in coded form in the arrays
  !     X and QRAUX.
  !
  !     On Entry
  !
  !        X      DOUBLE PRECISION(LDX,P).
  !               X contains the output of DQRDC.
  !
  !        LDX    INTEGER.
  !               LDX is the leading dimension of the array X.
  !
  !        N      INTEGER.
  !               N is the number of rows of the matrix XK.  It must
  !               have the same value as N in DQRDC.
  !
  !        K      INTEGER.
  !               K is the number of columns of the matrix XK.  K
  !               must not be greater than MIN(N,P), where P is the
  !               same as in the calling sequence to DQRDC.
  !
  !        QRAUX  DOUBLE PRECISION(P).
  !               QRAUX contains the auxiliary output from DQRDC.
  !
  !        Y      DOUBLE PRECISION(N)
  !               Y contains an N-vector that is to be manipulated
  !               by DQRSL.
  !
  !        JOB    INTEGER.
  !               JOB specifies what is to be computed.  JOB has
  !               the decimal expansion ABCDE, with the following
  !               meaning.
  !
  !                    If A .NE. 0, compute QY.
  !                    If B,C,D, or E .NE. 0, compute QTY.
  !                    If C .NE. 0, compute B.
  !                    If D .NE. 0, compute RSD.
  !                    If E .NE. 0, compute XB.
  !
  !               Note that a request to compute B, RSD, or XB
  !               automatically triggers the computation of QTY, for
  !               which an array must be provided in the calling
  !               sequence.
  !
  !     On Return
  !
  !        QY     DOUBLE PRECISION(N).
  !               QY contains Q*Y, if its computation has been
  !               requested.
  !
  !        QTY    DOUBLE PRECISION(N).
  !               QTY contains TRANS(Q)*Y, if its computation has
  !               been requested.  Here TRANS(Q) is the
  !               transpose of the matrix Q.
  !
  !        B      DOUBLE PRECISION(K)
  !               B contains the solution of the least squares problem
  !
  !                    minimize norm2(Y - XK*B),
  !
  !               if its computation has been requested.  (Note that
  !               if pivoting was requested in DQRDC, the J-th
  !               component of B will be associated with column JPVT(J)
  !               of the original matrix X that was input into DQRDC.)
  !
  !        RSD    DOUBLE PRECISION(N).
  !               RSD contains the least squares residual Y - XK*B,
  !               if its computation has been requested.  RSD is
  !               also the orthogonal projection of Y onto the
  !               orthogonal complement of the column space of XK.
  !
  !        XB     DOUBLE PRECISION(N).
  !               XB contains the least squares approximation XK*B,
  !               if its computation has been requested.  XB is also
  !               the orthogonal projection of Y onto the column space
  !               of X.
  !
  !        INFO   INTEGER.
  !               INFO is zero unless the computation of B has
  !               been requested and R is exactly singular.  In
  !               this case, INFO is the index of the first zero
  !               diagonal element of R and B is left unaltered.
  !
  !     The parameters QY, QTY, B, RSD, and XB are not referenced
  !     if their computation is not requested and in this case
  !     can be replaced by dummy variables in the calling program.
  !     To save storage, the user may in some cases use the same
  !     array for different parameters in the calling sequence.  A
  !     frequently occurring example is when one wishes to compute
  !     any of B, RSD, or XB and does not need Y or QTY.  In this
  !     case one may identify Y, QTY, and one of B, RSD, or XB, while
  !     providing separate arrays for anything else that is to be
  !     computed.  Thus the calling sequence
  !
  !          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
  !
  !     will result in the computation of B and RSD, with RSD
  !     overwriting Y.  More generally, each item in the following
  !     list contains groups of permissible identifications for
  !     a single calling sequence.
  !
  !          1. (Y,QTY,B) (RSD) (XB) (QY)
  !
  !          2. (Y,QTY,RSD) (B) (XB) (QY)
  !
  !          3. (Y,QTY,XB) (B) (RSD) (QY)
  !
  !          4. (Y,QY) (QTY,B) (RSD) (XB)
  !
  !          5. (Y,QY) (QTY,RSD) (B) (XB)
  !
  !          6. (Y,QY) (QTY,XB) (B) (RSD)
  !
  !     In any group the value returned in the array allocated to
  !     the group corresponds to the last member of the group.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DCOPY, DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Ldx, N, K, Job, Info
  REAL(DP) :: X(Ldx,*), Qraux(*), Y(*), Qy(*), Qty(*), B(*), Rsd(*), Xb(*)
  !
  INTEGER i, j, jj, ju, kp1
  REAL(DP) :: t, temp
  LOGICAL cb, cqy, cqty, cr, cxb
  !* FIRST EXECUTABLE STATEMENT  DQRSL
  !
  !     SET INFO FLAG.
  !
  Info = 0
  !
  !     DETERMINE WHAT IS TO BE COMPUTED.
  !
  cqy = Job/10000/=0
  cqty = MOD(Job,10000)/=0
  cb = MOD(Job,1000)/100/=0
  cr = MOD(Job,100)/10/=0
  cxb = MOD(Job,10)/=0
  ju = MIN(K,N-1)
  !
  !     SPECIAL ACTION WHEN N=1.
  !
  IF ( ju/=0 ) THEN
    !
    !        SET UP TO COMPUTE QY OR QTY.
    !
    IF ( cqy ) CALL DCOPY(N,Y,1,Qy,1)
    IF ( cqty ) CALL DCOPY(N,Y,1,Qty,1)
    IF ( cqy ) THEN
      !
      !           COMPUTE QY.
      !
      DO jj = 1, ju
        j = ju - jj + 1
        IF ( Qraux(j)/=0.0D0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          t = -DDOT(N-j+1,X(j,j),1,Qy(j),1)/X(j,j)
          CALL DAXPY(N-j+1,t,X(j,j),1,Qy(j),1)
          X(j,j) = temp
        END IF
      END DO
    END IF
    IF ( cqty ) THEN
      !
      !           COMPUTE TRANS(Q)*Y.
      !
      DO j = 1, ju
        IF ( Qraux(j)/=0.0D0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          t = -DDOT(N-j+1,X(j,j),1,Qty(j),1)/X(j,j)
          CALL DAXPY(N-j+1,t,X(j,j),1,Qty(j),1)
          X(j,j) = temp
        END IF
      END DO
    END IF
    !
    !        SET UP TO COMPUTE B, RSD, OR XB.
    !
    IF ( cb ) CALL DCOPY(K,Qty,1,B,1)
    kp1 = K + 1
    IF ( cxb ) CALL DCOPY(K,Qty,1,Xb,1)
    IF ( cr.AND.K<N ) CALL DCOPY(N-K,Qty(kp1),1,Rsd(kp1),1)
    IF ( .NOT.(.NOT.cxb.OR.kp1>N) ) THEN
      DO i = kp1, N
        Xb(i) = 0.0D0
      END DO
    END IF
    IF ( cr ) THEN
      DO i = 1, K
        Rsd(i) = 0.0D0
      END DO
    END IF
    IF ( cb ) THEN
      !
      !           COMPUTE B.
      !
      DO jj = 1, K
        j = K - jj + 1
        IF ( X(j,j)/=0.0D0 ) THEN
          B(j) = B(j)/X(j,j)
          IF ( j/=1 ) THEN
            t = -B(j)
            CALL DAXPY(j-1,t,X(1,j),1,B,1)
          END IF
        ELSE
          Info = j
          EXIT
        END IF
      END DO
    END IF
    IF ( .NOT.(.NOT.cr.AND..NOT.cxb) ) THEN
      !
      !           COMPUTE RSD OR XB AS REQUIRED.
      !
      DO jj = 1, ju
        j = ju - jj + 1
        IF ( Qraux(j)/=0.0D0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          IF ( cr ) THEN
            t = -DDOT(N-j+1,X(j,j),1,Rsd(j),1)/X(j,j)
            CALL DAXPY(N-j+1,t,X(j,j),1,Rsd(j),1)
          END IF
          IF ( cxb ) THEN
            t = -DDOT(N-j+1,X(j,j),1,Xb(j),1)/X(j,j)
            CALL DAXPY(N-j+1,t,X(j,j),1,Xb(j),1)
          END IF
          X(j,j) = temp
        END IF
      END DO
    END IF
  ELSE
    IF ( cqy ) Qy(1) = Y(1)
    IF ( cqty ) Qty(1) = Y(1)
    IF ( cxb ) Xb(1) = Y(1)
    IF ( cb ) THEN
      IF ( X(1,1)/=0.0D0 ) THEN
        B(1) = Y(1)/X(1,1)
      ELSE
        Info = 1
      END IF
    END IF
    IF ( cr ) Rsd(1) = 0.0D0
  END IF
END SUBROUTINE DQRSL
