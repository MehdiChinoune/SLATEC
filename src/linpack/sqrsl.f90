!DECK SQRSL
SUBROUTINE SQRSL(X,Ldx,N,K,Qraux,Y,Qy,Qty,B,Rsd,Xb,Job,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SQRSL
  !***PURPOSE  Apply the output of SQRDC to compute coordinate transfor-
  !            mations, projections, and least squares solutions.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D9, D2A1
  !***TYPE      SINGLE PRECISION (SQRSL-S, DQRSL-D, CQRSL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
  !             SOLVE
  !***AUTHOR  Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     SQRSL applies the output of SQRDC to compute coordinate
  !     transformations, projections, and least squares solutions.
  !     For K .LE. MIN(N,P), let XK be the matrix
  !
  !            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
  !
  !     formed from columns JPVT(1), ... ,JPVT(K) of the original
  !     N x P matrix X that was input to SQRDC (if no pivoting was
  !     done, XK consists of the first K columns of X in their
  !     original order).  SQRDC produces a factored orthogonal matrix Q
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
  !        X      REAL(LDX,P)
  !               X contains the output of SQRDC.
  !
  !        LDX    INTEGER
  !               LDX is the leading dimension of the array X.
  !
  !        N      INTEGER
  !               N is the number of rows of the matrix XK.  It must
  !               have the same value as N in SQRDC.
  !
  !        K      INTEGER
  !               K is the number of columns of the matrix XK.  K
  !               must not be greater than MIN(N,P), where P is the
  !               same as in the calling sequence to SQRDC.
  !
  !        QRAUX  REAL(P)
  !               QRAUX contains the auxiliary output from SQRDC.
  !
  !        Y      REAL(N)
  !               Y contains an N-vector that is to be manipulated
  !               by SQRSL.
  !
  !        JOB    INTEGER
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
  !        QY     REAL(N).
  !               QY contains Q*Y, if its computation has been
  !               requested.
  !
  !        QTY    REAL(N).
  !               QTY contains TRANS(Q)*Y, if its computation has
  !               been requested.  Here TRANS(Q) is the
  !               transpose of the matrix Q.
  !
  !        B      REAL(K)
  !               B contains the solution of the least squares problem
  !
  !                    minimize norm2(Y - XK*B),
  !
  !               if its computation has been requested.  (Note that
  !               if pivoting was requested in SQRDC, the J-th
  !               component of B will be associated with column JPVT(J)
  !               of the original matrix X that was input into SQRDC.)
  !
  !        RSD    REAL(N).
  !               RSD contains the least squares residual Y - XK*B,
  !               if its computation has been requested.  RSD is
  !               also the orthogonal projection of Y onto the
  !               orthogonal complement of the column space of XK.
  !
  !        XB     REAL(N).
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
  !          CALL SQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
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
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  SAXPY, SCOPY, SDOT
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SQRSL
  INTEGER Ldx, N, K, Job, Info
  REAL X(Ldx,*), Qraux(*), Y(*), Qy(*), Qty(*), B(*), Rsd(*), Xb(*)
  !
  INTEGER i, j, jj, ju, kp1
  REAL SDOT, t, temp
  LOGICAL cb, cqy, cqty, cr, cxb
  !***FIRST EXECUTABLE STATEMENT  SQRSL
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
    IF ( cqy ) CALL SCOPY(N,Y,1,Qy,1)
    IF ( cqty ) CALL SCOPY(N,Y,1,Qty,1)
    IF ( cqy ) THEN
      !
      !           COMPUTE QY.
      !
      DO jj = 1, ju
        j = ju - jj + 1
        IF ( Qraux(j)/=0.0E0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          t = -SDOT(N-j+1,X(j,j),1,Qy(j),1)/X(j,j)
          CALL SAXPY(N-j+1,t,X(j,j),1,Qy(j),1)
          X(j,j) = temp
        ENDIF
      ENDDO
    ENDIF
    IF ( cqty ) THEN
      !
      !           COMPUTE TRANS(Q)*Y.
      !
      DO j = 1, ju
        IF ( Qraux(j)/=0.0E0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          t = -SDOT(N-j+1,X(j,j),1,Qty(j),1)/X(j,j)
          CALL SAXPY(N-j+1,t,X(j,j),1,Qty(j),1)
          X(j,j) = temp
        ENDIF
      ENDDO
    ENDIF
    !
    !        SET UP TO COMPUTE B, RSD, OR XB.
    !
    IF ( cb ) CALL SCOPY(K,Qty,1,B,1)
    kp1 = K + 1
    IF ( cxb ) CALL SCOPY(K,Qty,1,Xb,1)
    IF ( cr.AND.K<N ) CALL SCOPY(N-K,Qty(kp1),1,Rsd(kp1),1)
    IF ( .NOT.(.NOT.cxb.OR.kp1>N) ) THEN
      DO i = kp1, N
        Xb(i) = 0.0E0
      ENDDO
    ENDIF
    IF ( cr ) THEN
      DO i = 1, K
        Rsd(i) = 0.0E0
      ENDDO
    ENDIF
    IF ( cb ) THEN
      !
      !           COMPUTE B.
      !
      DO jj = 1, K
        j = K - jj + 1
        IF ( X(j,j)/=0.0E0 ) THEN
          B(j) = B(j)/X(j,j)
          IF ( j/=1 ) THEN
            t = -B(j)
            CALL SAXPY(j-1,t,X(1,j),1,B,1)
          ENDIF
        ELSE
          Info = j
          EXIT
        ENDIF
      ENDDO
    ENDIF
    IF ( .NOT.(.NOT.cr.AND..NOT.cxb) ) THEN
      !
      !           COMPUTE RSD OR XB AS REQUIRED.
      !
      DO jj = 1, ju
        j = ju - jj + 1
        IF ( Qraux(j)/=0.0E0 ) THEN
          temp = X(j,j)
          X(j,j) = Qraux(j)
          IF ( cr ) THEN
            t = -SDOT(N-j+1,X(j,j),1,Rsd(j),1)/X(j,j)
            CALL SAXPY(N-j+1,t,X(j,j),1,Rsd(j),1)
          ENDIF
          IF ( cxb ) THEN
            t = -SDOT(N-j+1,X(j,j),1,Xb(j),1)/X(j,j)
            CALL SAXPY(N-j+1,t,X(j,j),1,Xb(j),1)
          ENDIF
          X(j,j) = temp
        ENDIF
      ENDDO
    ENDIF
  ELSE
    IF ( cqy ) Qy(1) = Y(1)
    IF ( cqty ) Qty(1) = Y(1)
    IF ( cxb ) Xb(1) = Y(1)
    IF ( cb ) THEN
      IF ( X(1,1)/=0.0E0 ) THEN
        B(1) = Y(1)/X(1,1)
      ELSE
        Info = 1
      ENDIF
    ENDIF
    IF ( cr ) Rsd(1) = 0.0E0
  ENDIF
END SUBROUTINE SQRSL
