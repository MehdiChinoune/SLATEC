!DECK SSVDC
SUBROUTINE SSVDC(X,Ldx,N,P,S,E,U,Ldu,V,Ldv,Work,Job,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SSVDC
  !***PURPOSE  Perform the singular value decomposition of a rectangular
  !            matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D6
  !***TYPE      SINGLE PRECISION (SSVDC-S, DSVDC-D, CSVDC-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX,
  !             SINGULAR VALUE DECOMPOSITION
  !***AUTHOR  Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     SSVDC is a subroutine to reduce a real NxP matrix X by orthogonal
  !     transformations U and V to diagonal form.  The elements S(I) are
  !     the singular values of X.  The columns of U are the corresponding
  !     left singular vectors, and the columns of V the right singular
  !     vectors.
  !
  !     On Entry
  !
  !         X         REAL(LDX,P), where LDX .GE. N.
  !                   X contains the matrix whose singular value
  !                   decomposition is to be computed.  X is
  !                   destroyed by SSVDC.
  !
  !         LDX       INTEGER
  !                   LDX is the leading dimension of the array X.
  !
  !         N         INTEGER
  !                   N is the number of rows of the matrix X.
  !
  !         P         INTEGER
  !                   P is the number of columns of the matrix X.
  !
  !         LDU       INTEGER
  !                   LDU is the leading dimension of the array U.
  !                   (See below).
  !
  !         LDV       INTEGER
  !                   LDV is the leading dimension of the array V.
  !                   (See below).
  !
  !         WORK      REAL(N)
  !                   work is a scratch array.
  !
  !         JOB       INTEGER
  !                   JOB controls the computation of the singular
  !                   vectors.  It has the decimal expansion AB
  !                   with the following meaning
  !
  !                        A .EQ. 0  Do not compute the left singular
  !                                  vectors.
  !                        A .EQ. 1  Return the N left singular vectors
  !                                  in U.
  !                        A .GE. 2  Return the first MIN(N,P) singular
  !                                  vectors in U.
  !                        B .EQ. 0  Do not compute the right singular
  !                                  vectors.
  !                        B .EQ. 1  Return the right singular vectors
  !                                  in V.
  !
  !     On Return
  !
  !         S         REAL(MM), where MM=MIN(N+1,P).
  !                   The first MIN(N,P) entries of S contain the
  !                   singular values of X arranged in descending
  !                   order of magnitude.
  !
  !         E         REAL(P).
  !                   E ordinarily contains zeros.  However, see the
  !                   discussion of INFO for exceptions.
  !
  !         U         REAL(LDU,K), where LDU .GE. N.  If JOBA .EQ. 1, then
  !                                   K .EQ. N.  If JOBA .GE. 2, then
  !                                   K .EQ. MIN(N,P).
  !                   U contains the matrix of right singular vectors.
  !                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
  !                   or if JOBA .EQ. 2, then U may be identified with X
  !                   in the subroutine call.
  !
  !         V         REAL(LDV,P), where LDV .GE. P.
  !                   V contains the matrix of right singular vectors.
  !                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
  !                   then V may be identified with X in the
  !                   subroutine call.
  !
  !         INFO      INTEGER.
  !                   the singular values (and their corresponding
  !                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
  !                   are correct (here M=MIN(N,P)).  Thus if
  !                   INFO .EQ. 0, all the singular values and their
  !                   vectors are correct.  In any event, the matrix
  !                   B = TRANS(U)*X*V is the bidiagonal matrix
  !                   with the elements of S on its diagonal and the
  !                   elements of E on its super-diagonal (TRANS(U)
  !                   is the transpose of U).  Thus the singular
  !                   values of X and B are the same.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SROT, SROTG, SSCAL, SSWAP
  !***REVISION HISTORY  (YYMMDD)
  !   790319  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SSVDC
  INTEGER Ldx, N, P, Ldu, Ldv, Job, Info
  REAL X(Ldx,*), S(*), E(*), U(Ldu,*), V(Ldv,*), Work(*)
  !
  !
  INTEGER i, iter, j, jobu, k, kase, kk, l, ll, lls, lm1, lp1, &
    ls, lu, m, maxit, mm, mm1, mp1, nct, nctp1, ncu, nrt, &
    nrtp1
  REAL SDOT, t
  REAL b, c, cs, el, emm1, f, g, SNRM2, scale, shift, sl, sm, &
    sn, smm1, t1, test, ztest
  LOGICAL wantu, wantv
  !***FIRST EXECUTABLE STATEMENT  SSVDC
  !
  !     SET THE MAXIMUM NUMBER OF ITERATIONS.
  !
  maxit = 30
  !
  !     DETERMINE WHAT IS TO BE COMPUTED.
  !
  wantu = .FALSE.
  wantv = .FALSE.
  jobu = MOD(Job,100)/10
  ncu = N
  IF ( jobu>1 ) ncu = MIN(N,P)
  IF ( jobu/=0 ) wantu = .TRUE.
  IF ( MOD(Job,10)/=0 ) wantv = .TRUE.
  !
  !     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
  !     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
  !
  Info = 0
  nct = MIN(N-1,P)
  nrt = MAX(0,MIN(P-2,N))
  lu = MAX(nct,nrt)
  IF ( lu>=1 ) THEN
    DO l = 1, lu
      lp1 = l + 1
      IF ( l<=nct ) THEN
        !
        !           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
        !           PLACE THE L-TH DIAGONAL IN S(L).
        !
        S(l) = SNRM2(N-l+1,X(l,l),1)
        IF ( S(l)/=0.0E0 ) THEN
          IF ( X(l,l)/=0.0E0 ) S(l) = SIGN(S(l),X(l,l))
          CALL SSCAL(N-l+1,1.0E0/S(l),X(l,l),1)
          X(l,l) = 1.0E0 + X(l,l)
        ENDIF
        S(l) = -S(l)
      ENDIF
      IF ( P>=lp1 ) THEN
        DO j = lp1, P
          IF ( l<=nct ) THEN
            IF ( S(l)/=0.0E0 ) THEN
              !
              !              APPLY THE TRANSFORMATION.
              !
              t = -SDOT(N-l+1,X(l,l),1,X(l,j),1)/X(l,l)
              CALL SAXPY(N-l+1,t,X(l,l),1,X(l,j),1)
            ENDIF
          ENDIF
          !
          !           PLACE THE L-TH ROW OF X INTO  E FOR THE
          !           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
          !
          E(j) = X(l,j)
        ENDDO
      ENDIF
      IF ( .NOT.(.NOT.wantu.OR.l>nct) ) THEN
        !
        !           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
        !           MULTIPLICATION.
        !
        DO i = l, N
          U(i,l) = X(i,l)
        ENDDO
      ENDIF
      IF ( l<=nrt ) THEN
        !
        !           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
        !           L-TH SUPER-DIAGONAL IN E(L).
        !
        E(l) = SNRM2(P-l,E(lp1),1)
        IF ( E(l)/=0.0E0 ) THEN
          IF ( E(lp1)/=0.0E0 ) E(l) = SIGN(E(l),E(lp1))
          CALL SSCAL(P-l,1.0E0/E(l),E(lp1),1)
          E(lp1) = 1.0E0 + E(lp1)
        ENDIF
        E(l) = -E(l)
        IF ( lp1<=N.AND.E(l)/=0.0E0 ) THEN
          !
          !              APPLY THE TRANSFORMATION.
          !
          DO i = lp1, N
            Work(i) = 0.0E0
          ENDDO
          DO j = lp1, P
            CALL SAXPY(N-l,E(j),X(lp1,j),1,Work(lp1),1)
          ENDDO
          DO j = lp1, P
            CALL SAXPY(N-l,-E(j)/E(lp1),Work(lp1),1,X(lp1,j),1)
          ENDDO
        ENDIF
        IF ( wantv ) THEN
          !
          !              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
          !              BACK MULTIPLICATION.
          !
          DO i = lp1, P
            V(i,l) = E(i)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  !
  !     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
  !
  m = MIN(P,N+1)
  nctp1 = nct + 1
  nrtp1 = nrt + 1
  IF ( nct<P ) S(nctp1) = X(nctp1,nctp1)
  IF ( N<m ) S(m) = 0.0E0
  IF ( nrtp1<m ) E(nrtp1) = X(nrtp1,m)
  E(m) = 0.0E0
  !
  !     IF REQUIRED, GENERATE U.
  !
  IF ( wantu ) THEN
    IF ( ncu>=nctp1 ) THEN
      DO j = nctp1, ncu
        DO i = 1, N
          U(i,j) = 0.0E0
        ENDDO
        U(j,j) = 1.0E0
      ENDDO
    ENDIF
    IF ( nct>=1 ) THEN
      DO ll = 1, nct
        l = nct - ll + 1
        IF ( S(l)==0.0E0 ) THEN
          DO i = 1, N
            U(i,l) = 0.0E0
          ENDDO
          U(l,l) = 1.0E0
        ELSE
          lp1 = l + 1
          IF ( ncu>=lp1 ) THEN
            DO j = lp1, ncu
              t = -SDOT(N-l+1,U(l,l),1,U(l,j),1)/U(l,l)
              CALL SAXPY(N-l+1,t,U(l,l),1,U(l,j),1)
            ENDDO
          ENDIF
          CALL SSCAL(N-l+1,-1.0E0,U(l,l),1)
          U(l,l) = 1.0E0 + U(l,l)
          lm1 = l - 1
          IF ( lm1>=1 ) THEN
            DO i = 1, lm1
              U(i,l) = 0.0E0
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  !
  !     IF IT IS REQUIRED, GENERATE V.
  !
  IF ( wantv ) THEN
    DO ll = 1, P
      l = P - ll + 1
      lp1 = l + 1
      IF ( l<=nrt ) THEN
        IF ( E(l)/=0.0E0 ) THEN
          DO j = lp1, P
            t = -SDOT(P-l,V(lp1,l),1,V(lp1,j),1)/V(lp1,l)
            CALL SAXPY(P-l,t,V(lp1,l),1,V(lp1,j),1)
          ENDDO
        ENDIF
      ENDIF
      DO i = 1, P
        V(i,l) = 0.0E0
      ENDDO
      V(l,l) = 1.0E0
    ENDDO
  ENDIF
  !
  !     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
  !
  mm = m
  iter = 0
  !
  !        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
  !
  DO WHILE ( m/=0 )
    !
    !        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
    !        FLAG AND RETURN.
    !
    IF ( iter<maxit ) THEN
      !
      !        THIS SECTION OF THE PROGRAM INSPECTS FOR
      !        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
      !        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
      !
      !           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
      !           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
      !           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
      !                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
      !           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
      !
      DO ll = 1, m
        l = m - ll
        IF ( l==0 ) EXIT
        test = ABS(S(l)) + ABS(S(l+1))
        ztest = test + ABS(E(l))
        IF ( ztest==test ) THEN
          E(l) = 0.0E0
          EXIT
        ENDIF
      ENDDO
      IF ( l/=m-1 ) THEN
        lp1 = l + 1
        mp1 = m + 1
        DO lls = lp1, mp1
          ls = m - lls + lp1
          IF ( ls==l ) EXIT
          test = 0.0E0
          IF ( ls/=m ) test = test + ABS(E(ls))
          IF ( ls/=l+1 ) test = test + ABS(E(ls-1))
          ztest = test + ABS(S(ls))
          IF ( ztest==test ) THEN
            S(ls) = 0.0E0
            EXIT
          ENDIF
        ENDDO
        IF ( ls==l ) THEN
          kase = 3
        ELSEIF ( ls/=m ) THEN
          kase = 2
          l = ls
        ELSE
          kase = 1
        ENDIF
      ELSE
        kase = 4
      ENDIF
      l = l + 1
      !
      !        PERFORM THE TASK INDICATED BY KASE.
      !
      SELECT CASE (kase)
        CASE (2)
          !
          !        SPLIT AT NEGLIGIBLE S(L).
          !
          f = E(l-1)
          E(l-1) = 0.0E0
          DO k = l, m
            t1 = S(k)
            CALL SROTG(t1,f,cs,sn)
            S(k) = t1
            f = -sn*E(k)
            E(k) = cs*E(k)
            IF ( wantu ) CALL SROT(N,U(1,k),1,U(1,l-1),1,cs,sn)
          ENDDO
        CASE (3)
          !
          !        PERFORM ONE QR STEP.
          !
          !
          !           CALCULATE THE SHIFT.
          !
          scale = MAX(ABS(S(m)),ABS(S(m-1)),ABS(E(m-1)),ABS(S(l)),ABS(E(l)))
          sm = S(m)/scale
          smm1 = S(m-1)/scale
          emm1 = E(m-1)/scale
          sl = S(l)/scale
          el = E(l)/scale
          b = ((smm1+sm)*(smm1-sm)+emm1**2)/2.0E0
          c = (sm*emm1)**2
          shift = 0.0E0
          IF ( b/=0.0E0.OR.c/=0.0E0 ) THEN
            shift = SQRT(b**2+c)
            IF ( b<0.0E0 ) shift = -shift
            shift = c/(b+shift)
          ENDIF
          f = (sl+sm)*(sl-sm) - shift
          g = sl*el
          !
          !           CHASE ZEROS.
          !
          mm1 = m - 1
          DO k = l, mm1
            CALL SROTG(f,g,cs,sn)
            IF ( k/=l ) E(k-1) = f
            f = cs*S(k) + sn*E(k)
            E(k) = cs*E(k) - sn*S(k)
            g = sn*S(k+1)
            S(k+1) = cs*S(k+1)
            IF ( wantv ) CALL SROT(P,V(1,k),1,V(1,k+1),1,cs,sn)
            CALL SROTG(f,g,cs,sn)
            S(k) = f
            f = cs*E(k) + sn*S(k+1)
            S(k+1) = -sn*E(k) + cs*S(k+1)
            g = sn*E(k+1)
            E(k+1) = cs*E(k+1)
            IF ( wantu.AND.k<N ) CALL SROT(N,U(1,k),1,U(1,k+1),1,cs,sn)
          ENDDO
          E(m-1) = f
          iter = iter + 1
        CASE (4)
          !
          !        CONVERGENCE.
          !
          !
          !           MAKE THE SINGULAR VALUE  POSITIVE.
          !
          IF ( S(l)<0.0E0 ) THEN
            S(l) = -S(l)
            IF ( wantv ) CALL SSCAL(P,-1.0E0,V(1,l),1)
          ENDIF
          !
          !           ORDER THE SINGULAR VALUE.
          !
          DO WHILE ( l/=mm )
            IF ( S(l)>=S(l+1) ) EXIT
            t = S(l)
            S(l) = S(l+1)
            S(l+1) = t
            IF ( wantv.AND.l<P ) CALL SSWAP(P,V(1,l),1,V(1,l+1),1)
            IF ( wantu.AND.l<N ) CALL SSWAP(N,U(1,l),1,U(1,l+1),1)
            l = l + 1
          ENDDO
          GOTO 50
        CASE DEFAULT
          !
          !        DEFLATE NEGLIGIBLE S(M).
          !
          mm1 = m - 1
          f = E(m-1)
          E(m-1) = 0.0E0
          DO kk = l, mm1
            k = mm1 - kk + l
            t1 = S(k)
            CALL SROTG(t1,f,cs,sn)
            S(k) = t1
            IF ( k/=l ) THEN
              f = -sn*E(k-1)
              E(k-1) = cs*E(k-1)
            ENDIF
            IF ( wantv ) CALL SROT(P,V(1,k),1,V(1,m),1,cs,sn)
          ENDDO
      END SELECT
      CYCLE
    ELSE
      Info = m
      EXIT
    ENDIF
    50     iter = 0
    m = m - 1
  ENDDO
END SUBROUTINE SSVDC
