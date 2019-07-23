!** DWNLSM
PURE SUBROUTINE DWNLSM(W,Mdw,Mme,Ma,N,L,Prgopt,X,Rnorm,Mode,Ipivot,Itype,Wd,H,&
    Scalee,Z,Temp,D)
  !> Subsidiary to DWNNLS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (WNLSM-S, DWNLSM-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to DWNNLS.
  !     The documentation for DWNNLS has complete usage instructions.
  !
  !     In addition to the parameters discussed in the prologue to
  !     subroutine DWNNLS, the following work arrays are used in
  !     subroutine DWNLSM  (they are passed through the calling
  !     sequence from DWNNLS for purposes of variable dimensioning).
  !     Their contents will in general be of no interest to the user.
  !
  !     Variables of type REAL are DOUBLE PRECISION.
  !
  !         IPIVOT(*)
  !            An array of length N.  Upon completion it contains the
  !         pivoting information for the cols of W(*,*).
  !
  !         ITYPE(*)
  !            An array of length M which is used to keep track
  !         of the classification of the equations.  ITYPE(I)=0
  !         denotes equation I as an equality constraint.
  !         ITYPE(I)=1 denotes equation I as a least squares
  !         equation.
  !
  !         WD(*)
  !            An array of length N.  Upon completion it contains the
  !         dual solution vector.
  !
  !         H(*)
  !            An array of length N.  Upon completion it contains the
  !         pivot scalars of the Householder transformations performed
  !         in the case KRANK<L.
  !
  !         SCALE(*)
  !            An array of length M which is used by the subroutine
  !         to store the diagonal matrix of weights.
  !         These are used to apply the modified Givens
  !         transformations.
  !
  !         Z(*),TEMP(*)
  !            Working arrays of length N.
  !
  !         D(*)
  !            An array of length N that contains the
  !         column scaling for the matrix (E).
  !                                       (A)
  !
  !***
  ! **See also:**  DWNNLS
  !***
  ! **Routines called:**  D1MACH, DASUM, DAXPY, DCOPY, DH12, DNRM2, DROTM,
  !                    DROTMG, DSCAL, DSWAP, DWNLIT, IDAMAX, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and revised.  (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Fixed an error message.  (RWC)
  !   900604  DP version created from SP version.  (RWC)
  !   900911  Restriction on value of ALAMDA included.  (WRB)
  USE service, ONLY : D1MACH
  USE blas, ONLY : DAXPY, DROTM, DROTMG, DSWAP
  USE linear, ONLY : DH12
  !
  INTEGER, INTENT(IN) :: L, Ma, Mdw, Mme, N
  INTEGER, INTENT(OUT) :: Mode, Ipivot(N), Itype(Mme+Ma)
  REAL(DP), INTENT(IN) :: Prgopt(:)
  REAL(DP), INTENT(INOUT) :: W(Mdw,N+1)
  REAL(DP), INTENT(OUT) :: Rnorm, D(N), H(N), Scalee(Mme+Ma), Temp(N), Wd(N), X(N), Z(N)
  !
  INTEGER :: i, idope(3), imax, isol, itemp, iter, itmax, iwmax, j, jcon, jp, key, &
    krank, l1, last, link, m, me, next, niv, nlink, nopt, nsoln, ntimes
  REAL(DP) :: alamda, alpha, alsq, amax, blowup, bnorm, dope(3), eanorm, fac, sm, &
    sparam(5), t, tau, wmax, z2, zz
  LOGICAL :: done, feasbl, hitcon, pos
  !
  REAL(DP), PARAMETER :: drelpr = D1MACH(4)
  !* FIRST EXECUTABLE STATEMENT  DWNLSM
  !
  !     Set the nominal tolerance used in the code.
  !
  tau = SQRT(drelpr)
  !
  m = Ma + Mme
  me = Mme
  Mode = 2
  !
  !     To process option vector
  !
  fac = 1.E-4_DP
  !
  !     Set the nominal blow up factor used in the code.
  !
  blowup = tau
  !
  !     The nominal column scaling used in the code is
  !     the identity scaling.
  !
  D(1:N) = 1._DP
  !
  !     Define bound for number of options to change.
  !
  nopt = 1000
  !
  !     Define bound for positive value of LINK.
  !
  nlink = 100000
  ntimes = 0
  last = 1
  link = INT( Prgopt(1) )
  IF( link<=0 .OR. link>nlink ) THEN
    ERROR STOP 'DWNLSM : IN DWNNLS, THE OPTION VECTOR IS UNDEFINED'
    RETURN
  END IF
  DO
    !
    IF( link>1 ) THEN
      ntimes = ntimes + 1
      IF( ntimes>nopt ) THEN
        ERROR STOP 'DWNLSM : IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.'
        RETURN
      END IF
      !
      key = INT( Prgopt(last+1) )
      IF( key==6 .AND. Prgopt(last+2)/=0._DP ) THEN
        DO j = 1, N
          t = NORM2(W(1:m,j))
          IF( t/=0._DP ) t = 1._DP/t
          D(j) = t
        END DO
      END IF
      !
      IF( key==7 ) D(1:N) = Prgopt(last+2:last+N+1)
      IF( key==8 ) tau = MAX(drelpr,Prgopt(last+2))
      IF( key==9 ) blowup = MAX(drelpr,Prgopt(last+2))
      !
      next = INT( Prgopt(link) )
      IF( next<=0 .OR. next>nlink ) THEN
        ERROR STOP 'DWNLSM : IN DWNNLS, THE OPTION VECTOR IS UNDEFINED'
        RETURN
      END IF
      !
      last = link
      link = next
      CYCLE
    END IF
    !
    DO j = 1, N
      W(1:m,j) = W(1:m,j)*D(j)
    END DO
    !
    !     Process option vector
    !
    done = .FALSE.
    iter = 0
    itmax = 3*(N-L)
    Mode = 0
    nsoln = L
    l1 = MIN(m,L)
    !
    !     Compute scalee factor to apply to equality constraint equations.
    !
    DO j = 1, N
      Wd(j) = SUM(ABS(W(1:m,j)))
    END DO
    !
    imax = MAXLOC(Wd(1:N),1)
    eanorm = Wd(imax)
    bnorm = SUM(ABS(W(1:m,N+1)))
    alamda = eanorm/(drelpr*fac)
    !
    !     On machines, such as the VAXes using D floating, with a very
    !     limited exponent range for double precision values, the previously
    !     computed value of ALAMDA may cause an overflow condition.
    !     Therefore, this code further limits the value of ALAMDA.
    !
    alamda = MIN(alamda,SQRT(D1MACH(2)))
    !
    !     Define scaling diagonal matrix for modified Givens usage and
    !     classify equation types.
    !
    alsq = alamda**2
    DO i = 1, m
      !
      !        When equation I is heavily weighted ITYPE(I)=0,
      !        else ITYPE(I)=1.
      !
      IF( i<=me ) THEN
        t = alsq
        itemp = 0
      ELSE
        t = 1._DP
        itemp = 1
      END IF
      Scalee(i) = t
      Itype(i) = itemp
    END DO
    !
    !     Set the solution vector X(*) to zero and the column interchange
    !     matrix to the identity.
    !
    X(1:N) = 0._DP
    DO i = 1, N
      Ipivot(i) = i
    END DO
    !
    !     Perform initial triangularization in the submatrix
    !     corresponding to the unconstrained variables.
    !     Set first L components of dual vector to zero because
    !     these correspond to the unconstrained variables.
    !
    Wd(1:L) = 0._DP
    !
    !     The arrays IDOPE(*) and DOPE(*) are used to pass
    !     information to DWNLIT().  This was done to avoid
    !     a long calling sequence or the use of COMMON.
    !
    idope(1) = me
    idope(2) = nsoln
    idope(3) = l1
    !
    dope(1) = alsq
    dope(2) = eanorm
    dope(3) = tau
    CALL DWNLIT(W,Mdw,m,N,L,Ipivot,Itype,H,Scalee,Rnorm,idope,dope,done)
    me = idope(1)
    krank = idope(2)
    niv = idope(3)
    EXIT
  END DO
  !
  !     Perform WNNLS algorithm using the following steps.
  !
  !     Until(DONE)
  !        compute search direction and feasible point
  !        when (HITCON) add constraints
  !        else perform multiplier test and drop a constraint
  !        fin
  !     Compute-Final-Solution
  !
  !     To compute search direction and feasible point,
  !     solve the triangular system of currently non-active
  !     variables and store the solution in Z(*).
  !
  !     To solve system
  !     Copy right hand side into TEMP vector to use overwriting method.
  !
  DO WHILE( .NOT. (done) )
    isol = L + 1
    IF( nsoln>=isol ) THEN
      Temp(1:niv) = W(1:niv,N+1)
      DO j = nsoln, isol, -1
        IF( j>krank ) THEN
          i = niv - nsoln + j
        ELSE
          i = j
        END IF
        !
        IF( j>krank .AND. j<=L ) THEN
          Z(j) = 0._DP
        ELSE
          Z(j) = Temp(i)/W(i,j)
          CALL DAXPY(i-1,-Z(j),W(1,j),1,Temp,1)
        END IF
      END DO
    END IF
    !
    !     Increment iteration counter and check against maximum number
    !     of iterations.
    !
    iter = iter + 1
    IF( iter>itmax ) THEN
      Mode = 1
      done = .TRUE.
    END IF
    !
    !     Check to see if any constraints have become active.
    !     If so, calculate an interpolation factor so that all
    !     active constraints are removed from the basis.
    !
    alpha = 2._DP
    hitcon = .FALSE.
    DO j = L + 1, nsoln
      zz = Z(j)
      IF( zz<=0._DP ) THEN
        t = X(j)/(X(j)-zz)
        IF( t<alpha ) THEN
          alpha = t
          jcon = j
        END IF
        hitcon = .TRUE.
      END IF
    END DO
    !
    !     Compute search direction and feasible point
    !
    IF( hitcon ) THEN
      !
      !        To add constraints, use computed ALPHA to interpolate between
      !        last feasible solution X(*) and current unconstrained (and
      !        infeasible) solution Z(*).
      !
      DO j = L + 1, nsoln
        X(j) = X(j) + alpha*(Z(j)-X(j))
      END DO
      feasbl = .FALSE.
      !
      !        Remove column JCON and shift columns JCON+1 through N to the
      !        left.  Swap column JCON into the N th position.  This achieves
      !        upper Hessenberg form for the nonactive constraints and
      !        leaves an upper Hessenberg matrix to retriangularize.
      !
      20 CONTINUE
      DO i = 1, m
        t = W(i,jcon)
        W(i,jcon:N-1) = W(i,jcon+1:N)
        W(i,N) = t
      END DO
      !
      !        Update permuted index vector to reflect this shift and swap.
      !
      itemp = Ipivot(jcon)
      DO i = jcon, N - 1
        Ipivot(i) = Ipivot(i+1)
      END DO
      Ipivot(N) = itemp
      !
      !        Similarly permute X(*) vector.
      !
      X(jcon:N-1) = X(jcon+1:N)
      X(N) = 0._DP
      nsoln = nsoln - 1
      niv = niv - 1
      !
      !        Retriangularize upper Hessenberg matrix after adding
      !        constraints.
      !
      i = krank + jcon - L
      DO j = jcon, nsoln
        IF( Itype(i)==0 .AND. Itype(i+1)==0 ) THEN
          !
          !              Zero IP1 to I in column J
          !
          IF( W(i+1,j)/=0._DP ) THEN
            CALL DROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0._DP
            CALL DROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF( Itype(i)==1 .AND. Itype(i+1)==1 ) THEN
          !
          !              Zero IP1 to I in column J
          !
          IF( W(i+1,j)/=0._DP ) THEN
            CALL DROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0._DP
            CALL DROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF( Itype(i)==1 .AND. Itype(i+1)==0 ) THEN
          CALL DSWAP(N+1,W(i,1),Mdw,W(i+1,1),Mdw)
          CALL DSWAP(1,Scalee(i),1,Scalee(i+1),1)
          itemp = Itype(i+1)
          Itype(i+1) = Itype(i)
          Itype(i) = itemp
          !
          !              Swapped row was formerly a pivot element, so it will
          !              be large enough to perform elimination.
          !              Zero IP1 to I in column J.
          !
          IF( W(i+1,j)/=0._DP ) THEN
            CALL DROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0._DP
            CALL DROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        ELSEIF( Itype(i)==0 .AND. Itype(i+1)==1 ) THEN
          IF( Scalee(i)*W(i,j)**2/alsq<=(tau*eanorm)**2 ) THEN
            CALL DSWAP(N+1,W(i,1),Mdw,W(i+1,1),Mdw)
            CALL DSWAP(1,Scalee(i),1,Scalee(i+1),1)
            itemp = Itype(i+1)
            Itype(i+1) = Itype(i)
            Itype(i) = itemp
            W(i+1,j) = 0._DP
            !
            !                 Zero IP1 to I in column J
            !
          ELSEIF( W(i+1,j)/=0._DP ) THEN
            CALL DROTMG(Scalee(i),Scalee(i+1),W(i,j),W(i+1,j),sparam)
            W(i+1,j) = 0._DP
            CALL DROTM(N+1-j,W(i,j+1),Mdw,W(i+1,j+1),Mdw,sparam)
          END IF
        END IF
        i = i + 1
      END DO
      !
      !        See if the remaining coefficients in the solution set are
      !        feasible.  They should be because of the way ALPHA was
      !        determined.  If any are infeasible, it is due to roundoff
      !        error.  Any that are non-positive will be set to zero and
      !        removed from the solution set.
      !
      DO jcon = L + 1, nsoln
        IF( X(jcon)<=0._DP ) GOTO 40
      END DO
      feasbl = .TRUE.
      40  IF( .NOT. feasbl ) GOTO 20
    ELSE
      !
      !        To perform multiplier test and drop a constraint.
      !
      X(1:nsoln) = Z(1:nsoln)
      IF( nsoln<N ) X(nsoln+1:N) = 0._DP
      !
      !        Reclassify least squares equations as equalities as necessary.
      !
      i = niv + 1
      DO
        IF( i<=me ) THEN
          IF( Itype(i)==0 ) THEN
            i = i + 1
          ELSE
            CALL DSWAP(N+1,W(i,1),Mdw,W(me,1),Mdw)
            CALL DSWAP(1,Scalee(i),1,Scalee(me),1)
            itemp = Itype(i)
            Itype(i) = Itype(me)
            Itype(me) = itemp
            me = me - 1
          END IF
          CYCLE
        END IF
        !
        !        Form inner product vector WD(*) of dual coefficients.
        !
        DO j = nsoln + 1, N
          sm = 0._DP
          DO i = nsoln + 1, m
            sm = sm + Scalee(i)*W(i,j)*W(i,N+1)
          END DO
          Wd(j) = sm
        END DO
        EXIT
      END DO
      DO
        !
        !        Find J such that WD(J)=WMAX is maximum.  This determines
        !        that the incoming column J will reduce the residual vector
        !        and be positive.
        !
        wmax = 0._DP
        iwmax = nsoln + 1
        DO j = nsoln + 1, N
          IF( Wd(j)>wmax ) THEN
            wmax = Wd(j)
            iwmax = j
          END IF
        END DO
        IF( wmax<=0._DP ) GOTO 100
        !
        !        Set dual coefficients to zero for incoming column.
        !
        Wd(iwmax) = 0._DP
        !
        !        WMAX > 0.D0, so okay to move column IWMAX to solution set.
        !        Perform transformation to retriangularize, and test for near
        !        linear dependence.
        !
        !        Swap column IWMAX into NSOLN-th position to maintain upper
        !        Hessenberg form of adjacent columns, and add new column to
        !        triangular decomposition.
        !
        nsoln = nsoln + 1
        niv = niv + 1
        IF( nsoln/=iwmax ) THEN
          CALL DSWAP(m,W(1,nsoln),1,W(1,iwmax),1)
          Wd(iwmax) = Wd(nsoln)
          Wd(nsoln) = 0._DP
          itemp = Ipivot(nsoln)
          Ipivot(nsoln) = Ipivot(iwmax)
          Ipivot(iwmax) = itemp
        END IF
        !
        !        Reduce column NSOLN so that the matrix of nonactive constraints
        !        variables is triangular.
        !
        DO j = m, niv + 1, -1
          jp = j - 1
          !
          !           When operating near the ME line, test to see if the pivot
          !           element is near zero.  If so, use the largest element above
          !           it as the pivot.  This is to maintain the sharp interface
          !           between weighted and non-weighted rows in all cases.
          !
          IF( j==me+1 ) THEN
            imax = me
            amax = Scalee(me)*W(me,nsoln)**2
            DO jp = j - 1, niv, -1
              t = Scalee(jp)*W(jp,nsoln)**2
              IF( t>amax ) THEN
                imax = jp
                amax = t
              END IF
            END DO
            jp = imax
          END IF
          !
          IF( W(j,nsoln)/=0._DP ) THEN
            CALL DROTMG(Scalee(jp),Scalee(j),W(jp,nsoln),W(j,nsoln),sparam)
            W(j,nsoln) = 0._DP
            CALL DROTM(N+1-nsoln,W(jp,nsoln+1),Mdw,W(j,nsoln+1),Mdw,sparam)
          END IF
        END DO
        !
        !        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
        !        this is nonpositive or too large.  If this was true or if the
        !        pivot term was zero, reject the column as dependent.
        !
        IF( W(niv,nsoln)/=0._DP ) THEN
          isol = niv
          z2 = W(isol,N+1)/W(isol,nsoln)
          Z(nsoln) = z2
          pos = z2>0._DP
          IF( z2*eanorm>=bnorm .AND. pos )&
            pos = .NOT. (blowup*z2*eanorm>=bnorm)
          !
          !           Try to add row ME+1 as an additional equality constraint.
          !           Check size of proposed new solution component.
          !           Reject it if it is too large.
          !
        ELSEIF( niv<=me .AND. W(me+1,nsoln)/=0._DP ) THEN
          isol = me + 1
          IF( pos ) THEN
            !
            !              Swap rows ME+1 and NIV, and scale factors for these rows.
            !
            CALL DSWAP(N+1,W(me+1,1),Mdw,W(niv,1),Mdw)
            CALL DSWAP(1,Scalee(me+1),1,Scalee(niv),1)
            itemp = Itype(me+1)
            Itype(me+1) = Itype(niv)
            Itype(niv) = itemp
            me = me + 1
          END IF
        ELSE
          pos = .FALSE.
        END IF
        !
        IF( .NOT. pos ) THEN
          nsoln = nsoln - 1
          niv = niv - 1
        END IF
        IF( pos .OR. done ) EXIT
      END DO
    END IF
  END DO
  !
  !     Else perform multiplier test and drop a constraint.  To compute
  !     final solution.  Solve system, store results in X(*).
  !
  !     Copy right hand side into TEMP vector to use overwriting method.
  !
  100  isol = 1
  IF( nsoln>=isol ) THEN
    Temp(1:niv) = W(1:niv,N+1)
    DO j = nsoln, isol, -1
      IF( j>krank ) THEN
        i = niv - nsoln + j
      ELSE
        i = j
      END IF
      !
      IF( j>krank .AND. j<=L ) THEN
        Z(j) = 0._DP
      ELSE
        Z(j) = Temp(i)/W(i,j)
        CALL DAXPY(i-1,-Z(j),W(1,j),1,Temp,1)
      END IF
    END DO
  END IF
  !
  !     Solve system.
  !
  X(1:nsoln) = Z(1:nsoln)
  !
  !     Apply Householder transformations to X(*) if KRANK<L
  !
  IF( krank<L ) THEN
    DO i = 1, krank
      CALL DH12(2,i,krank+1,L,W(i,1),Mdw,H(i),X,1,1,1)
    END DO
  END IF
  !
  !     Fill in trailing zeroes for constrained variables not in solution.
  !
  IF( nsoln<N ) X(nsoln+1:N) = 0._DP
  !
  !     Permute solution vector to natural order.
  !
  DO i = 1, N
    j = i
    DO WHILE( Ipivot(j)/=i )
      j = j + 1
    END DO
    !
    Ipivot(j) = Ipivot(i)
    Ipivot(i) = j
    CALL DSWAP(1,X(j),1,X(i),1)
  END DO
  !
  !     Rescale the solution using the column scaling.
  !
  DO j = 1, N
    X(j) = X(j)*D(j)
  END DO
  !
  DO i = nsoln + 1, m
    t = W(i,N+1)
    IF( i<=me ) t = t/alamda
    t = (Scalee(i)*t)*t
    Rnorm = Rnorm + t
  END DO
  !
  Rnorm = SQRT(Rnorm)
  !
END SUBROUTINE DWNLSM