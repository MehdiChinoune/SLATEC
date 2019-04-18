!** SSORT
SUBROUTINE SSORT(X,Y,N,Kflag)
  !>
  !***
  !  Sort an array and optionally make the same interchanges in
  !            an auxiliary array.  The array may be sorted in increasing
  !            or decreasing order.  A slightly modified QUICKSORT
  !            algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A2B
  !***
  ! **Type:**      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   SSORT sorts array X and optionally makes the same interchanges in
  !   array Y.  The array X may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      X - array of values to be sorted   (usually abscissas)
  !      Y - array to be (optionally) carried along
  !      N - number of values in array X to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort X in increasing order and carry Y along.
  !            =  1  means sort X in increasing order (ignoring Y)
  !            = -1  means sort X in decreasing order (ignoring Y)
  !            = -2  means sort X in decreasing order and carry Y along.
  !
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An efficient algorithm
  !                 for sorting with minimal storage, Communications of
  !                 the ACM, 12, 3 (1969), pp. 185-187.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761101  DATE WRITTEN
  !   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced statement labels.  (WRB)
  !   891024  Changed category.  (WRB)
  !   891024  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
  !   920501  Reformatted the REFERENCES section.  (DWL, WRB)
  !   920519  Clarified error messages.  (DWL)
  !   920801  Declarations section rebuilt and code restructured to use
  !           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
  USE service, ONLY : XERMSG
  !     .. Scalar Arguments ..
  INTEGER Kflag, N
  !     .. Array Arguments ..
  REAL X(*), Y(*)
  !     .. Local Scalars ..
  REAL r, t, tt, tty, ty
  INTEGER i, ij, j, k, kk, l, m, nn
  !     .. Local Arrays ..
  INTEGER il(21), iu(21)
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
  !* FIRST EXECUTABLE STATEMENT  SSORT
  nn = N
  IF ( nn<1 ) THEN
    CALL XERMSG('SLATEC','SSORT',&
      'The number of values to be sorted is not positive.',1,1)
    RETURN
  END IF
  !
  kk = ABS(Kflag)
  IF ( kk/=1.AND.kk/=2 ) THEN
    CALL XERMSG('SLATEC','SSORT',&
      'The sort control parameter, K, is not 2, 1, -1, or -2.',2,1)
    RETURN
  END IF
  !
  !     Alter array X to get decreasing order if needed
  !
  IF ( Kflag<=-1 ) THEN
    DO i = 1, nn
      X(i) = -X(i)
    END DO
  END IF
  !
  IF ( kk==2 ) THEN
    !
    !     Sort X and carry Y along
    !
    m = 1
    i = 1
    j = nn
    r = 0.375E0
    GOTO 500
  ELSE
    !
    !     Sort X only
    !
    m = 1
    i = 1
    j = nn
    r = 0.375E0
  END IF
  !
  100 CONTINUE
  IF ( i==j ) GOTO 300
  IF ( r<=0.5898437E0 ) THEN
    r = r + 3.90625E-2
  ELSE
    r = r - 0.21875E0
  END IF
  !
  200  k = i
  !
  !     Select a central element of the array and save it in location T
  !
  ij = i + INT((j-i)*r)
  t = X(ij)
  !
  !     If first element of array is greater than T, interchange with T
  !
  IF ( X(i)>t ) THEN
    X(ij) = X(i)
    X(i) = t
    t = X(ij)
  END IF
  l = j
  !
  !     If last element of array is less than than T, interchange with T
  !
  IF ( X(j)<t ) THEN
    X(ij) = X(j)
    X(j) = t
    t = X(ij)
    !
    !        If first element of array is greater than T, interchange with T
    !
    IF ( X(i)>t ) THEN
      X(ij) = X(i)
      X(i) = t
      t = X(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
    l = l - 1
    IF ( X(l)<=t ) THEN
      DO
        !
        !     Find an element in the first half of the array which is greater
        !     than T
        !
        k = k + 1
        IF ( X(k)>=t ) THEN
          !
          !     Interchange these elements
          !
          IF ( k<=l ) THEN
            tt = X(l)
            X(l) = X(k)
            X(k) = tt
            EXIT
          END IF
          !
          !     Save upper and lower subscripts of the array yet to be sorted
          !
          IF ( l-i>j-k ) THEN
            il(m) = i
            iu(m) = l
            i = k
            m = m + 1
          ELSE
            il(m) = k
            iu(m) = j
            j = l
            m = m + 1
          END IF
          GOTO 400
        END IF
      END DO
    END IF
  END DO
  !
  !     Begin again on another portion of the unsorted array
  !
  300  m = m - 1
  IF ( m==0 ) GOTO 900
  i = il(m)
  j = iu(m)
  !
  400 CONTINUE
  IF ( j-i>=1 ) GOTO 200
  IF ( i==1 ) GOTO 100
  i = i - 1
  DO
    !
    i = i + 1
    IF ( i==j ) GOTO 300
    t = X(i+1)
    IF ( X(i)>t ) THEN
      k = i
      DO
        !
        X(k+1) = X(k)
        k = k - 1
        IF ( t>=X(k) ) THEN
          X(k+1) = t
          EXIT
        END IF
      END DO
    END IF
  END DO
  !
  500 CONTINUE
  IF ( i==j ) GOTO 700
  IF ( r<=0.5898437E0 ) THEN
    r = r + 3.90625E-2
  ELSE
    r = r - 0.21875E0
  END IF
  !
  600  k = i
  !
  !     Select a central element of the array and save it in location T
  !
  ij = i + INT((j-i)*r)
  t = X(ij)
  ty = Y(ij)
  !
  !     If first element of array is greater than T, interchange with T
  !
  IF ( X(i)>t ) THEN
    X(ij) = X(i)
    X(i) = t
    t = X(ij)
    Y(ij) = Y(i)
    Y(i) = ty
    ty = Y(ij)
  END IF
  l = j
  !
  !     If last element of array is less than T, interchange with T
  !
  IF ( X(j)<t ) THEN
    X(ij) = X(j)
    X(j) = t
    t = X(ij)
    Y(ij) = Y(j)
    Y(j) = ty
    ty = Y(ij)
    !
    !        If first element of array is greater than T, interchange with T
    !
    IF ( X(i)>t ) THEN
      X(ij) = X(i)
      X(i) = t
      t = X(ij)
      Y(ij) = Y(i)
      Y(i) = ty
      ty = Y(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
    l = l - 1
    IF ( X(l)<=t ) THEN
      DO
        !
        !     Find an element in the first half of the array which is greater
        !     than T
        !
        k = k + 1
        IF ( X(k)>=t ) THEN
          !
          !     Interchange these elements
          !
          IF ( k<=l ) THEN
            tt = X(l)
            X(l) = X(k)
            X(k) = tt
            tty = Y(l)
            Y(l) = Y(k)
            Y(k) = tty
            EXIT
          END IF
          !
          !     Save upper and lower subscripts of the array yet to be sorted
          !
          IF ( l-i>j-k ) THEN
            il(m) = i
            iu(m) = l
            i = k
            m = m + 1
          ELSE
            il(m) = k
            iu(m) = j
            j = l
            m = m + 1
          END IF
          GOTO 800
        END IF
      END DO
    END IF
  END DO
  !
  !     Begin again on another portion of the unsorted array
  !
  700  m = m - 1
  IF ( m==0 ) GOTO 900
  i = il(m)
  j = iu(m)
  !
  800 CONTINUE
  IF ( j-i>=1 ) GOTO 600
  IF ( i==1 ) GOTO 500
  i = i - 1
  DO
    !
    i = i + 1
    IF ( i==j ) GOTO 700
    t = X(i+1)
    ty = Y(i+1)
    IF ( X(i)>t ) THEN
      k = i
      DO
        !
        X(k+1) = X(k)
        Y(k+1) = Y(k)
        k = k - 1
        IF ( t>=X(k) ) THEN
          X(k+1) = t
          Y(k+1) = ty
          EXIT
        END IF
      END DO
    END IF
  END DO
  !
  !     Clean up
  !
  900 CONTINUE
  IF ( Kflag<=-1 ) THEN
    DO i = 1, nn
      X(i) = -X(i)
    END DO
  END IF
END SUBROUTINE SSORT
