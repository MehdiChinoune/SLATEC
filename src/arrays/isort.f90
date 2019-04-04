!** ISORT
SUBROUTINE ISORT(Ix,Iy,N,Kflag)
  IMPLICIT NONE
  !>
  !***
  !  Sort an array and optionally make the same interchanges in
  !            an auxiliary array.  The array may be sorted in increasing
  !            or decreasing order.  A slightly modified QUICKSORT
  !            algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A2A
  !***
  ! **Type:**      INTEGER (SSORT-S, DSORT-D, ISORT-I)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Kahaner, D. K., (NBS)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   ISORT sorts array IX and optionally makes the same interchanges in
  !   array IY.  The array IX may be sorted in increasing order or
  !   decreasing order.  A slightly modified quicksort algorithm is used.
  !
  !   Description of Parameters
  !      IX - integer array of values to be sorted
  !      IY - integer array to be (optionally) carried along
  !      N  - number of values in integer array IX to be sorted
  !      KFLAG - control parameter
  !            =  2  means sort IX in increasing order and carry IY along.
  !            =  1  means sort IX in increasing order (ignoring IY)
  !            = -1  means sort IX in decreasing order (ignoring IY)
  !            = -2  means sort IX in decreasing order and carry IY along.
  !
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An efficient algorithm
  !                 for sorting with minimal storage, Communications of
  !                 the ACM, 12, 3 (1969), pp. 185-187.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761118  DATE WRITTEN
  !   810801  Modified by David K. Kahaner.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced statement labels.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
  !   920501  Reformatted the REFERENCES section.  (DWL, WRB)
  !   920519  Clarified error messages.  (DWL)
  !   920801  Declarations section rebuilt and code restructured to use
  !           IF-THEN-ELSE-ENDIF.  (RWC, WRB)

  !     .. Scalar Arguments ..
  INTEGER Kflag, N
  !     .. Array Arguments ..
  INTEGER Ix(*), Iy(*)
  !     .. Local Scalars ..
  REAL r
  INTEGER i, ij, j, k, kk, l, m, nn, t, tt, tty, ty
  !     .. Local Arrays ..
  INTEGER il(21), iu(21)
  !     .. External Subroutines ..
  EXTERNAL :: XERMSG
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
  !* FIRST EXECUTABLE STATEMENT  ISORT
  nn = N
  IF ( nn<1 ) THEN
    CALL XERMSG('SLATEC','ISORT',&
      'The number of values to be sorted is not positive.',1,1)
    RETURN
  END IF
  !
  kk = ABS(Kflag)
  IF ( kk/=1.AND.kk/=2 ) THEN
    CALL XERMSG('SLATEC','ISORT',&
      'The sort control parameter, K, is not 2, 1, -1, or -2.',2,1)
    RETURN
  END IF
  !
  !     Alter array IX to get decreasing order if needed
  !
  IF ( Kflag<=-1 ) THEN
    DO i = 1, nn
      Ix(i) = -Ix(i)
    END DO
  END IF
  !
  IF ( kk==2 ) THEN
    !
    !     Sort IX and carry IY along
    !
    m = 1
    i = 1
    j = nn
    r = 0.375E0
    GOTO 500
  ELSE
    !
    !     Sort IX only
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
  t = Ix(ij)
  !
  !     If first element of array is greater than T, interchange with T
  !
  IF ( Ix(i)>t ) THEN
    Ix(ij) = Ix(i)
    Ix(i) = t
    t = Ix(ij)
  END IF
  l = j
  !
  !     If last element of array is less than than T, interchange with T
  !
  IF ( Ix(j)<t ) THEN
    Ix(ij) = Ix(j)
    Ix(j) = t
    t = Ix(ij)
    !
    !        If first element of array is greater than T, interchange with T
    !
    IF ( Ix(i)>t ) THEN
      Ix(ij) = Ix(i)
      Ix(i) = t
      t = Ix(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
    l = l - 1
    IF ( Ix(l)<=t ) THEN
      DO
        !
        !     Find an element in the first half of the array which is greater
        !     than T
        !
        k = k + 1
        IF ( Ix(k)>=t ) THEN
          !
          !     Interchange these elements
          !
          IF ( k<=l ) THEN
            tt = Ix(l)
            Ix(l) = Ix(k)
            Ix(k) = tt
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
    t = Ix(i+1)
    IF ( Ix(i)>t ) THEN
      k = i
      DO
        !
        Ix(k+1) = Ix(k)
        k = k - 1
        IF ( t>=Ix(k) ) THEN
          Ix(k+1) = t
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
  t = Ix(ij)
  ty = Iy(ij)
  !
  !     If first element of array is greater than T, interchange with T
  !
  IF ( Ix(i)>t ) THEN
    Ix(ij) = Ix(i)
    Ix(i) = t
    t = Ix(ij)
    Iy(ij) = Iy(i)
    Iy(i) = ty
    ty = Iy(ij)
  END IF
  l = j
  !
  !     If last element of array is less than T, interchange with T
  !
  IF ( Ix(j)<t ) THEN
    Ix(ij) = Ix(j)
    Ix(j) = t
    t = Ix(ij)
    Iy(ij) = Iy(j)
    Iy(j) = ty
    ty = Iy(ij)
    !
    !        If first element of array is greater than T, interchange with T
    !
    IF ( Ix(i)>t ) THEN
      Ix(ij) = Ix(i)
      Ix(i) = t
      t = Ix(ij)
      Iy(ij) = Iy(i)
      Iy(i) = ty
      ty = Iy(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
    l = l - 1
    IF ( Ix(l)<=t ) THEN
      DO
        !
        !     Find an element in the first half of the array which is greater
        !     than T
        !
        k = k + 1
        IF ( Ix(k)>=t ) THEN
          !
          !     Interchange these elements
          !
          IF ( k<=l ) THEN
            tt = Ix(l)
            Ix(l) = Ix(k)
            Ix(k) = tt
            tty = Iy(l)
            Iy(l) = Iy(k)
            Iy(k) = tty
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
    t = Ix(i+1)
    ty = Iy(i+1)
    IF ( Ix(i)>t ) THEN
      k = i
      DO
        !
        Ix(k+1) = Ix(k)
        Iy(k+1) = Iy(k)
        k = k - 1
        IF ( t>=Ix(k) ) THEN
          Ix(k+1) = t
          Iy(k+1) = ty
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
      Ix(i) = -Ix(i)
    END DO
  END IF
END SUBROUTINE ISORT
