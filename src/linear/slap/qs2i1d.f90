!** QS2I1D
SUBROUTINE QS2I1D(Ia,Ja,A,N,Kflag)
  !>
  !  Sort an integer array, moving an integer and DP array.
  !            This routine sorts the integer array IA and makes the same
  !            interchanges in the integer array JA and the double pre-
  !            cision array A.  The array IA may be sorted in increasing
  !            order or decreasing order.  A slightly modified QUICKSORT
  !            algorithm is used.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  N6A2A
  !***
  ! **Type:**      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
  !***
  ! **Keywords:**  SINGLETON QUICKSORT, SLAP, SORT, SORTING
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Kahaner, D. K., (NBS)
  !           Seager, M. K., (LLNL) seager@llnl.gov
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !     Written by Rondall E Jones
  !     Modified by John A. Wisniewski to use the Singleton QUICKSORT
  !     algorithm. date 18 November 1976.
  !
  !     Further modified by David K. Kahaner
  !     National Bureau of Standards
  !     August, 1981
  !
  !     Even further modification made to bring the code up to the
  !     Fortran 77 level and make it more readable and to carry
  !     along one integer array and one double precision array during
  !     the sort by
  !     Mark K. Seager
  !     Lawrence Livermore National Laboratory
  !     November, 1987
  !     This routine was adapted from the ISORT routine.
  !
  !     ABSTRACT
  !         This routine sorts an integer array IA and makes the same
  !         interchanges in the integer array JA and the double precision
  !         array A.
  !         The array IA may be sorted in increasing order or decreasing
  !         order.  A slightly modified quicksort algorithm is used.
  !
  !     DESCRIPTION OF PARAMETERS
  !        IA - Integer array of values to be sorted.
  !        JA - Integer array to be carried along.
  !         A - Double Precision array to be carried along.
  !         N - Number of values in integer array IA to be sorted.
  !     KFLAG - Control parameter
  !           = 1 means sort IA in INCREASING order.
  !           =-1 means sort IA in DECREASING order.
  !
  !***
  ! **See also:**  DS2Y
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An Efficient Algorithm
  !                 for Sorting With Minimal Storage, Communications ACM
  !                 12:3 (1969), pp.185-7.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761118  DATE WRITTEN
  !   890125  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !   920929  Corrected format of reference.  (FNF)
  !   921012  Corrected all f.p. constants to double precision.  (FNF)
  USE service, ONLY : XERMSG
  !     .. Scalar Arguments ..
  INTEGER Kflag, N
  !     .. Array Arguments ..
  REAL(8) :: A(N)
  INTEGER Ia(N), Ja(N)
  !     .. Local Scalars ..
  REAL(8) :: r, ta, tta
  INTEGER i, iit, ij, it, j, jjt, jt, k, kk, l, m, nn
  !     .. Local Arrays ..
  INTEGER il(21), iu(21)
  !* FIRST EXECUTABLE STATEMENT  QS2I1D
  nn = N
  IF ( nn<1 ) THEN
    CALL XERMSG('SLATEC','QS2I1D',&
      'The number of values to be sorted was not positive.',1,1)
    RETURN
  END IF
  IF ( N==1 ) RETURN
  kk = ABS(Kflag)
  IF ( kk/=1 ) THEN
    CALL XERMSG('SLATEC','QS2I1D',&
      'The sort control parameter, K, was not 1 or -1.',2,1)
    RETURN
  END IF
  !
  !     Alter array IA to get decreasing order if needed.
  !
  IF ( Kflag<1 ) THEN
    DO i = 1, nn
      Ia(i) = -Ia(i)
    END DO
  END IF
  !
  !     Sort IA and carry JA and A along.
  !     And now...Just a little black magic...
  m = 1
  i = 1
  j = nn
  r = .375D0
  100 CONTINUE
  IF ( r<=0.5898437D0 ) THEN
    r = r + 3.90625D-2
  ELSE
    r = r - .21875D0
  END IF
  200  k = i
  !
  !     Select a central element of the array and save it in location
  !     it, jt, at.
  !
  ij = i + INT((j-i)*r)
  it = Ia(ij)
  jt = Ja(ij)
  ta = A(ij)
  !
  !     If first element of array is greater than it, interchange with it.
  !
  IF ( Ia(i)>it ) THEN
    Ia(ij) = Ia(i)
    Ia(i) = it
    it = Ia(ij)
    Ja(ij) = Ja(i)
    Ja(i) = jt
    jt = Ja(ij)
    A(ij) = A(i)
    A(i) = ta
    ta = A(ij)
  END IF
  l = j
  !
  !     If last element of array is less than it, swap with it.
  !
  IF ( Ia(j)<it ) THEN
    Ia(ij) = Ia(j)
    Ia(j) = it
    it = Ia(ij)
    Ja(ij) = Ja(j)
    Ja(j) = jt
    jt = Ja(ij)
    A(ij) = A(j)
    A(j) = ta
    ta = A(ij)
    !
    !     If first element of array is greater than it, swap with it.
    !
    IF ( Ia(i)>it ) THEN
      Ia(ij) = Ia(i)
      Ia(i) = it
      it = Ia(ij)
      Ja(ij) = Ja(i)
      Ja(i) = jt
      jt = Ja(ij)
      A(ij) = A(i)
      A(i) = ta
      ta = A(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is
    !     smaller than it.
    !
    l = l - 1
    IF ( Ia(l)<=it ) THEN
      DO
        !
        !     Find an element in the first half of the array which is
        !     greater than it.
        !
        k = k + 1
        IF ( Ia(k)>=it ) THEN
          !
          !     Interchange these elements.
          !
          IF ( k<=l ) THEN
            iit = Ia(l)
            Ia(l) = Ia(k)
            Ia(k) = iit
            jjt = Ja(l)
            Ja(l) = Ja(k)
            Ja(k) = jjt
            tta = A(l)
            A(l) = A(k)
            A(k) = tta
            EXIT
          END IF
          !
          !     Save upper and lower subscripts of the array yet to be sorted.
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
  !     Begin again on another portion of the unsorted array.
  !
  300  m = m - 1
  IF ( m==0 ) THEN
    !
    !     Clean up, if necessary.
    !
    IF ( Kflag<1 ) THEN
      DO i = 1, nn
        Ia(i) = -Ia(i)
      END DO
    END IF
    RETURN
  ELSE
    i = il(m)
    j = iu(m)
  END IF
  400 CONTINUE
  IF ( j-i>=1 ) GOTO 200
  IF ( i==j ) GOTO 300
  IF ( i==1 ) GOTO 100
  i = i - 1
  DO
    i = i + 1
    IF ( i==j ) GOTO 300
    it = Ia(i+1)
    jt = Ja(i+1)
    ta = A(i+1)
    IF ( Ia(i)>it ) THEN
      k = i
      DO
        Ia(k+1) = Ia(k)
        Ja(k+1) = Ja(k)
        A(k+1) = A(k)
        k = k - 1
        IF ( it>=Ia(k) ) THEN
          Ia(k+1) = it
          Ja(k+1) = jt
          A(k+1) = ta
          EXIT
        END IF
      END DO
    END IF
  END DO
  !------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
  RETURN
END SUBROUTINE QS2I1D
