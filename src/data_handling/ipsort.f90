!** IPSORT
SUBROUTINE IPSORT(Ix,N,Iperm,Kflag,Ier)
  !>
  !  Return the permutation vector generated by sorting a given
  !            array and, optionally, rearrange the elements of the array.
  !            The array may be sorted in increasing or decreasing order.
  !            A slightly modified quicksort algorithm is used.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N6A1A, N6A2A
  !***
  ! **Type:**      INTEGER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
  !***
  ! **Keywords:**  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !           Kahaner, D. K., (NBS)
  !           Rhoads, G. S., (NBS)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !   IPSORT returns the permutation vector IPERM generated by sorting
  !   the array IX and, optionally, rearranges the values in IX.  IX may
  !   be sorted in increasing or decreasing order.  A slightly modified
  !   quicksort algorithm is used.
  !
  !   IPERM is such that IX(IPERM(I)) is the Ith value in the
  !   rearrangement of IX.  IPERM may be applied to another array by
  !   calling IPPERM, SPPERM, DPPERM or HPPERM.
  !
  !   The main difference between IPSORT and its active sorting equivalent
  !   ISORT is that the data are referenced indirectly rather than
  !   directly.  Therefore, IPSORT should require approximately twice as
  !   long to execute as ISORT.  However, IPSORT is more general.
  !
  !   Description of Parameters
  !      IX - input/output -- integer array of values to be sorted.
  !           If ABS(KFLAG) = 2, then the values in IX will be
  !           rearranged on output; otherwise, they are unchanged.
  !      N  - input -- number of values in array IX to be sorted.
  !      IPERM - output -- permutation array such that IPERM(I) is the
  !              index of the value in the original order of the
  !              IX array that is in the Ith location in the sorted
  !              order.
  !      KFLAG - input -- control parameter:
  !            =  2  means return the permutation vector resulting from
  !                  sorting IX in increasing order and sort IX also.
  !            =  1  means return the permutation vector resulting from
  !                  sorting IX in increasing order and do not sort IX.
  !            = -1  means return the permutation vector resulting from
  !                  sorting IX in decreasing order and do not sort IX.
  !            = -2  means return the permutation vector resulting from
  !                  sorting IX in decreasing order and sort IX also.
  !      IER - output -- error indicator:
  !          =  0  if no error,
  !          =  1  if N is zero or negative,
  !          =  2  if KFLAG is not 2, 1, -1, or -2.
  !***
  ! **References:**  R. C. Singleton, Algorithm 347, An efficient algorithm
  !                 for sorting with minimal storage, Communications of
  !                 the ACM, 12, 3 (1969), pp. 185-187.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   761101  DATE WRITTEN
  !   761118  Modified by John A. Wisniewski to use the Singleton
  !           quicksort algorithm.
  !   810801  Further modified by David K. Kahaner.
  !   870423  Modified by Gregory S. Rhoads for passive sorting with the
  !           option for the rearrangement of the original data.
  !   890620  Algorithm for rearranging the data vector corrected by R.
  !           Boisvert.
  !   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
  !   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
  !   920507  Modified by M. McClain to revise prologue text.
  !   920818  Declarations section rebuilt and code restructured to use
  !           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
  USE service, ONLY : XERMSG
  !     .. Scalar Arguments ..
  INTEGER Ier, Kflag, N
  !     .. Array Arguments ..
  INTEGER Iperm(N), Ix(N)
  !     .. Local Scalars ..
  REAL r
  INTEGER i, ij, indx, indx0, istrt, itemp, j, k, kk, l, lm, lmt, m, nn
  !     .. Local Arrays ..
  INTEGER il(21), iu(21)  !     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
  !* FIRST EXECUTABLE STATEMENT  IPSORT
  Ier = 0
  nn = N
  IF ( nn<1 ) THEN
    Ier = 1
    CALL XERMSG('SLATEC','IPSORT',&
      'The number of values to be sorted, N, is not positive.',Ier,1)
    RETURN
  END IF
  kk = ABS(Kflag)
  IF ( kk/=1.AND.kk/=2 ) THEN
    Ier = 2
    CALL XERMSG('SLATEC','IPSORT',&
      'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.',Ier,1)
    RETURN
  END IF
  !
  !     Initialize permutation vector
  !
  DO i = 1, nn
    Iperm(i) = i
  END DO
  !
  !     Return if only one value is to be sorted
  !
  IF ( nn==1 ) RETURN
  !
  !     Alter array IX to get decreasing order if needed
  !
  IF ( Kflag<=-1 ) THEN
    DO i = 1, nn
      Ix(i) = -Ix(i)
    END DO
  END IF
  !
  !     Sort IX only
  !
  m = 1
  i = 1
  j = nn
  r = .375E0
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
  !     Select a central element of the array and save it in location L
  !
  ij = i + INT((j-i)*r)
  lm = Iperm(ij)
  !
  !     If first element of array is greater than LM, interchange with LM
  !
  IF ( Ix(Iperm(i))>Ix(lm) ) THEN
    Iperm(ij) = Iperm(i)
    Iperm(i) = lm
    lm = Iperm(ij)
  END IF
  l = j
  !
  !     If last element of array is less than LM, interchange with LM
  !
  IF ( Ix(Iperm(j))<Ix(lm) ) THEN
    Iperm(ij) = Iperm(j)
    Iperm(j) = lm
    lm = Iperm(ij)
    !
    !        If first element of array is greater than LM, interchange
    !        with LM
    !
    IF ( Ix(Iperm(i))>Ix(lm) ) THEN
      Iperm(ij) = Iperm(i)
      Iperm(i) = lm
      lm = Iperm(ij)
    END IF
  END IF
  DO
    !
    !     Find an element in the second half of the array which is smaller
    !     than LM
    !
    l = l - 1
    IF ( Ix(Iperm(l))<=Ix(lm) ) THEN
      DO
        !
        !     Find an element in the first half of the array which is greater
        !     than LM
        !
        k = k + 1
        IF ( Ix(Iperm(k))>=Ix(lm) ) THEN
          !
          !     Interchange these elements
          !
          IF ( k<=l ) THEN
            lmt = Iperm(l)
            Iperm(l) = Iperm(k)
            Iperm(k) = lmt
            EXIT
          ELSE
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
        END IF
      END DO
    END IF
  END DO
  !
  !     Begin again on another portion of the unsorted array
  !
  300  m = m - 1
  IF ( m==0 ) THEN
    !
    !     Clean up
    !
    IF ( Kflag<=-1 ) THEN
      DO i = 1, nn
        Ix(i) = -Ix(i)
      END DO
    END IF
    !
    !     Rearrange the values of IX if desired
    !
    IF ( kk==2 ) THEN
      !
      !        Use the IPERM vector as a flag.
      !        If IPERM(I) < 0, then the I-th value is in correct location
      !
      DO istrt = 1, nn
        IF ( Iperm(istrt)>=0 ) THEN
          indx = istrt
          indx0 = indx
          itemp = Ix(istrt)
          DO
            IF ( Iperm(indx)>0 ) THEN
              Ix(indx) = Ix(Iperm(indx))
              indx0 = indx
              Iperm(indx) = -Iperm(indx)
              indx = ABS(Iperm(indx))
              CYCLE
            END IF
            Ix(indx0) = itemp
            EXIT
          END DO
        END IF
      END DO
      !
      !        Revert the signs of the IPERM values
      !
      DO i = 1, nn
        Iperm(i) = -Iperm(i)
      END DO
      !
    END IF
    RETURN
  ELSE
    i = il(m)
    j = iu(m)
  END IF
  !
  400 CONTINUE
  IF ( j-i>=1 ) GOTO 200
  IF ( i==1 ) GOTO 100
  i = i - 1
  DO
    !
    i = i + 1
    IF ( i==j ) GOTO 300
    lm = Iperm(i+1)
    IF ( Ix(Iperm(i))>Ix(lm) ) THEN
      k = i
      DO
        !
        Iperm(k+1) = Iperm(k)
        k = k - 1
        !
        IF ( Ix(lm)>=Ix(Iperm(k)) ) THEN
          Iperm(k+1) = lm
          EXIT
        END IF
      END DO
    END IF
  END DO
  !
  RETURN
END SUBROUTINE IPSORT
