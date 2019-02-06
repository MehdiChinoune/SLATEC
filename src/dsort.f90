!*==DSORT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DSORT
      SUBROUTINE DSORT(Dx,Dy,N,Kflag)
      IMPLICIT NONE
!*--DSORT5
!***BEGIN PROLOGUE  DSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   DSORT sorts array DX and optionally makes the same interchanges in
!   array DY.  The array DX may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      DX - array of values to be sorted   (usually abscissas)
!      DY - array to be (optionally) carried along
!      N  - number of values in array DX to be sorted
!      KFLAG - control parameter
!            =  2  means sort DX in increasing order and carry DY along.
!            =  1  means sort DX in increasing order (ignoring DY)
!            = -1  means sort DX in decreasing order (ignoring DY)
!            = -2  means sort DX in decreasing order and carry DY along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to DX,DY; changed
!           code to parallel SSORT. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  DSORT
!     .. Scalar Arguments ..
      INTEGER Kflag , N
!     .. Array Arguments ..
      DOUBLE PRECISION Dx(*) , Dy(*)
!     .. Local Scalars ..
      DOUBLE PRECISION r , t , tt , tty , ty
      INTEGER i , ij , j , k , kk , l , m , nn
!     .. Local Arrays ..
      INTEGER il(21) , iu(21)
!     .. External Subroutines ..
      EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT
!***FIRST EXECUTABLE STATEMENT  DSORT
      nn = N
      IF ( nn<1 ) THEN
        CALL XERMSG('SLATEC','DSORT',
     &              'The number of values to be sorted is not positive.',1,1)
        RETURN
      ENDIF
!
      kk = ABS(Kflag)
      IF ( kk/=1.AND.kk/=2 ) THEN
        CALL XERMSG('SLATEC','DSORT',
     &              'The sort control parameter, K, is not 2, 1, -1, or -2.',2,
     &              1)
        RETURN
      ENDIF
!
!     Alter array DX to get decreasing order if needed
!
      IF ( Kflag<=-1 ) THEN
        DO i = 1 , nn
          Dx(i) = -Dx(i)
        ENDDO
      ENDIF
!
      IF ( kk==2 ) THEN
!
!     Sort DX and carry DY along
!
        m = 1
        i = 1
        j = nn
        r = 0.375D0
        GOTO 500
      ELSE
!
!     Sort DX only
!
        m = 1
        i = 1
        j = nn
        r = 0.375D0
      ENDIF
!
 100  IF ( i==j ) GOTO 300
      IF ( r<=0.5898437D0 ) THEN
        r = r + 3.90625D-2
      ELSE
        r = r - 0.21875D0
      ENDIF
!
 200  k = i
!
!     Select a central element of the array and save it in location T
!
      ij = i + INT((j-i)*r)
      t = Dx(ij)
!
!     If first element of array is greater than T, interchange with T
!
      IF ( Dx(i)>t ) THEN
        Dx(ij) = Dx(i)
        Dx(i) = t
        t = Dx(ij)
      ENDIF
      l = j
!
!     If last element of array is less than than T, interchange with T
!
      IF ( Dx(j)<t ) THEN
        Dx(ij) = Dx(j)
        Dx(j) = t
        t = Dx(ij)
!
!        If first element of array is greater than T, interchange with T
!
        IF ( Dx(i)>t ) THEN
          Dx(ij) = Dx(i)
          Dx(i) = t
          t = Dx(ij)
        ENDIF
      ENDIF
      DO
!
!     Find an element in the second half of the array which is smaller
!     than T
!
        l = l - 1
        IF ( Dx(l)<=t ) THEN
          DO
!
!     Find an element in the first half of the array which is greater
!     than T
!
            k = k + 1
            IF ( Dx(k)>=t ) THEN
!
!     Interchange these elements
!
              IF ( k<=l ) THEN
                tt = Dx(l)
                Dx(l) = Dx(k)
                Dx(k) = tt
                EXIT
              ENDIF
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
              ENDIF
              GOTO 400
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!     Begin again on another portion of the unsorted array
!
 300  m = m - 1
      IF ( m==0 ) GOTO 900
      i = il(m)
      j = iu(m)
!
 400  IF ( j-i>=1 ) GOTO 200
      IF ( i==1 ) GOTO 100
      i = i - 1
      DO
!
        i = i + 1
        IF ( i==j ) GOTO 300
        t = Dx(i+1)
        IF ( Dx(i)>t ) THEN
          k = i
          DO
!
            Dx(k+1) = Dx(k)
            k = k - 1
            IF ( t>=Dx(k) ) THEN
              Dx(k+1) = t
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
 500  IF ( i==j ) GOTO 700
      IF ( r<=0.5898437D0 ) THEN
        r = r + 3.90625D-2
      ELSE
        r = r - 0.21875D0
      ENDIF
!
 600  k = i
!
!     Select a central element of the array and save it in location T
!
      ij = i + INT((j-i)*r)
      t = Dx(ij)
      ty = Dy(ij)
!
!     If first element of array is greater than T, interchange with T
!
      IF ( Dx(i)>t ) THEN
        Dx(ij) = Dx(i)
        Dx(i) = t
        t = Dx(ij)
        Dy(ij) = Dy(i)
        Dy(i) = ty
        ty = Dy(ij)
      ENDIF
      l = j
!
!     If last element of array is less than T, interchange with T
!
      IF ( Dx(j)<t ) THEN
        Dx(ij) = Dx(j)
        Dx(j) = t
        t = Dx(ij)
        Dy(ij) = Dy(j)
        Dy(j) = ty
        ty = Dy(ij)
!
!        If first element of array is greater than T, interchange with T
!
        IF ( Dx(i)>t ) THEN
          Dx(ij) = Dx(i)
          Dx(i) = t
          t = Dx(ij)
          Dy(ij) = Dy(i)
          Dy(i) = ty
          ty = Dy(ij)
        ENDIF
      ENDIF
      DO
!
!     Find an element in the second half of the array which is smaller
!     than T
!
        l = l - 1
        IF ( Dx(l)<=t ) THEN
          DO
!
!     Find an element in the first half of the array which is greater
!     than T
!
            k = k + 1
            IF ( Dx(k)>=t ) THEN
!
!     Interchange these elements
!
              IF ( k<=l ) THEN
                tt = Dx(l)
                Dx(l) = Dx(k)
                Dx(k) = tt
                tty = Dy(l)
                Dy(l) = Dy(k)
                Dy(k) = tty
                EXIT
              ENDIF
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
              ENDIF
              GOTO 800
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!     Begin again on another portion of the unsorted array
!
 700  m = m - 1
      IF ( m==0 ) GOTO 900
      i = il(m)
      j = iu(m)
!
 800  IF ( j-i>=1 ) GOTO 600
      IF ( i==1 ) GOTO 500
      i = i - 1
      DO
!
        i = i + 1
        IF ( i==j ) GOTO 700
        t = Dx(i+1)
        ty = Dy(i+1)
        IF ( Dx(i)>t ) THEN
          k = i
          DO
!
            Dx(k+1) = Dx(k)
            Dy(k+1) = Dy(k)
            k = k - 1
            IF ( t>=Dx(k) ) THEN
              Dx(k+1) = t
              Dy(k+1) = ty
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!
!     Clean up
!
 900  IF ( Kflag<=-1 ) THEN
        DO i = 1 , nn
          Dx(i) = -Dx(i)
        ENDDO
      ENDIF
      END SUBROUTINE DSORT
