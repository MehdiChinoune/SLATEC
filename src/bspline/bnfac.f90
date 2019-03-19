!** BNFAC
SUBROUTINE BNFAC(W,Nroww,Nrow,Nbandl,Nbandu,Iflag)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BINT4 and BINTK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BNFAC-S, DBNFAC-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  BNFAC is the BANFAC routine from
  !        * A Practical Guide to Splines *  by C. de Boor
  !
  !  Returns in  W  the lu-factorization (without pivoting) of the banded
  !  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
  !  onals in the work array  W .
  !
  !- ****  I N P U T  ******
  !  W.....Work array of size  (NROWW,NROW)  containing the interesting
  !        part of a banded matrix  A, with the diagonals or bands of  A
  !        stored in the rows of  W, while columns of  A  correspond to
  !        columns of  W . This is the storage mode used in  LINPACK  and
  !        results in efficient innermost loops.
  !           Explicitly,  A  has  NBANDL  bands below the diagonal
  !                            +     1     (main) diagonal
  !                            +   NBANDU  bands above the diagonal
  !        and thus, with    MIDDLE = NBANDU + 1,
  !          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
  !                                              J=1,...,NROW .
  !        For example, the interesting entries of A (1,2)-banded matrix
  !        of order  9  would appear in the first  1+1+2 = 4  rows of  W
  !        as follows.
  !                          13 24 35 46 57 68 79
  !                       12 23 34 45 56 67 78 89
  !                    11 22 33 44 55 66 77 88 99
  !                    21 32 43 54 65 76 87 98
  !
  !        All other entries of  W  not identified in this way with an en-
  !        try of  A  are never referenced .
  !  NROWW.....Row dimension of the work array  W .
  !        must be  .GE.  NBANDL + 1 + NBANDU  .
  !  NBANDL.....Number of bands of  A  below the main diagonal
  !  NBANDU.....Number of bands of  A  above the main diagonal .
  !
  !- ****  O U T P U T  ******
  !  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
  !     If  IFLAG = 1, then
  !  W.....contains the LU-factorization of  A  into a unit lower triangu-
  !        lar matrix  L  and an upper triangular matrix  U (both banded)
  !        and stored in customary fashion over the corresponding entries
  !        of  A . This makes it possible to solve any particular linear
  !        system  A*X = B  for  X  by A
  !              CALL BNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
  !        with the solution X  contained in  B  on return .
  !     If  IFLAG = 2, then
  !        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
  !        one of the potential pivots was found to be zero indicating
  !        that  A  does not have an LU-factorization. This implies that
  !        A  is singular in case it is totally positive .
  !
  !- ****  M E T H O D  ******
  !     Gauss elimination  W I T H O U T  pivoting is used. The routine is
  !  intended for use with matrices  A  which do not require row inter-
  !  changes during factorization, especially for the  T O T A L L Y
  !  P O S I T I V E  matrices which occur in spline calculations.
  !     The routine should not be used for an arbitrary banded matrix.
  !
  !***
  ! **See also:**  BINT4, BINTK
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  
  !
  INTEGER Iflag, Nbandl, Nbandu, Nrow, Nroww, i, ipk, j, jmax, k, &
    kmax, middle, midmk, nrowm1
  REAL W(Nroww,*), factor, pivot
  !
  !* FIRST EXECUTABLE STATEMENT  BNFAC
  Iflag = 1
  middle = Nbandu + 1
  !                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
  nrowm1 = Nrow - 1
  IF ( nrowm1<0 ) GOTO 100
  IF ( nrowm1/=0 ) THEN
    IF ( Nbandl<=0 ) THEN
      !                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO i = 1, nrowm1
        IF ( W(middle,i)==0.0E0 ) GOTO 100
      ENDDO
    ELSEIF ( Nbandu>0 ) THEN
      !
      !        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
      DO i = 1, nrowm1
        !                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        pivot = W(middle,i)
        IF ( pivot==0.0E0 ) GOTO 100
        !                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
        !                     BELOW THE DIAGONAL .
        jmax = MIN(Nbandl,Nrow-i)
        !              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO j = 1, jmax
          W(middle+j,i) = W(middle+j,i)/pivot
        ENDDO
        !                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
        !                     THE RIGHT OF THE DIAGONAL .
        kmax = MIN(Nbandu,Nrow-i)
        !                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
        !                  (BELOW ROW  I ) .
        DO k = 1, kmax
          ipk = i + k
          midmk = middle - k
          factor = W(midmk,ipk)
          DO j = 1, jmax
            W(midmk+j,ipk) = W(midmk+j,ipk) - W(middle+j,i)*factor
          ENDDO
        ENDDO
      ENDDO
    ELSE
      !              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
      !                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO i = 1, nrowm1
        pivot = W(middle,i)
        IF ( pivot==0.0E0 ) GOTO 100
        jmax = MIN(Nbandl,Nrow-i)
        DO j = 1, jmax
          W(middle+j,i) = W(middle+j,i)/pivot
        ENDDO
      ENDDO
      RETURN
    ENDIF
  ENDIF
  !                                       CHECK THE LAST DIAGONAL ENTRY .
  IF ( W(middle,Nrow)/=0.0E0 ) RETURN
  100 CONTINUE
  IFlag = 2
END SUBROUTINE BNFAC
