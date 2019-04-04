!** CVTCAT
CHARACTER(10) FUNCTION CVTCAT(Categ)
  IMPLICIT NONE
  !>
  !***
  !  Expand a GAMS category
  !***
  ! **Library:**   (NONE)
  !***
  ! **Author:**  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***
  ! **Description:**
  !
  !   Function to expand a GAMS category name by adding a zero before a
  !   one digit character.
  !   E.G.,  D2D1A   becomes   D02D01A
  !
  !***
  ! **See also:**  SLADOC
  !***
  ! **Routines called:**  LENSTR

  !* REVISION HISTORY  (YYMMDD)
  !   871201  DATE WRITTEN
  !   891208  Prologue converted to Version 4.0 format.  (BAB)
  !   900202  Corrected error which occurred when expression I+1 exceeded
  !           length of input category (L).  (WRB)
  !   920911  Declarations section restructured.  (WRB)

  !     .. Scalar Arguments ..
  CHARACTER*(*) Categ
  !     .. Local Scalars ..
  INTEGER i, ii, l
  !     .. External Functions ..
  INTEGER, EXTERNAL :: LENSTR
  !* FIRST EXECUTABLE STATEMENT  CVTCAT
  l = LENSTR(Categ)
  CVTCAT(1:1) = Categ(1:1)
  CVTCAT(2:10) = '        '
  i = 2
  ii = 2
  DO WHILE ( i<=l )
    IF ( i==l ) THEN
      CVTCAT(ii:ii) = '0'
      CVTCAT(ii+1:ii+1) = Categ(i:i)
      i = i + 1
      CYCLE
    END IF
    IF ( Categ(i+1:i+1)>'9'.OR.Categ(i+1:i+1)<'0' ) THEN
      CVTCAT(ii:ii) = '0'
      CVTCAT(ii+1:ii+1) = Categ(i:i)
      CVTCAT(ii+2:ii+2) = Categ(i+1:i+1)
      i = i + 2
      ii = ii + 3
    ELSE
      CVTCAT(ii:ii) = Categ(i:i)
      CVTCAT(ii+1:ii+1) = Categ(i+1:i+1)
      IF ( i+1<l ) CVTCAT(ii+2:ii+2) = Categ(i+2:i+2)
      i = i + 3
      ii = ii + 3
    END IF
  END DO
END FUNCTION CVTCAT
!** FIND
INTEGER FUNCTION FIND(X,N,T)
  IMPLICIT NONE
  !>
  !***
  !  Use a binary search to locate a phrase within a character
  !            array.
  !***
  ! **Library:**   (NONE)
  !***
  ! **Category:**  (NONE)
  !***
  ! **Keywords:**  BINARY SEARCH
  !***
  ! **Author:**  Bacon, Barbara A., C-10, Los Alamos National Laboratory
  !***
  ! **Description:**
  !
  !   This function attempts to find the member of the array matching
  !   the input phrase.  If the phrase is located, the subscript of the
  !   array element matching the input phrase is returned; if the table
  !   is empty (N=0) or the input phrase is alphabetically less than
  !   the first entry in the table, a 0 (zero) is returned; otherwise,
  !   the negative of the subscript of the location in the table of the
  !   phrase alphabetically less than the input phrase is returned.
  !
  !   INPUT
  !
  !      X      - the array to be searched.
  !      N      - the number of elements in the array X.
  !      T      - the phrase being searched for.
  !
  !   The routine uses a binary search strategy to implement the following
  !   Fortran sequence:
  !
  !      DO 10 I = 1,N
  !        IF (T .EQ. X(I)) THEN
  !          FIND = I
  !          RETURN
  !        END IF
  !   10 CONTINUE
  !      FIND = 0
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   891101  DATE WRITTEN
  !   891208  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)

  !     .. Parameters ..
  REAL, PARAMETER :: B = 0.6931
  !     .. Scalar Arguments ..
  INTEGER N
  CHARACTER*(*) T
  !     .. Array Arguments ..
  CHARACTER*(*) X(*)
  !     .. Local Scalars ..
  INTEGER i, k, m
  !     .. Intrinsic Functions ..
  INTRINSIC LOG, REAL
  !* FIRST EXECUTABLE STATEMENT  FIND
  IF ( N==0 ) THEN
    FIND = 0
    RETURN
  END IF
  m = INT( LOG(REAL(N))/B )
  k = 2**m
  i = N - k + 1
  DO WHILE ( T/=X(i) )
    IF ( T<X(i) ) THEN
      DO
        m = m - 1
        IF ( m<0 ) GOTO 100
        k = k/2
        i = i - k
        IF ( i>=1 ) EXIT
        i = i + k
      END DO
    ELSE
      m = m - 1
      IF ( m<0 ) GOTO 100
      k = k/2
      i = i + k
    END IF
  END DO
  FIND = i
  RETURN
  100  IF ( T<X(i) ) i = i - 1
  FIND = -i
  RETURN
END FUNCTION FIND
!** LENSTR
INTEGER FUNCTION LENSTR(String)
  IMPLICIT NONE
  !>
  !***
  !  Find the position of the last non-blank character in a
  !            string.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N3
  !***
  ! **Type:**      CHARACTER (LENSTR-H)
  !***
  ! **Keywords:**  POSITION OF LAST NON-BLANK CHARACTER
  !***
  ! **Author:**  Berkbigler, Kathryn, C-3, Los Alamos National Laboratory
  !***
  ! **Description:**
  !
  !   This function finds the position of the last non-blank character in
  !   a string.
  !
  !   INPUT
  !
  !      STRING - a character variable which contains the string to be
  !               examined.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830124  DATE WRITTEN
  !   830324  REVISION DATE from the pre-1990 prologue.
  !   891208  Prologue converted to the 1990 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)

  !     .. Scalar Arguments ..
  CHARACTER*(*) String
  !     .. Local Scalars ..
  INTEGER i, l
  !     .. Intrinsic Functions ..
  INTRINSIC LEN
  !* FIRST EXECUTABLE STATEMENT  LENSTR
  l = LEN(String)
  DO i = l, 1, -1
    IF ( String(i:i)/=' ' ) THEN
      LENSTR = i
      RETURN
    END IF
  END DO
  LENSTR = 0
END FUNCTION LENSTR
!** UPCASE
SUBROUTINE UPCASE(Str1,Str2)
  IMPLICIT NONE
  !>
  !***
  !  Convert all lower case alphabetic characters in a character
  !            string to upper case.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  N3
  !***
  ! **Type:**      CHARACTER (UPCASE-H)
  !***
  ! **Keywords:**  CONVERT CHARACTERS, UPPER CASE
  !***
  ! **Author:**  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***
  ! **Description:**
  !
  !   This routine converts all the lower case alphabetic characters
  !   in a character string to UPPER case.  If the output string is
  !   longer than the input string, it is padded with blanks.  If the
  !   output string is not long enough, an error message is written and
  !   the input string is truncated.
  !
  !   INPUT
  !
  !     STR1   - a character string which is to be converted to
  !              upper case.
  !
  !   OUTPUT
  !
  !     STR2   - a character string which has all alphabetic
  !              characters in upper case.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  LENSTR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   830308  DATE WRITTEN
  !   870902  REVISION DATE from the pre-1990 prologue.
  !   891208  Prologue converted to the 1990 format.  (BAB)
  !   901101  Correctly check the lengths of STR1 and STR2.  (WRB)
  !   920911  Declarations section restructured.  (WRB)

  !     .. Parameters ..
  CHARACTER(41), PARAMETER :: MSG = 'Input string is longer than output string'
  !     .. Scalar Arguments ..
  CHARACTER*(*) Str1, Str2
  !     .. Local Scalars ..
  INTEGER i, l, m, n
  !     .. External Functions ..
  INTEGER, EXTERNAL :: LENSTR
  !     .. External Subroutines ..
  EXTERNAL :: XERMSG
  !     .. Intrinsic Functions ..
  INTRINSIC CHAR, ICHAR, LEN
  !* FIRST EXECUTABLE STATEMENT  UPCASE
  m = LENSTR(Str1)
  n = LEN(Str2)
  IF ( m>n ) THEN
    CALL XERMSG(' ','UPCASE',MSG,1,2)
    m = n
  END IF
  DO i = 1, m
    Str2(i:i) = Str1(i:i)
    l = ICHAR(Str2(i:i))
    IF ( l>=97.AND.l<=122 ) Str2(i:i) = CHAR(l-32)
  END DO
  DO i = m + 1, n
    Str2(i:i) = ' '
  END DO
END SUBROUTINE UPCASE
