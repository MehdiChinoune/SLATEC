!DECK CVTCAT
CHARACTER(10) FUNCTION CVTCAT(Categ)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CVTCAT
  !***SUBSIDIARY
  !***PURPOSE  Expand a GAMS category
  !***LIBRARY   (NONE)
  !***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   Function to expand a GAMS category name by adding a zero before a
  !   one digit character.
  !   E.G.,  D2D1A   becomes   D02D01A
  !
  !***SEE ALSO  SLADOC
  !***ROUTINES CALLED  LENSTR
  !***REVISION HISTORY  (YYMMDD)
  !   871201  DATE WRITTEN
  !   891208  Prologue converted to Version 4.0 format.  (BAB)
  !   900202  Corrected error which occurred when expression I+1 exceeded
  !           length of input category (L).  (WRB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  CVTCAT
  !     .. Scalar Arguments ..
  CHARACTER*(*) Categ
  !     .. Local Scalars ..
  INTEGER i, ii, l
  !     .. External Functions ..
  INTEGER LENSTR
  EXTERNAL LENSTR
  !***FIRST EXECUTABLE STATEMENT  CVTCAT
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
    ENDIF
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
    ENDIF
  ENDDO
END FUNCTION CVTCAT
!DECK FIND
INTEGER FUNCTION FIND(X,N,T)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  FIND
  !***PURPOSE  Use a binary search to locate a phrase within a character
  !            array.
  !***LIBRARY   (NONE)
  !***CATEGORY  (NONE)
  !***KEYWORDS  BINARY SEARCH
  !***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory
  !***DESCRIPTION
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
  !        ENDIF
  !   10 CONTINUE
  !      FIND = 0
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891101  DATE WRITTEN
  !   891208  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  FIND
  !     .. Parameters ..
  REAL B
  PARAMETER (B=0.6931)
  !     .. Scalar Arguments ..
  INTEGER N
  CHARACTER*(*) T
  !     .. Array Arguments ..
  CHARACTER*(*) X(*)
  !     .. Local Scalars ..
  INTEGER i, k, m
  !     .. Intrinsic Functions ..
  INTRINSIC LOG, REAL
  !***FIRST EXECUTABLE STATEMENT  FIND
  IF ( N==0 ) THEN
    FIND = 0
    RETURN
  ENDIF
  m = LOG(REAL(N))/B
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
      ENDDO
    ELSE
      m = m - 1
      IF ( m<0 ) GOTO 100
      k = k/2
      i = i + k
    ENDIF
  ENDDO
  FIND = i
  GOTO 99999
  100  IF ( T<X(i) ) i = i - 1
  FIND = -i
  RETURN
  99999 CONTINUE
END FUNCTION FIND
!DECK LENSTR
INTEGER FUNCTION LENSTR(String)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  LENSTR
  !***PURPOSE  Find the position of the last non-blank character in a
  !            string.
  !***LIBRARY   SLATEC
  !***CATEGORY  N3
  !***TYPE      CHARACTER (LENSTR-H)
  !***KEYWORDS  POSITION OF LAST NON-BLANK CHARACTER
  !***AUTHOR  Berkbigler, Kathryn, C-3, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   This function finds the position of the last non-blank character in
  !   a string.
  !
  !   INPUT
  !
  !      STRING - a character variable which contains the string to be
  !               examined.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830124  DATE WRITTEN
  !   830324  REVISION DATE from the pre-1990 prologue.
  !   891208  Prologue converted to the 1990 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  LENSTR
  !     .. Scalar Arguments ..
  CHARACTER*(*) String
  !     .. Local Scalars ..
  INTEGER i, l
  !     .. Intrinsic Functions ..
  INTRINSIC LEN
  !***FIRST EXECUTABLE STATEMENT  LENSTR
  l = LEN(String)
  DO i = l, 1, -1
    IF ( String(i:i)/=' ' ) THEN
      LENSTR = i
      RETURN
    ENDIF
  ENDDO
  LENSTR = 0
END FUNCTION LENSTR
!DECK UPCASE
SUBROUTINE UPCASE(Str1,Str2)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  UPCASE
  !***PURPOSE  Convert all lower case alphabetic characters in a character
  !            string to upper case.
  !***LIBRARY   SLATEC
  !***CATEGORY  N3
  !***TYPE      CHARACTER (UPCASE-H)
  !***KEYWORDS  CONVERT CHARACTERS, UPPER CASE
  !***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  LENSTR, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   830308  DATE WRITTEN
  !   870902  REVISION DATE from the pre-1990 prologue.
  !   891208  Prologue converted to the 1990 format.  (BAB)
  !   901101  Correctly check the lengths of STR1 and STR2.  (WRB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  UPCASE
  !     .. Parameters ..
  CHARACTER(41) :: MSG
  PARAMETER (MSG='Input string is longer than output string')
  !     .. Scalar Arguments ..
  CHARACTER*(*) Str1, Str2
  !     .. Local Scalars ..
  INTEGER i, l, m, n
  !     .. External Functions ..
  INTEGER LENSTR
  EXTERNAL LENSTR
  !     .. External Subroutines ..
  EXTERNAL XERMSG
  !     .. Intrinsic Functions ..
  INTRINSIC CHAR, ICHAR, LEN
  !***FIRST EXECUTABLE STATEMENT  UPCASE
  m = LENSTR(Str1)
  n = LEN(Str2)
  IF ( m>n ) THEN
    CALL XERMSG(' ','UPCASE',MSG,1,2)
    m = n
  ENDIF
  DO i = 1, m
    Str2(i:i) = Str1(i:i)
    l = ICHAR(Str2(i:i))
    IF ( l>=97.AND.l<=122 ) Str2(i:i) = CHAR(l-32)
  ENDDO
  DO i = m + 1, n
    Str2(i:i) = ' '
  ENDDO
END SUBROUTINE UPCASE
