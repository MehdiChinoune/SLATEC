*DECK CVTCAT
      CHARACTER * 10 FUNCTION CVTCAT (CATEG)    
C***BEGIN PROLOGUE  CVTCAT      
C***SUBSIDIARY  
C***PURPOSE  Expand a GAMS category     
C***LIBRARY   (NONE)    
C***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory      
C***DESCRIPTION 
C       
C   Function to expand a GAMS category name by adding a zero before a   
C   one digit character.
C   E.G.,  D2D1A   becomes   D02D01A    
C       
C***SEE ALSO  SLADOC    
C***ROUTINES CALLED  LENSTR     
C***REVISION HISTORY  (YYMMDD)  
C   871201  DATE WRITTEN
C   891208  Prologue converted to Version 4.0 format.  (BAB)    
C   900202  Corrected error which occurred when expression I+1 exceeded 
C           length of input category (L).  (WRB)
C   920911  Declarations section restructured.  (WRB)   
C***END PROLOGUE  CVTCAT
C     .. Scalar Arguments ..    
      CHARACTER*(*) CATEG       
C     .. Local Scalars ..       
      INTEGER I, II, L  
C     .. External Functions ..  
      INTEGER LENSTR    
      EXTERNAL LENSTR   
C***FIRST EXECUTABLE STATEMENT  CVTCAT  
      L = LENSTR(CATEG) 
      CVTCAT(1:1) = CATEG(1:1)  
      CVTCAT(2:10) = '        ' 
      I = 2     
      II = 2    
   10 IF (I .LE. L) THEN
        IF (I .EQ. L) THEN      
          CVTCAT(II:II) = '0'   
          CVTCAT(II+1:II+1) = CATEG(I:I)
          I = I+1       
          GO TO 20      
        ENDIF   
        IF (CATEG(I+1:I+1).GT.'9' .OR. CATEG(I+1:I+1).LT.'0') THEN      
          CVTCAT(II:II) = '0'   
          CVTCAT(II+1:II+1) = CATEG(I:I)
          CVTCAT(II+2:II+2) = CATEG(I+1:I+1)    
          I = I+2       
          II = II+3     
        ELSE    
          CVTCAT(II:II) = CATEG(I:I)    
          CVTCAT(II+1:II+1) = CATEG(I+1:I+1)    
          IF (I+1 .LT. L) THEN  
            CVTCAT(II+2:II+2) = CATEG(I+2:I+2)  
          ENDIF 
          I = I+3       
          II = II+3     
        ENDIF   
   20   GO TO 10
      ENDIF     
      RETURN    
      END       
*DECK FIND
      INTEGER FUNCTION FIND (X, N, T)   
C***BEGIN PROLOGUE  FIND
C***PURPOSE  Use a binary search to locate a phrase within a character  
C            array.     
C***LIBRARY   (NONE)    
C***CATEGORY  (NONE)    
C***KEYWORDS  BINARY SEARCH     
C***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory     
C***DESCRIPTION 
C       
C   This function attempts to find the member of the array matching     
C   the input phrase.  If the phrase is located, the subscript of the   
C   array element matching the input phrase is returned; if the table   
C   is empty (N=0) or the input phrase is alphabetically less than      
C   the first entry in the table, a 0 (zero) is returned; otherwise,    
C   the negative of the subscript of the location in the table of the   
C   phrase alphabetically less than the input phrase is returned.       
C       
C   INPUT       
C       
C      X      - the array to be searched.       
C      N      - the number of elements in the array X.  
C      T      - the phrase being searched for.  
C       
C   The routine uses a binary search strategy to implement the following
C   Fortran sequence:   
C       
C      DO 10 I = 1,N    
C        IF (T .EQ. X(I)) THEN  
C          FIND = I     
C          RETURN       
C        ENDIF  
C   10 CONTINUE 
C      FIND = 0 
C       
C***REFERENCES  (NONE)  
C***ROUTINES CALLED  (NONE)     
C***REVISION HISTORY  (YYMMDD)  
C   891101  DATE WRITTEN
C   891208  Prologue converted to Version 4.0 format.  (BAB)    
C   920911  Declarations section restructured.  (WRB)   
C***END PROLOGUE  FIND  
C     .. Parameters ..  
      REAL B    
      PARAMETER (B = 0.6931)    
C     .. Scalar Arguments ..    
      INTEGER N 
      CHARACTER*(*) T   
C     .. Array Arguments ..     
      CHARACTER*(*) X(*)
C     .. Local Scalars ..       
      INTEGER I, K, M   
C     .. Intrinsic Functions .. 
      INTRINSIC LOG, REAL       
C***FIRST EXECUTABLE STATEMENT  FIND    
      IF (N .EQ. 0) THEN
        FIND = 0
        RETURN  
      ENDIF     
      M = LOG(REAL(N))/B
      K = 2**M  
      I = N-K+1 
   10 IF (T .EQ. X(I)) GO TO 40 
      IF (T .LT. X(I)) GO TO 20 
      M = M-1   
      IF (M .LT. 0) GO TO 30    
      K = K/2   
      I = I+K   
      GO TO 10  
   20 M = M-1   
      IF (M .LT. 0) GO TO 30    
      K = K/2   
      I = I-K   
      IF (I .GE. 1) GO TO 10    
      I = I+K   
      GO TO 20  
   30 CONTINUE  
      IF (T .LT. X(I)) I = I-1  
      FIND = -I 
      RETURN    
   40 FIND = I  
      RETURN    
      END       
*DECK LENSTR
      INTEGER FUNCTION LENSTR (STRING)  
C***BEGIN PROLOGUE  LENSTR      
C***PURPOSE  Find the position of the last non-blank character in a     
C            string.    
C***LIBRARY   SLATEC    
C***CATEGORY  N3
C***TYPE      CHARACTER (LENSTR-H)      
C***KEYWORDS  POSITION OF LAST NON-BLANK CHARACTER      
C***AUTHOR  Berkbigler, Kathryn, C-3, Los Alamos National Laboratory    
C***DESCRIPTION 
C       
C   This function finds the position of the last non-blank character in 
C   a string.   
C       
C   INPUT       
C       
C      STRING - a character variable which contains the string to be    
C               examined.       
C       
C***REFERENCES  (NONE)  
C***ROUTINES CALLED  (NONE)     
C***REVISION HISTORY  (YYMMDD)  
C   830124  DATE WRITTEN
C   830324  REVISION DATE from the pre-1990 prologue.   
C   891208  Prologue converted to the 1990 format.  (BAB)       
C   920911  Declarations section restructured.  (WRB)   
C***END PROLOGUE  LENSTR
C     .. Scalar Arguments ..    
      CHARACTER*(*) STRING      
C     .. Local Scalars ..       
      INTEGER I, L      
C     .. Intrinsic Functions .. 
      INTRINSIC LEN     
C***FIRST EXECUTABLE STATEMENT  LENSTR  
      L = LEN(STRING)   
      DO 10 I = L,1,-1  
        IF (STRING(I:I) .NE. ' ') THEN  
          LENSTR = I    
          RETURN
        ENDIF   
   10 CONTINUE  
      LENSTR = 0
      RETURN    
      END       
*DECK UPCASE
      SUBROUTINE UPCASE (STR1, STR2)    
C***BEGIN PROLOGUE  UPCASE      
C***PURPOSE  Convert all lower case alphabetic characters in a character
C            string to upper case.      
C***LIBRARY   SLATEC    
C***CATEGORY  N3
C***TYPE      CHARACTER (UPCASE-H)      
C***KEYWORDS  CONVERT CHARACTERS, UPPER CASE    
C***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory      
C***DESCRIPTION 
C       
C   This routine converts all the lower case alphabetic characters      
C   in a character string to UPPER case.  If the output string is       
C   longer than the input string, it is padded with blanks.  If the     
C   output string is not long enough, an error message is written and   
C   the input string is truncated.      
C       
C   INPUT       
C       
C     STR1   - a character string which is to be converted to   
C              upper case.      
C       
C   OUTPUT      
C       
C     STR2   - a character string which has all alphabetic      
C              characters in upper case.
C       
C***REFERENCES  (NONE)  
C***ROUTINES CALLED  LENSTR, XERMSG     
C***REVISION HISTORY  (YYMMDD)  
C   830308  DATE WRITTEN
C   870902  REVISION DATE from the pre-1990 prologue.   
C   891208  Prologue converted to the 1990 format.  (BAB)       
C   901101  Correctly check the lengths of STR1 and STR2.  (WRB)
C   920911  Declarations section restructured.  (WRB)   
C***END PROLOGUE  UPCASE
C     .. Parameters ..  
      CHARACTER*41 MSG  
      PARAMETER (MSG = 'Input string is longer than output string')     
C     .. Scalar Arguments ..    
      CHARACTER*(*) STR1, STR2  
C     .. Local Scalars ..       
      INTEGER I, L, M, N
C     .. External Functions ..  
      INTEGER LENSTR    
      EXTERNAL LENSTR   
C     .. External Subroutines ..
      EXTERNAL XERMSG   
C     .. Intrinsic Functions .. 
      INTRINSIC CHAR, ICHAR, LEN
C***FIRST EXECUTABLE STATEMENT  UPCASE  
      M = LENSTR(STR1)  
      N = LEN(STR2)     
      IF (M .GT. N) THEN
        CALL XERMSG (' ', 'UPCASE', MSG, 1, 2)  
        M = N   
      ENDIF     
      DO 10 I = 1,M     
        STR2(I:I) = STR1(I:I)   
        L = ICHAR(STR2(I:I))    
        IF (L.GE.97 .AND. L.LE.122) STR2(I:I) = CHAR(L-32)      
   10 CONTINUE  
      DO 20 I = M+1,N   
        STR2(I:I) = ' ' 
   20 CONTINUE  
      RETURN    
      END       
