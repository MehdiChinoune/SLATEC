*DECK SLPREP
      PROGRAM SLPREP
C***BEGIN PROLOGUE  SLPREP
C***PURPOSE  Prepare one direct access file and three sequential files
C            for the SLATEC documentation program.
C***LIBRARY   (NONE)
C***CATEGORY  R4
C***KEYWORDS  DOCUMENTATION, SLADOC, SLATEC
C***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
C           Bacon, Barbara A., C-10, Los Alamos National Laboratory
C***DESCRIPTION
C
C   This program reads a sequential documentation file where each module
C   in the file consists of the complete subprogram statement and a
C   SLATEC-style prologue.  The program uses the information obtained to
C   generate four files.  These are
C
C   1)  a direct access file of the subprogram statement and prologue
C       for each routine in the library.
C   2)  a sequential file of routine names, categories, etc.,
C   3)  a sequential file of keywords and pointers to the routines.
C   4)  a sequential file of expanded categories and messages.
C
C   These four files constitute the database for the SLADOC
C   documentation program.
C
C   There are a number of system and library dependent parameters
C   which the user of this program may have to change before compiling
C   and running the code.  All parameters are defined in the records
C   which immediately follow this prologue.  In the discussion here, we
C   refer to the default values which are distributed with this code;
C   in order to assist others using this code, we give values for
C   several different machine/operating system configurations.
C
C      MXLFN  - the maximum length of a file name to be used.  The value
C               used is highly user and system dependent.  Set the value
C               to the length of longest of the 8 file names FOUT, FERR,
C               FINP, FDAF, FKWD, FTBL, FCLASS, and FCAT described
C               below.
C      FINP   - the name of the input file which contains the prologues.
C               Each prologue must be preceeded by a *DECK name record,
C               where name is the name of the subprogram, and a
C               subprogram procedure declaration statement, i.e. the
C               SUBROUTINE or FUNCTION statement.
C      FCLASS - the name of the input file which contains the GAMS
C               (SLATEC) classification file.
C      FCAT   - the name of the output sequential file which will
C               contain a linked list of all classifications used
C               by the SLATEC routines and their associated messages.
C      FDAF   - the name of the output direct access file which will
C               contain the documentation.
C      FKWD   - the name of the output sequential file which will
C               contain:
C               1) the number of keyword phrases in section 2.
C               2) the alphabetised KEYWORD phrases
C               3) a table of two pointer arrays:
C                    IL = {0 or next cell pointing to this keyword}
C                    IR = {location of the routine containing this
C                          keyword}
C               The first MXNKWD cells of this section correspond to
C               the entries in section 2 of this file.  Cells after that
C               are "continuation" cells.
C      FTBL   - the name of the output sequential file which will
C               contain:
C               1) name (first) category number
C               2) deck or routine name
C               3) line number of the subprogram statement
C               4) line number of the "END PROLOGUE" statement
C               5) line number of the PURPOSE statement
C      FOUT   - the name of the output file.  Some typical names are
C               tty (CTSS), OUTPUT (NOS), /dev/tty (UNIX) and
C               SYS$OUTPUT (VMS).
C      FLOG   - the name of the transaction log file.
C      FERR   - the name of the file which is to contain error
C               information.  All errors are processed by the XERMSG
C               package.  Some typical names are  tty (CTSS),
C               OUTPUT (NOS), /dev/tty (UNIX) and SYS$OUTPUT (VMS).
C
C      MXLRN  - the maximum length of a routine name.  For most Fortran
C               based libraries, including SLATEC, the value must be at
C               least 6.  If your library uses names longer than 6, you
C               should set the value of this parameter to the maximum
C               length.
C      MXNRN  - the maximum number of routine names in the library.
C      MXLCAT - the maximum length of a category number.  For the GAMS
C               classification scheme which is used by the SLATEC
C               Common Mathematical Library, the value is 10.
C      MXNCAT - the maximum number of categories in the entire library.
C      MXNKWD - the maximum number of keyword phrases in the entire
C               library.
C      MXNCL  - the maximum number of lines in the GAMS classification
C               file.
C      KMAXI  - the maximum number of characters in a keyword phrase.
C      KMAXJ  - the maximum number of keyword phrases in a subroutine.
C      LLN    - the maximum number of characters in an input line.
C
C***REFERENCES  Guide to the SLATEC Common Mathematical Library.
C***ROUTINES CALLED  CVTCAT, FIND, FTRUN, I1MACH, IFDECK, LENSTR, PSCAT,
C                    SETSIZ, SORT, UPCASE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   870818  DATE WRITTEN
C   880324  REVISION DATE from Version 3.2
C   891215  Prologue converted to Version 4.0 format.  (BAB)
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  SLPREP
C
C     System dependent parameter definitions.
C
      INTEGER MXLFN
      PARAMETER (MXLFN = 32)
      CHARACTER * (MXLFN) FINP, FCLASS, FCAT, FDAF, FKWD, FTBL, FOUT,
     +                    FLOG, FERR, FINPUT
      CHARACTER * (MXLFN) DFINP, DFCLAS, DFCAT, DFDAF, DFKWD, DFTBL,
     +                    DFLOG
      PARAMETER (DFINP   = 'slainp',
     +           DFCLAS  = 'class',
     +           DFCAT   = 'slacat',
     +           DFDAF   = 'sladaf',
     +           DFKWD   = 'slakwd',
     +           DFTBL   = 'slatbl',
     +           DFLOG   = 'slalog',
     +            FOUT   = '/dev/tty',
     +            FINPUT = '/dev/tty',
     +            FERR   = '/dev/tty')
C
C     Library dependent parameter definitions.
C
      INTEGER MXLRN, MXNRN, MXLCAT, MXNCAT, MXNKWD, MXNCL
      PARAMETER (MXLRN = 6, MXNRN = 1600, MXLCAT = 10, MXNCAT = 750,
     +           MXNKWD = 500, MXNCL = 751)
      INTEGER KMAXI, KMAXJ, LLN
      PARAMETER (KMAXI = 60, KMAXJ = 40, LLN = 72)
C
C     Other declarations.
C
      INTEGER I, IB, IC, ICHNG, ICOM, ID, IENTRY, IFIND, ILEN, INEXT,
     +        INFO, IPE, IPS, IRD, IWR, J, JJ, MNCL, MXLKW, MXLR, MXNCA,
     +        MXNKW, NCAT, NCC, NCLASS, NERR, NEXTL, NKWD, NSTMTS,
     +        NTCAT, NTKWD, NUMR, NUMRR
C
      INTEGER LU5, LU6, LU12, LU13, LU14, LU15, LU17, LU18, LU19
      PARAMETER (LU12 = 12, LU13 = 13, LU14 = 14, LU15 = 15, LU6 = 6,
     +           LU17 = 17, LU18 = 18, LU19 = 19, LU5 = 5)
      CHARACTER * (LLN) LINE
      CHARACTER * 80 CLLINE
      CHARACTER * 80 MSG
      CHARACTER * (KMAXI) KWRDS(KMAXJ), TKWD(MXNKWD)
      CHARACTER * (MXLRN) RTNAME
      CHARACTER * (MXLCAT) TCAT(MXNCAT), ETCAT(MXNCAT), CATEG(15),
     +                     TCLASS(MXNCAT)
      INTEGER IPTR(MXNCAT), JPTR(MXNCAT), KPTR(MXNCAT)
      CHARACTER * 80 CLASS(MXNCL+1), STMTS(MXNCL)
      INTEGER IPTRL(7*MXNRN), IPTRR(7*MXNRN)
C
C     External functions.
C
      INTEGER FIND, LENSTR
      LOGICAL IFDECK, IFIF, IFSID
      CHARACTER*10 CVTCAT
      EXTERNAL CVTCAT, FIND, IFDECK, IFIF, IFSID, LENSTR
C
C     Intrinsic functions.
C
      INTRINSIC ABS, INDEX, MAX, MOD
C***FIRST EXECUTABLE STATEMENT  SLPREP
C     OPEN (UNIT=LU6, FILE=FOUT, STATUS='UNKNOWN', FORM='FORMATTED',
C    +      IOSTAT = INFO)
C     IF (INFO .NE. 0) THEN
C       MSG = 'Failure in attempting to open ' // FOUT
C       NERR = 1
C       GO TO 240
C     ENDIF
C
C     OPEN (UNIT=LU5, FILE=FINPUT, STATUS='UNKNOWN', FORM='FORMATTED',
C    +      IOSTAT = INFO)
C     IF (INFO .NE. 0) THEN
C       MSG = 'Failure in attempting to open ' // FINPUT
C       NERR = 1
C       GO TO 240
C     ENDIF
C
      FINP = ' '
      WRITE (UNIT=LU6, FMT=9110) DFINP
      READ (UNIT=LU5, FMT=9000, END=12) FINP
   12 IF (LENSTR(FINP) .EQ. 0) FINP = DFINP
      OPEN (UNIT=LU15, FILE=FINP, STATUS='OLD', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FINP
        NERR = 1
        GO TO 240
      ENDIF
C
      FCLASS = ' '
      WRITE (UNIT=LU6, FMT=9090) DFCLAS
      READ (UNIT=LU5, FMT=9000, END=8) FCLASS
    8 IF (LENSTR(FCLASS) .EQ. 0) FCLASS = DFCLAS
      OPEN (UNIT=LU13, FILE=FCLASS, STATUS='OLD', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FCLASS
        NERR = 1
        GO TO 240
      ENDIF
C
      FCAT = ' '
      WRITE (UNIT=LU6, FMT=9100) DFCAT
      READ (UNIT=LU5, FMT=9000, END=10) FCAT
   10 IF (LENSTR(FCAT) .EQ. 0) FCAT = DFCAT
      OPEN (UNIT=LU14, FILE=FCAT, STATUS='NEW', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FCAT
        NERR = 1
        GO TO 240
      ENDIF
C
      FDAF = ' '
      WRITE (UNIT=LU6, FMT=9120) DFDAF
      READ (UNIT=LU5, FMT=9000, END=14) FDAF
   14 IF (LENSTR(FDAF) .EQ. 0) FDAF = DFDAF
      OPEN (UNIT=LU17, FILE=FDAF, STATUS='NEW', ACCESS='DIRECT',
     +      FORM='FORMATTED', RECL=LLN, IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FDAF
        NERR = 1
        GO TO 240
      ENDIF
C
      FTBL = ' '
      WRITE (UNIT=LU6,FMT=9130) DFTBL
      READ (UNIT=LU5,FMT=9000,END=16) FTBL
   16 IF (LENSTR(FTBL) .EQ. 0) FTBL = DFTBL
      OPEN (UNIT=LU18, FILE=FTBL, STATUS='NEW', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FTBL
        NERR = 1
        GO TO 240
      ENDIF
C
      FKWD = ' '
      WRITE (UNIT=LU6,FMT=9140) DFKWD
      READ (UNIT=LU5,FMT=9000,END=18) FKWD
   18 IF (LENSTR(FKWD) .EQ. 0) FKWD = DFKWD
      OPEN (UNIT=LU19, FILE=FKWD, STATUS='NEW', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FKWD
        NERR = 1
        GO TO 240
      ELSE
        NTKWD = 0
      ENDIF
C
      FLOG = ' '
      WRITE (UNIT=LU6, FMT=9080) DFLOG
      READ (UNIT=LU5, FMT=9000, END=6) FLOG
    6 IF (LENSTR(FLOG) .EQ. 0) FLOG = DFLOG
      OPEN (UNIT=LU12, FILE=FLOG, STATUS='NEW', FORM='FORMATTED',
     +      IOSTAT=INFO)
      IF (INFO .NE. 0) THEN
        MSG = 'Failure in attempting to open ' // FLOG
        NERR = 1
        GO TO 240
      ENDIF
C
      WRITE (UNIT=LU6, FMT=9070) FLOG
C
C     Write the names of all files to the transaction log file.
C
      WRITE (UNIT=LU12, FMT=9060) FINP, FOUT, FCLASS, FCAT, FDAF, FTBL,
     +                            FKWD
C
C     IRD is the "READ" line number and IWR is the "WRITE" line number.
C
      IRD = 0
      IWR = 0
C
C     NUMR is the number of routines in the library and MXLR is the
C     maximum length of all these routines.
C
      NUMR = 0
      NUMRR = 0
      MXLR = 0
      MXNCA = 0
      NTCAT = 0
      INEXT = MXNKWD+1
      MXNKW = 0
      MXLKW = 0
C
   20 READ (UNIT=LU15, FMT=9000, END=260) LINE
      IRD = IRD+1
      IF (IFDECK(LINE) .OR. IFIF(LINE)) GO TO 20
      IWR = IWR+1
      WRITE (UNIT=LU17, FMT=9000, REC=IWR) LINE
C
      IF (LINE(1:6) .EQ. '      ' .OR. IFSID(LINE)) THEN
C
C       Subprogram statement.  Save starting record number.
C
        IB = IWR
        NCAT = 1
        CATEG(1) = ' '
        NKWD = 0
      ENDIF
C
      IF (LINE(2:6) .EQ. '***BE') THEN
C
C       BEGIN PROLOGUE record.  Extract routine name.
C
        RTNAME = LINE(21:20+MXLRN)
C
C       Check case of routine name;  if not UPPER, terminate.
C
        CALL UPCASE (RTNAME, RTNAME)
        IF (RTNAME .NE. LINE(21:20+MXLRN)) THEN
          MSG = 'Routine name not in UPPER case'
          NERR = 2
          GO TO 240
        ENDIF
        NUMR = NUMR+1
        IF (NUMR .GT. MXNRN) THEN
          MSG = 'Too many routine names.  Recompile code with ' //
     +          'larger MXNRN.'
          NERR = 2
          GO TO 240
        ENDIF
C
C       Write every 10th routine name to standard output.
C
        IF (MOD(NUMR,10) .EQ. 0) WRITE (UNIT=LU6, FMT=9020) RTNAME, NUMR
      ENDIF
C
C     Effective with version 4.0, every routine MUST have a PURPOSE.
C
      IF (LINE(2:6) .EQ. '***PU') THEN
C
C       Record starting line number of Purpose section.
C
        IPS = IWR
   30   READ (UNIT=LU15, FMT=9000, END=260) LINE
        IRD = IRD + 1
        IWR = IWR + 1
        WRITE (UNIT=LU17, FMT=9000, REC=IWR) LINE
        IF (LINE(2:4) .NE. '***') GO TO 30
C
C       Record ending line number of Purpose section.
C
        IPE = IWR - 1
      ENDIF
C
      IF (LINE(2:6) .EQ. '***CA') THEN
C
C       "CATEGORY" record.
C
        NCAT = 0
C
C       Initialize the pointer to point to the first category.
C
   40   IC = 15
        ILEN = LENSTR(LINE)
        IF (IC .GT. ILEN) THEN
          MSG = 'No category on CATEGORY record'
          NERR = 3
          GO TO 240
        ENDIF
C
C       Get the offset location of the next delimiter.
C
   60   ICOM = INDEX(LINE(IC:ILEN),',')
        IF (ICOM .EQ. 0) THEN
          ID = ILEN+1
        ELSE
          ID = IC+ICOM-1
        ENDIF
        NCAT = NCAT+1
        IF (NCAT .GT. 15) THEN
          MSG = 'Too many categories in a routine.  Recompile code ' //
     +          'with larger MXNCAT.'
          NERR = 3
          GO TO 240
        ENDIF
C
C       Put category into table in UPPER case.
C
        CATEG(NCAT) = ' '
        CALL UPCASE (LINE(IC:ID-1), CATEG(NCAT)(1:ID-IC))
        IF (LINE(IC:ID-1) .NE. CATEG(NCAT)(1:ID-IC)) THEN
          MSG = 'Category not in UPPER case'
          NERR = 3
          GO TO 240
        ENDIF
C
C       Check to see if this category is already in the table of
C       categories.
C
        IFIND = FIND (TCAT, NTCAT, CATEG(NCAT))
        IF (IFIND .LE. 0) THEN
C
C         Check to see if there is room in the table.
C
          NTCAT = NTCAT+1
          IF (NTCAT .GT. MXNCAT) THEN
            MSG = 'Too many categories in the library.  Recompile ' //
     +            'code with larger MXNCAT.'
            NERR = 3
            GO TO 240
          ENDIF
C
C         Put category in table and check case.  If not UPPER case,
C         terminate.
C
          DO 70 I = NTCAT-1,ABS(IFIND)+1,-1
            TCAT(I+1) = TCAT(I)
   70     CONTINUE
          TCAT(ABS(IFIND)+1) = CATEG(NCAT)
        ENDIF
C
C       Is this the last category for this subprogram?
C
        IF (ICOM .EQ. 0) GO TO 80
C
C       Is this the last category on the current line?
C
        IF (ID .GE. ILEN) THEN
C
C         We are at the end of the line and need to read again.
C
          READ (UNIT=LU15, FMT=9000, END=260) LINE
          IRD = IRD+1
          IF (LINE(1:5).NE.'C   ' .AND. LINE(1:5).NE.'*   ') THEN
            MSG = 'CATEGORY section not in correct form'
            NERR = 3
            GO TO 240
          ENDIF
          IWR = IWR+1
          WRITE (UNIT=LU17, FMT=9000, REC=IWR) LINE
          GO TO 40
        ENDIF
C
C       Check to see that there is a blank character following the
C       delimiter after the category just processed.
C
        IF (LINE(ID+1:ID+1) .NE. ' ') THEN
          MSG = 'CATEGORY section not in correct form'
          NERR = 3
          GO TO 240
        ENDIF
C
C       Set the pointer IC to the next category.
C
        IC = ID+2
        GO TO 60
   80   MXNCA = MAX(MXNCA,NCAT)
      ENDIF
C
      IF (LINE(2:6) .EQ. '***KE') THEN
C
C       "KEYWORD" section.
C
C
C       Initialize the pointer to point to the first keyword phrase.
C
  100   IC = 15
        ILEN = LENSTR(LINE)
        IF (IC .GT. ILEN) THEN
          MSG = 'No keyword phrase on KEYWORD record'
          NERR = 4
          GO TO 240
        ENDIF
C
C       Get the offset location of the next delimiter.
C
  120   ICOM = INDEX(LINE(IC:ILEN),',')
        IF (ICOM .EQ. 0) THEN
          ID = ILEN+1
        ELSE
          ID = IC+ICOM-1
        ENDIF
        NKWD = NKWD+1
        IF (NKWD .GT. KMAXJ) THEN
          MSG = 'Too many keyword phrases in a routine.  Recompile ' //
     +          'code with larger KMAXJ.'
          NERR = 4
          GO TO 240
        ENDIF
        KWRDS(NKWD) = ' '
        IF (ID-IC .GT. KMAXI) THEN
          MSG = 'Keyword phrase is too long'
          NERR = 4
          GO TO 240
        ENDIF
C
C       Put keyword phrase in table and check case.  If not UPPER case,
C       terminate.
C
        CALL UPCASE (LINE(IC:ID-1), KWRDS(NKWD)(1:ID-IC))
        IF (LINE(IC:ID-1) .NE. KWRDS(NKWD)(1:ID-IC)) THEN
          MSG = 'Keyword phrase not in UPPER case'
          NERR = 4
          GO TO 240
        ENDIF
        MXLKW = MAX(MXLKW,ID-IC)
C
C       Is this the last keyword for this subprogram?
C
        IF (ICOM .EQ. 0) GO TO 140
        IF (ID .GE. ILEN) THEN
C
C         We are at the end of the line and need to read again.
C
          READ (UNIT=LU15, FMT=9000, END=260) LINE
          IRD = IRD+1
          IF (LINE(1:5).NE.'C   ' .AND. LINE(1:5).NE.'*   ') THEN
            MSG = 'KEYWORD section not in correct form'
            NERR = 4
            GO TO 240
          ENDIF
          IWR = IWR+1
          WRITE (UNIT=LU17, FMT=9000, REC=IWR) LINE
          GO TO 100
        ENDIF
C
C       Check to see that there is a blank character following the
C       delimiter after the category just processed.
C
        IF (LINE(ID+1:ID+1) .NE. ' ') THEN
          MSG = 'KEYWORD section not in correct form'
          NERR = 4
          GO TO 240
        ENDIF
        IC = ID+2
        GO TO 120
  140   MXNKW = MAX(MXNKW,NKWD)
      ENDIF
C
      IF (LINE(2:6) .EQ. '***EN') THEN
C
C       "END PROLOGUE" statement.
C
        IF (LINE(19:18+MXLRN) .EQ. RTNAME(1:MXLRN)) THEN
          MXLR = MAX(MXLR,IWR-IB+1)
C
C         Write the categories for this routine, the routine name, etc.
C         to the file FTBL.
C
          WRITE (UNIT=LU18, FMT=9010) (CATEG(J), RTNAME, IB, IWR, IPS,
     +                               IPE, J=1,NCAT)
          NUMRR = NUMRR+NCAT
C
          IF (NKWD .GT. 0) THEN
            MXNKW = MAX(MXNKW,NKWD)
            DO 200 J = NKWD,1,-1
              IFIND = FIND(TKWD(1), NTKWD, KWRDS(J))
              IF (IFIND .LE. 0) THEN
                NTKWD = NTKWD+1
                IF (NTKWD .GT. MXNKWD) THEN
                  MSG = 'Too many keyword phrases.  Recompile code ' //
     +                  'with larger MXNKWD.'
                  NERR=2
                  GO TO 240
                ENDIF
                DO 170 I = NTKWD-1,ABS(IFIND)+1,-1
                  TKWD(I+1) = TKWD(I)
                  IPTRR(I+1) = IPTRR(I)
                  IPTRL(I+1) = IPTRL(I)
  170           CONTINUE
                TKWD(ABS(IFIND)+1) = KWRDS(J)
                IPTRR(ABS(IFIND)+1) = NUMRR
                IPTRL(ABS(IFIND)+1) = 0
              ELSE
                IENTRY = IPTRL(IFIND)
                IF (IENTRY .EQ. 0) THEN
C
C                 No entries in the chain.
C
                  IPTRL(IFIND) = INEXT
                ELSE
  180             NEXTL = IPTRL(IENTRY)
                  IF (NEXTL .NE. 0) THEN
C
C                   Not at end of chain.
C
                    IENTRY = NEXTL
                    GO TO 180
                  ENDIF
                  IPTRL(IENTRY) = INEXT
                ENDIF
                IPTRL(INEXT) = 0
                IPTRR(INEXT) = NUMRR
                INEXT = INEXT+1
              ENDIF
  200       CONTINUE
          ENDIF
        ELSE
          MSG = 'BEGIN name and END name differ'
          NERR = 5
          GO TO 240
        ENDIF
C
C       Read to the next *DECK record.
C
  220   READ (UNIT=LU15, FMT=9000, END=260) LINE
        IRD = IRD+1
        IF (.NOT.IFDECK(LINE)) GO TO 220
      ENDIF
C
C     GOTO read.
C
      GO TO 20
C
C     Abnormal termination.
C
  240 CONTINUE
      CALL XERMSG (' ', 'SLPREP', MSG, NERR, 1)
      WRITE (UNIT=LU6, FMT=9050) RTNAME, IRD, LINE
C
C     Write the job abort message to the transaction log file also.
C
      WRITE (UNIT=LU12, FMT=9050) RTNAME, IRD, LINE
      CLOSE (UNIT=LU14, STATUS='DELETE')
      CLOSE (UNIT=LU17, STATUS='DELETE')
      CLOSE (UNIT=LU18, STATUS='DELETE')
      CLOSE (UNIT=LU19, STATUS='DELETE')
      GO TO 420
C
  260 DO 280 J = 1,NTCAT
        ETCAT(J) = CVTCAT(TCAT(J))
  280 CONTINUE
      CALL SORT (ETCAT, NTCAT, 1, TCAT)
C
C     Read in the GAMS classification file.
C
      I = 0
  300 READ (UNIT=LU13, FMT=9000, END = 320) CLLINE
      IF (I .GE. MXNCL) THEN
        MSG = 'Too many lines in file ' // FCLASS
        NERR = 6
        GO TO 240
      ENDIF
      I = I+1
      CLASS(I) = CLLINE
      GO TO 300
  320 NCLASS = I
      CLASS(I+1) = '@'
      MNCL = MXNCL
      CALL PSCAT (ETCAT, NTCAT, CLASS, MNCL, NCC, TCLASS, IPTR, JPTR,
     +            KPTR, NSTMTS, STMTS, NERR)
      IF (NERR .NE. 0) THEN
        MSG = 'Could not locate ' // TCLASS(NERR) // ' in GAMS file.'
        LINE = 'Input file completely read.'
        NERR = 7
        GO TO 240
      ENDIF
C
C     The reason NCC+1 is being written is to call attention to
C     the fact that an extra line is included in the KPTR array
C     which gives the line beyond the last line of the STMTS array.
C
      WRITE (UNIT=LU14, FMT=9003) NCC+1
      WRITE (UNIT=LU14, FMT=9005) (IPTR(I), JPTR(I), KPTR(I),
     +        TCLASS(I)(1:LENSTR(TCLASS(I))),I=1,NCC)
      WRITE (UNIT=LU14, FMT=9007) KPTR(NCC+1)
      WRITE (UNIT=LU14, FMT=9000)
     +        (STMTS(I)(1:LENSTR(STMTS(I))),I=1,NSTMTS)
C
      WRITE (UNIT=LU19, FMT=9003) NTKWD
      DO 340 J = 1,NTKWD
        ILEN = LENSTR(TKWD(J))
        WRITE (UNIT=LU19, FMT=9000) TKWD(J)(1:ILEN)
  340 CONTINUE
C
C     We now need to compress the IPTRL and IPTRR tables and remove
C     any unneeded cells between the original allocated cells (MXNKWD)
C     and the final number (NTKWD) of pointer cells.
C
      ICHNG=MXNKWD-NTKWD
      IF (ICHNG .GT. 0) THEN
        DO 360 J = 1,NTKWD
          IF (IPTRL(J) .NE. 0) THEN
            IPTRL(J)=IPTRL(J)-ICHNG
          ENDIF
  360   CONTINUE
        JJ=NTKWD+1
        DO 380 J = MXNKWD+1,INEXT-1
          IF (IPTRL(J) .NE. 0) THEN
            IPTRL(JJ)=IPTRL(J)-ICHNG
          ELSE
C           IPTRL(JJ)=IPTRL(J)
            IPTRL(JJ)=0
          ENDIF
          IPTRR(JJ)=IPTRR(J)
          JJ=JJ+1
  380   CONTINUE
      ENDIF
C
C     Set last set of pointers to 0.
C
      IPTRL(INEXT-ICHNG) = 0
      IPTRR(INEXT-ICHNG) = 0
      DO 400 J = 1,INEXT-ICHNG
        WRITE (UNIT=LU19, FMT=9003) IPTRL(J),IPTRR(J)
  400 CONTINUE
C
C     Normal termination.
C
      WRITE (UNIT=LU6, FMT=9030) NUMR, IRD, IWR, NTCAT, NTKWD, MXLR,
     +                           MXNCA, MXNKW, MXLKW, NUMRR, NTKWD,
     +                           INEXT-ICHNG
      WRITE (UNIT=LU6, FMT=9040) NCLASS, NCC, NSTMTS
C
C     Write summary information to the transaction log file also.
C
      WRITE (UNIT=LU12, FMT=9030) NUMR, IRD, IWR, NTCAT, NTKWD, MXLR,
     +                            MXNCA, MXNKW, MXLKW, NUMRR, NTKWD,
     +                            INEXT-ICHNG
      WRITE (UNIT=LU12, FMT=9040) NCLASS, NCC, NSTMTS
      CLOSE (UNIT=LU14)
C
C     Call to FTRUN for CTSS/CFTLIB systems only.
C
      CLOSE (UNIT=LU17)
      CLOSE (UNIT=LU18)
      CLOSE (UNIT=LU19)
C
  420 CLOSE (UNIT=LU12)
      CLOSE (UNIT=LU13)
      CLOSE (UNIT=LU15)
C     CLOSE (UNIT=LU5)
C     CLOSE (UNIT=LU6)
      STOP
C
 9000 FORMAT (A)
 9003 FORMAT (I5, 2X, I5)
 9005 FORMAT (3I5, 3X, A)
 9007 FORMAT (I15)
 9010 FORMAT (1X, 2A, 4I8)
 9020 FORMAT (' Processing routine ', A, ',  number ', I4)
 9030 FORMAT (// ' S U M M A R Y' //
     +        ' There were ', I4, ' documentation modules and a total ',
     +        'of ', I6, ' lines read ' /
     +        ' from the sequential input file and a total of ',
     +        I6, ' lines written ' /
     +        ' to the direct access documentation file.' /
     +        ' The library has a total of ', I3, ' distinct ',
     +        'categories and ', I3, ' distinct keywords.' /
     +        ' At least one module has ', I4, ' lines, one module ',
     +        'has ', I2, ' categories,' /
     +        ' one module has ', I2, ' keyword phrases and ',
     +        'the longest keyword phrase is ' /
     +        1X, I2, ' characters.' /
     +        ' There were ', I4, ' category/routine lines written to ',
     +         'the category file.' /
     +        ' There were ', I3, ' distinct keyword phrases written ',
     +        'to the keyword file, and ' /
     +        1X,  I4, ' keyword pointers written to the file.')
 9040 FORMAT (' There were ', I3, ' lines read from the GAMS ',
     +        'classification scheme file.' /
     +        ' There were ', I3, ' classification categories written ',
     +        'to the classification file' /
     +        ' and a total of ', I3, ' classification descriptors ',
     +        'written to the file.' //)
 9050 FORMAT (// ' J O B   A B O R T' //
     +        ' The job aborted in routine ', A, ' after processing ',
     +        'approximately ', I6, ' lines' / ' of the input file.' /
     +        ' The last line read from the input file was ' / 1X, A //)
 9060 FORMAT (// 5X, 'T R A N S A C T I O N   L O G' //
     +        2X, 'Input source file: ', A /
     +        2X, 'Output source file: ', A /
     +        2X, 'GAMS classification file: ', A /
     +        2X, 'Sequential file containing linked list of all ',
     +            'classifications used: ', A /
     +        2X, 'Direct access file containing the documentation: ',
     +             A /
     +        2X, 'Sequential file containing the deck names and ',
     +            'prologue line numbers: ', A /
     +        2X, 'Sequential file of keyword phrases and pointer ',
     +            'arrays: ', A/)
 9070 FORMAT (2X, 'Transaction log file: ', A /)
 9080 FORMAT (//, ' Enter the name of your transaction log file or',
     +            ' <cr>', /,
     +            ' (The default is ''', A, ''')')
 9090 FORMAT (//, ' Enter the name of the GAMS (SLATEC) classification',
     +            ' file or <cr>', /,
     =            ' (The default is ''', A, ''')')
 9100 FORMAT (//, ' Enter the name of the output file to contain the',
     +            ' classifications used ', /, ' or <cr>', /,
     +            ' (The default is ''', A, ''')')
 9110 FORMAT (//, ' Enter the name of your input source file or <cr>',/,
     +            ' (The default is ''', A, ''')')
 9120 FORMAT (//, ' Enter the name of the output direct access file',
     +            ' which will contain the ', /, ' documentation or',
     +            ' <cr>', /,
     +            ' (The default is ''', A, ''')')
 9130 FORMAT (//, ' Enter the name of the output file which will',
     +            ' contain the line numbers', /, ' of data in the',
     +            ' documentation file or <cr>', /,
     +            ' (The default is ''', A, ''')')
 9140 FORMAT (//, ' Enter the name of the output file which will',
     +            ' contain the keyword phrases ', /, ' and their',
     +            ' pointers or <cr>', /,
     +            ' (The default is ''', A, ''')')
      END
*DECK IFDECK
      LOGICAL FUNCTION IFDECK (LINE)
C***BEGIN PROLOGUE  IFDECK
C***PURPOSE  Determine if a line is an *DECK or *DK statement.
C***LIBRARY   SLATEC
C***CATEGORY  R3C
C***TYPE      LOGICAL (IFDECK-L)
C***KEYWORDS  DECK
C***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
C***DESCRIPTION
C
C   This FUNCTION determines if a record is an  *DECK or *DK  statement.
C   If it is, the value  .TRUE.  is returned; otherwise, the value
C   .FALSE. is returned.
C
C   INPUT
C
C      LINE   - a character variable containing the string to be
C               examined.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  UPCASE
C***REVISION HISTORY  (YYMMDD)
C   840427  DATE WRITTEN
C   840430  REVISION DATE from the pre-1990 prologue.
C   891215  Prologue converted to Version 4.0 format.  (BAB)
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  IFDECK
C     .. Scalar Arguments ..
      CHARACTER*(*) LINE
C     .. Local Scalars ..
      CHARACTER*6 TEMP
C     .. External Subroutines ..
      EXTERNAL UPCASE
C***FIRST EXECUTABLE STATEMENT  IFDECK
      CALL UPCASE (LINE(1:6), TEMP(1:6))
      IFDECK = .TRUE.
      IF (TEMP(1:6).NE.'*DECK ' .AND. TEMP(1:4).NE.'*DK ')
     +    IFDECK = .FALSE.
      RETURN
      END
*DECK IFIF
      LOGICAL FUNCTION IFIF (LINE)
C***BEGIN PROLOGUE  IFIF
C***PURPOSE  Determine if a line is an  *IF  statement.
C***LIBRARY   (NONE)
C***CATEGORY  (NONE)
C***KEYWORDS  (NONE)
C***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
C***DESCRIPTION
C
C   This FUNCTION determines if a record is an  *IF  statement.
C   If it is, the value  .TRUE.  is returned; otherwise, the value
C   .FALSE. is returned.
C
C   INPUT
C
C      LINE   - a character variable containing the string to be
C               examined.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  UPCASE
C***REVISION HISTORY  (YYMMDD)
C   910208  DATE WRITTEN
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  IFIF
C     .. Scalar Arguments ..
      CHARACTER*(*) LINE
C     .. Local Scalars ..
      CHARACTER*3 TEMP
C     .. External Subroutines ..
      EXTERNAL UPCASE
C***FIRST EXECUTABLE STATEMENT  IFIF
      CALL UPCASE (LINE(1:3), TEMP(1:3))
      IFIF = .TRUE.
      IF (TEMP(1:3).NE.'*IF') IFIF = .FALSE.
      RETURN
      END
*DECK IFSID
      LOGICAL FUNCTION IFSID (LINE)
C***BEGIN PROLOGUE  IFSID
C***PURPOSE  Determine if a line is a special *...IDENT statement.
C***LIBRARY   (NONE)
C***CATEGORY  (NONE)
C***KEYWORDS  (NONE)
C***AUTHOR  Bacon, Barbara, C-10, (LANL)
C***DESCRIPTION
C
C   This FUNCTION determines if a record is a  *...IDENT  statement.
C   If it is, the value  .TRUE.  is returned; otherwise, the value
C   .FALSE. is returned.
C
C   INPUT
C
C     LINE   - a character variable containing the string to be
C              examined.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  UPCASE
C***REVISION HISTORY  (YYMMDD)
C   910211  DATE WRITTEN
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  IFSID
C     .. Scalar Arguments ..
      CHARACTER*(*) LINE
C     .. Local Scalars ..
      INTEGER I, ID
      CHARACTER*20 TEMP
C     .. External Subroutines ..
      EXTERNAL UPCASE
C     .. Intrinsic Functions ..
      INTRINSIC INDEX
C***FIRST EXECUTABLE STATEMENT  IFSID
      CALL UPCASE (LINE(1:20), TEMP(1:20))
      IFSID = .FALSE.
      IF (LINE(1:1) .EQ. '*') THEN
        ID = INDEX(LINE,'IDENT')
        IF (ID .NE. 0) THEN
C         IF (LINE(2:ID-1) .EQ. ' ') IFSID = .TRUE.
          DO 10 I=2,ID-1
            IF (LINE(I:I) .NE. ' ') GO TO 20
   10     CONTINUE
C
C         We have a line with *^^...^^IDENT
C
          IFSID = .TRUE.
        ENDIF
      ENDIF
   20 RETURN
      END
*DECK PSCAT
      SUBROUTINE PSCAT (ECAT, NCAT, CLASS, MNCL, NCC, TCLASS, IPTR,
     +   JPTR, KPTR, ISTMT, STMTS, NERR)
C***BEGIN PROLOGUE  PSCAT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SLPREP
C***LIBRARY   (NONE)
C***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory
C***DESCRIPTION
C
C     Creates the data for the FCAT file consisting of two parts:
C     1) linked lists enabling one to obtain the various categories
C        and subcategories present in the library being examined
C     2) statements from the classification file associated with
C        the categories.
C
C   INPUT
C
C     ECAT   - an array of "extended" categories found in the library
C
C     NCAT   - the number of categories in the above arrays
C
C   OUTPUT
C
C     NCC    - the number of elements in TCLASS, IPTR, etc.
C
C     TCLASS - the expanded table of categories for this library;
C              e.g., if A1C were a category in the library, then
C              TCLASS would contain A, A1, and A1C.
C
C     IPTR   - the linked list allowing one to travel down categories
C              of the same length (i.e., A to C to ... to Z;
C              or  A1 to A3 to ... to An).
C
C     JPTR   - the linked list allowing one to get the next larger
C              element in a category (i.e., A to A1)
C
C     KPTR   - the location in the STMTS array of the statement
C              associated with the element in TCLASS.
C
C     ISTMT  - the number of statements in STMTS
C
C     STMTS  - the statements from the classification file associated
C              with TCLASS.
C
C     NERR   - the number of the element in TCLASS which could not
C              be located in the GAMS file.  This value is set to
C              zero upon entry to PSCAT and will indicate an error
C              to the calling routine if non-zero upon exit.
C
C     The easiest way to describe how this program works is to
C     use examples.
C     1) A02B05
C        At the beginning of the DO loop, the following values
C        will be set:
C        OPART1=OPART2=OPART3=OPART4=OPART5=OPART6=OPART7=' '
C        PART1='A'   PART2='02'   PART3='B'   PART4='05'
C        PART5=' '   PART6='  '   PART7=' '
C
C        Since PART1 .NE. OPART1, we will go to 100, clear all
C        the OPART's.  Then since OPART1 is not equal to PART1,
C        OPART1 will be set to PART1 and the first member of
C        TCLASS will be set to 'A'.  Then we will loop back to
C        statement 160 to build the rest of A02B05.
C        Now OPART1 is equal to PART1, so go to the next IF block.
C        Since OPART2=' ' and PART2='02', we will set OPART2=PART2,
C        and build the next member of TCLASS to be "A02" by
C        concatenating OPART1 and OPART2.
C        Then loop back to statement 160.
C        This process continues down through PART4, so that we will
C        build A02B and A02B05 in TCLASS.  So on the fifth time
C        down we will have OPART1='A', OPART2='02', OPART3='B',
C        OPART4='05', OPART5=OPART6=OPART7=' '.
C        Then since PART5=' ', we will transfer to the end of the
C        do 180 loop and return to the loop's beginning to pick up
C        a new member of ECAT.
C     2) A02B05C
C        At this point, the following values will be set when
C        we reach the first IF statement:
C        OPART1='A'  OPART2='02'  OPART3='B'  OPART4='05'
C        OPART5=' '  OPART6='  '  OPART7=' '
C        PART1 ='A'  PART2 ='02'  PART3 ='B'  PART4 ='05'
C        PART5 ='C'  PART6 =' '   PART7 =' '
C        So we will transfer to 100 and be sure that OPART5, OPART6,
C        and OPART7 are ' ' (it is possible for those to be non-blank
C        if the previous value of ECAT had been longer).
C        Then we will bypass the first four IF blocks and do number
C        five.  Since PART5 is not blank and OPART5.NE.PART5, we
C        will set OPART5=PART5 which is 'C' and build the next
C        member of TCLASS as "A02B05C' by concatenating OPART1,
C        OPART2, OPART3, OPART4, and OPART5.
C        Loop back to 160 where we will completely fall through all
C        the IF blocks to the end of the loop.
C     3) A03
C        Now PART1='A'  PART2='03' and PART3 through PART7 will be ' '.
C        But OPART1='A', OPART2='02', OPART3='B', OPART4='05',
C        OPART5='C' and OPART6=OPART7=' '.
C        This time we will branch to statement 40 since OPART1=PART1,
C        but OPART2 is not equal to PART2.
C        OPART2 through OPART7 will be set to ' '.
C        We will bypass the first IF block and do the second one,
C        since PART2 is not blank and OPART2.NE.PART2.
C        OPART2 will be set to PART2 (which is '03') and we will
C        build the next member of TCLASS by concatenating OPART1 and
C        OPART2 ('A03').  Then go back to statement 160.  All the IF
C        blocks will fail and we will drop through to the end of the
C        loop and go back to the beginning of the loop to get a new
C        value of ECAT.
C     4) This process continues until all the members of ECAT have been
C        processed.  The complete table is built in TCLASS and the
C        number of elements in it is NCC.
C
C
C     The next two DO blocks (the DO 280 and the DO 320) build the
C     IPTR and JPTR arrays, such that the members of IPTR point to
C     the next member of TCLASS which has the same length of "CLASS"
C     and the (n-1)th part of the member is the same but the n-th part
C     is different.  For example, the IPTR pointer associated with 'A'
C     would point to where 'C' is in TCLASS; then the IPTR pointer
C     associated with 'C' would point to where 'D' is, etc.  The final
C     IPTR pointer associated with a single letter in TCLASS will be
C     zero.  The IPTR pointer for A04 would point to where A06 is in
C     TCLASS and the IPTR pointer for A06 would be zero because there
C     are no more An.  The IPTR pointer for C01 points to where C02
C     is in TCLASS.  The IPTR pointer for C02 points to where C03 is.
C     ...  C14 is the last member of this type and its IPTR value is
C     zero.
C
C     The JPTR pointers are similar to the IPTR pointers, but point to
C     where the next sized member of TCLASS will be located.  For
C     example, the JPTR pointer for A points to where A04 is in TCLASS,
C     the JPTR pointer for A04 points to where A04A is in TCLASS, and
C     the JPTR pointer for A04A is zero since there is nothing more in
C     TCLASS which starts with A04A.
C
C     The final section puts into the array STMTS the statements
C     associated with each member of TCLASS.  Since a statement may take
C     up more than one line, we need another pointer table (KPTR) which
C     will point to the member of STMTS associated with each member of
C     TCLASS.  E.g. The statement for A is in STMTS(1), ..., the
C     statement for C01 is in STMTS(7), but the statement for C02 is
C     STMTS(9) because the one for C01 is two lines long.
C     By setting up KPTR this way (and putting a dummy value in the
C     cell beyond the last entry corresponding to a member of TCLASS),
C     we can find out how long a statement is by taking the KPTR pointer
C     for the member of TCLASS in which we are interested and subtract
C     it from the next member of KPTR.
C
C     The following is an excerpt from a file produced by this
C     code to show how it works.  Column 1 are the IPTR pointers,
C     column 2 are the JPTR pointers, column 3 are the KPTR pointers,
C     and column 4 are the members of TCLASS.
C
C       |<-------6       <----2       1   A
C       | |<-----4       |<-->3       2   A04
C       | |      0       |--->0       3   A04A
C       | |----->0       <----5       4   A06
C       |        0       |--->0       5   A06B
C       |<----->37       <----7       6   C
C       | |<-----8       |--->0       7   C01
C       | |<---->9            0       9   C02
C       | |<--->12       <---10      10   C03
C       | |      0       |<->11      11   C03A
C       | |      0       |--->0      12   C03A02
C       | |<--->16       <---13      13   C04
C       | | |--<14       |--->0      14   C04A
C       | | |<->15            0      15   C04B
C       | | |--->0            0      16   C04C
C       | |<--->17            0      17   C05
C       | |<--->23       <---18      18   C07
C       | | |<--19       |--->0      19   C07A
C       | | |<->20            0      20   C07B
C       | | |<->21            0      21   C07C
C       | | |<->22            0      22   C07E
C       | | |--->0            0      23   C07F
C       | |<--->26       <---24      24   C08
C       | | |--<25       |--->0      25   C08A
C       | | |--->0            0      27   C08C
C       | |<--->35       <---27      28   C10
C       | | |--<30       |<->28      29   C10A
C       | | | |<29       |--->0      30   C10A01
C       | | | |->0            0      31   C10A03
C       | | |<->33       <---31      32   C10B
C       | | | |<32       |--->0      33   C10B01
C       | | | |->0            0      34   C10B03
C       | | |<->34            0      35   C10D
C       | | |--->0            0      36   C10F
C       | |<--->36            0      37   C11
C       | |----->0            0      38   C14
C       |<---->123       <---38      39   D
C
C***ROUTINES CALLED  CVTCAT, LENSTR
C***REVISION HISTORY  (YYMMDD)
C   891215  DATE WRITTEN
C   891215  Prologue converted to Version 4.0 format.  (BAB)
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  PSCAT
C     .. Scalar Arguments ..
      INTEGER ISTMT, MNCL, NCAT, NCC, NERR
C     .. Array Arguments ..
      INTEGER IPTR(*), JPTR(*), KPTR(*)
      CHARACTER*(*) CLASS(*), ECAT(*), STMTS(*), TCLASS(*)
C     .. Local Scalars ..
      INTEGER I, ICLASS, ILEN, IPER, ISTART, J, K, NLEN
      CHARACTER*1 OPART1, OPART3, OPART5, OPART7, PART1, PART3,
     +            PART5, PART7
      CHARACTER*2 OPART2, OPART4, OPART6, PART2, PART4, PART6
C     .. Local Arrays ..
      INTEGER SIZE(7)
C     .. External Functions ..
      INTEGER LENSTR
      CHARACTER*10 CVTCAT
      EXTERNAL LENSTR,CVTCAT
C     .. Intrinsic Functions ..
      INTRINSIC INDEX
C     .. Data statements ..
      DATA SIZE /1, 3, 4, 6, 7, 9, 10/
C***FIRST EXECUTABLE STATEMENT  PSCAT
      NERR = 0
      NCC = 0
      OPART1 = ' '
      OPART2 = ' '
      OPART3 = ' '
      OPART4 = ' '
      OPART5 = ' '
      OPART6 = ' '
      OPART7 = ' '
      DO 180 I = 1,NCAT
        PART1 = ECAT(I)(1:1)
        PART2 = ECAT(I)(2:3)
        PART3 = ECAT(I)(4:4)
        PART4 = ECAT(I)(5:6)
        PART5 = ECAT(I)(7:7)
        PART6 = ECAT(I)(8:9)
        PART7 = ECAT(I)(10:10)
        IF (PART1 .NE. OPART1) GO TO 20
        IF (PART2 .NE. OPART2) GO TO 40
        IF (PART3 .NE. OPART3) GO TO 60
        IF (PART4 .NE. OPART4) GO TO 80
        IF (PART5 .NE. OPART5) GO TO 100
        IF (PART6 .NE. OPART6) GO TO 120
        IF (PART7 .NE. OPART7) GO TO 140
   20   OPART1 = ' '
   40   OPART2 = ' '
   60   OPART3 = ' '
   80   OPART4 = ' '
  100   OPART5 = ' '
  120   OPART6 = ' '
  140   OPART7 = ' '
  160   IF (OPART1 .NE. PART1) THEN
          OPART1 = PART1
          NCC = NCC+1
          TCLASS(NCC) = OPART1
        ENDIF
        IF (PART2 .EQ. ' ') GO TO 180
        IF (OPART2 .NE. PART2) THEN
          OPART2 = PART2
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2
        GO TO 160
        ENDIF
        IF (PART3 .EQ. ' ') GO TO 180
        IF (OPART3 .NE. PART3) THEN
          OPART3 = PART3
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2 // OPART3
        GO TO 160
        ENDIF
        IF (PART4 .EQ. ' ') GO TO 180
        IF (OPART4 .NE. PART4) THEN
          OPART4 = PART4
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2 // OPART3 // OPART4
        GO TO 160
        ENDIF
        IF (PART5 .EQ. ' ') GO TO 180
        IF (OPART5 .NE. PART5) THEN
          OPART5 = PART5
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2 // OPART3 // OPART4 // OPART5
        GO TO 160
        ENDIF
        IF (PART6 .EQ. ' ') GO TO 180
        IF (OPART6 .NE. PART6) THEN
          OPART6 = PART6
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2 // OPART3 // OPART4 //
     +               OPART5 // OPART6
        GO TO 160
        ENDIF
        IF (PART7 .EQ. ' ') GO TO 180
        IF (OPART7 .NE. PART7) THEN
          OPART7 = PART7
          NCC = NCC+1
          TCLASS(NCC) = OPART1 // OPART2 // OPART3 // OPART4 //
     +                  OPART5 // OPART6 // OPART7
        ENDIF
  180 CONTINUE
      DO 200 I = 1,NCC
        IPTR(I) = 0
        JPTR(I) = 0
  200 CONTINUE
      DO 280 J = 1,7
        ISTART = 1
  220   DO 260 I = ISTART,NCC
          ILEN = LENSTR(TCLASS(I))
          IF (ILEN .EQ. SIZE(J)) THEN
            DO 240 K = I+1,NCC
              NLEN = LENSTR(TCLASS(K))
              IF (ILEN .EQ. NLEN) THEN
                IF (J .EQ. 1) THEN
                  IPTR(I) = K
                ELSEIF (TCLASS(I)(1:SIZE(J-1)) .EQ.
     +                  TCLASS(K)(1:SIZE(J-1)))
     +                 THEN
                  IPTR(I) = K
                ELSE
                  IPTR(I) = 0
                ENDIF
                ISTART = K
                GO TO 220
              ENDIF
  240       CONTINUE
          ENDIF
  260   CONTINUE
  280 CONTINUE
      DO 320 J = 1,7
        DO 300 I = 1,NCC
          ILEN = LENSTR(TCLASS(I))
          IF (ILEN .EQ. SIZE(J)) THEN
            NLEN = LENSTR(TCLASS(I+1))
            IF (ILEN .LT. NLEN) THEN
              IF (TCLASS(I)(1:SIZE(J)) .EQ. TCLASS(I+1)(1:SIZE(J))) THEN
                JPTR(I) = I+1
              ENDIF
            ENDIF
          ENDIF
  300   CONTINUE
  320 CONTINUE
      DO 340 I = 1,NCC
        KPTR(I) = 0
  340 CONTINUE
C
      ICLASS = 1
      ISTMT = 1
      DO 400 I = 1,NCC
C
  360   IPER = INDEX(CLASS(ICLASS)(1:8),'.')
        IF (IPER .EQ. 0) THEN
          ICLASS = ICLASS+1
          IF (ICLASS.GT.MNCL) THEN
            NERR = I
            RETURN
          ENDIF
          GO TO 360
        ENDIF
C
        IF (TCLASS(I) .EQ. CVTCAT(CLASS(ICLASS)(1:IPER-1))) THEN
          KPTR(I) = ISTMT
  380     STMTS(ISTMT) = CLASS(ICLASS)(IPER+2:80)
          ISTMT = ISTMT+1
          ICLASS = ICLASS+1
          IF (CLASS(ICLASS)(1:1) .NE. ' ') GO TO 400
          GO TO 380
        ELSE
          ICLASS = ICLASS+1
          IF (ICLASS.GT.MNCL) THEN
            NERR = I
            RETURN
          ENDIF
          GO TO 360
        ENDIF
  400 CONTINUE
      KPTR(NCC+1) = ISTMT
      RETURN
      END
*DECK SORT
      SUBROUTINE SORT (R, N, NR, CR)
C***BEGIN PROLOGUE  SORT
C***PURPOSE  Alphabetise a character array and re-order one dependent
C            character array.
C***LIBRARY   (NONE)
C***CATEGORY  (NONE)
C***KEYWORDS  (NONE)
C***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory
C***DESCRIPTION
C
C   This subroutine alphabetises a character array.  It then accordingly
C   reorders up to one dependent character array.
C
C   INPUT
C
C      R      - the character array to be alphabetised.
C      N      - the number of elements in the array R.
C      NR     - the number (0 or 1) of dependent arrays.
C      CR     - the dependent array of length at most 15 characters.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   891101  DATE WRITTEN
C   891215  Prologue converted to Version 4.0 format.  (BAB)
C   920911  Declarations section restructured.  (WRB)
C***END PROLOGUE  SORT
C     .. Scalar Arguments ..
      INTEGER N, NR
C     .. Array Arguments ..
      CHARACTER*(*) CR(*), R(*)
C     .. Local Scalars ..
      INTEGER I, J, K, M
      CHARACTER*15 IT
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C***FIRST EXECUTABLE STATEMENT  SORT
      M = N
C     DO 5 I = 1,N
C       IR(I) = I
C   5 CONTINUE
   10 IF (M .GT. 15) GO TO 80
      M = M/2
      IF (M .EQ. 0) RETURN
   20 IF (MOD(M,2) .EQ. 0) M = M+1
      K = N-M
      J = 1
   30 I = J
   40 IF (R(I) .GT. R(I+M)) GO TO 60
   50 J = J+1
      IF (J .GT. K) GO TO 10
      GO TO 30
   60 IT = R(I)
      R(I) = R(I+M)
      R(I+M) = IT
      IF (NR .GT. 0) THEN
        IT = CR(I)
        CR(I) = CR(I+M)
        CR(I+M) = IT
      ENDIF
      I = I-M
      IF (I .LT. 1) GO TO 50
      GO TO 40
   80 M = M/4
      GO TO 20
      END
