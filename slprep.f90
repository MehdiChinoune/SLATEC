!DECK SLPREP
PROGRAM SLPREP
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SLPREP
  !***PURPOSE  Prepare one direct access file and three sequential files
  !            for the SLATEC documentation program.
  !***LIBRARY   (NONE)
  !***CATEGORY  R4
  !***KEYWORDS  DOCUMENTATION, SLADOC, SLATEC
  !***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !           Bacon, Barbara A., C-10, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   This program reads a sequential documentation file where each module
  !   in the file consists of the complete subprogram statement and a
  !   SLATEC-style prologue.  The program uses the information obtained to
  !   generate four files.  These are
  !
  !   1)  a direct access file of the subprogram statement and prologue
  !       for each routine in the library.
  !   2)  a sequential file of routine names, categories, etc.,
  !   3)  a sequential file of keywords and pointers to the routines.
  !   4)  a sequential file of expanded categories and messages.
  !
  !   These four files constitute the database for the SLADOC
  !   documentation program.
  !
  !   There are a number of system and library dependent parameters
  !   which the user of this program may have to change before compiling
  !   and running the code.  All parameters are defined in the records
  !   which immediately follow this prologue.  In the discussion here, we
  !   refer to the default values which are distributed with this code;
  !   in order to assist others using this code, we give values for
  !   several different machine/operating system configurations.
  !
  !      MXLFN  - the maximum length of a file name to be used.  The value
  !               used is highly user and system dependent.  Set the value
  !               to the length of longest of the 8 file names FOUT, FERR,
  !               FINP, FDAF, FKWD, FTBL, FCLASS, and FCAT described
  !               below.
  !      FINP   - the name of the input file which contains the prologues.
  !               Each prologue must be preceeded by a *DECK name record,
  !               where name is the name of the subprogram, and a
  !               subprogram procedure declaration statement, i.e. the
  !               SUBROUTINE or FUNCTION statement.
  !      FCLASS - the name of the input file which contains the GAMS
  !               (SLATEC) classification file.
  !      FCAT   - the name of the output sequential file which will
  !               contain a linked list of all classifications used
  !               by the SLATEC routines and their associated messages.
  !      FDAF   - the name of the output direct access file which will
  !               contain the documentation.
  !      FKWD   - the name of the output sequential file which will
  !               contain:
  !               1) the number of keyword phrases in section 2.
  !               2) the alphabetised KEYWORD phrases
  !               3) a table of two pointer arrays:
  !                    IL = {0 or next cell pointing to this keyword}
  !                    IR = {location of the routine containing this
  !                          keyword}
  !               The first MXNKWD cells of this section correspond to
  !               the entries in section 2 of this file.  Cells after that
  !               are "continuation" cells.
  !      FTBL   - the name of the output sequential file which will
  !               contain:
  !               1) name (first) category number
  !               2) deck or routine name
  !               3) line number of the subprogram statement
  !               4) line number of the "END PROLOGUE" statement
  !               5) line number of the PURPOSE statement
  !      FOUT   - the name of the output file.  Some typical names are
  !               tty (CTSS), OUTPUT (NOS), /dev/tty (UNIX) and
  !               SYS$OUTPUT (VMS).
  !      FLOG   - the name of the transaction log file.
  !      FERR   - the name of the file which is to contain error
  !               information.  All errors are processed by the XERMSG
  !               package.  Some typical names are  tty (CTSS),
  !               OUTPUT (NOS), /dev/tty (UNIX) and SYS$OUTPUT (VMS).
  !
  !      MXLRN  - the maximum length of a routine name.  For most Fortran
  !               based libraries, including SLATEC, the value must be at
  !               least 6.  If your library uses names longer than 6, you
  !               should set the value of this parameter to the maximum
  !               length.
  !      MXNRN  - the maximum number of routine names in the library.
  !      MXLCAT - the maximum length of a category number.  For the GAMS
  !               classification scheme which is used by the SLATEC
  !               Common Mathematical Library, the value is 10.
  !      MXNCAT - the maximum number of categories in the entire library.
  !      MXNKWD - the maximum number of keyword phrases in the entire
  !               library.
  !      MXNCL  - the maximum number of lines in the GAMS classification
  !               file.
  !      KMAXI  - the maximum number of characters in a keyword phrase.
  !      KMAXJ  - the maximum number of keyword phrases in a subroutine.
  !      LLN    - the maximum number of characters in an input line.
  !
  !***REFERENCES  Guide to the SLATEC Common Mathematical Library.
  !***ROUTINES CALLED  CVTCAT, FIND, FTRUN, I1MACH, IFDECK, LENSTR, PSCAT,
  !                    SETSIZ, SORT, UPCASE, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   870818  DATE WRITTEN
  !   880324  REVISION DATE from Version 3.2
  !   891215  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  SLPREP
  !
  !     System dependent parameter definitions.
  !
  INTEGER MXLFN
  PARAMETER (MXLFN=32)
  CHARACTER(MXLFN) :: finp, fclass, fcat, fdaf, fkwd, ftbl, FOUT ,&
    flog, FERR, FINPUT
  CHARACTER(MXLFN) :: DFINP, DFCLAS, DFCAT, DFDAF, DFKWD, DFTBL, DFLOG
  PARAMETER (DFINP='slainp',DFCLAS='class',DFCAT='slacat',DFDAF='sladaf',&
    DFKWD='slakwd',DFTBL='slatbl',DFLOG='slalog',FOUT='/dev/tty',&
    FINPUT='/dev/tty',FERR='/dev/tty')
  !
  !     Library dependent parameter definitions.
  !
  INTEGER MXLRN, MXNRN, MXLCAT, MXNCAT, MXNKWD, MXNCL
  PARAMETER (MXLRN=6,MXNRN=1600,MXLCAT=10,MXNCAT=750,MXNKWD=500,MXNCL=751)
  INTEGER KMAXI, KMAXJ, LLN
  PARAMETER (KMAXI=60,KMAXJ=40,LLN=72)
  !
  !     Other declarations.
  !
  INTEGER i, ib, ic, ichng, icom, id, ientry, ifind, ilen, inext ,&
    info, ipe, ips, ird, iwr, j, jj, mncl, mxlkw, mxlr ,&
    mxnca, mxnkw, ncat, ncc, nclass, nerr, nextl, nkwd ,&
    nstmts, ntcat, ntkwd, numr, numrr
  !
  INTEGER LU5, LU6, LU12, LU13, LU14, LU15, LU17, LU18, LU19
  PARAMETER (LU12=12,LU13=13,LU14=14,LU15=15,LU6=6,LU17=17,LU18=18,LU19=19,&
    LU5=5)
  CHARACTER(LLN) :: line
  CHARACTER(80) :: clline
  CHARACTER(80) :: msg
  CHARACTER(KMAXI) :: kwrds(KMAXJ), tkwd(MXNKWD)
  CHARACTER(MXLRN) :: rtname
  CHARACTER(MXLCAT) :: tcat(MXNCAT), etcat(MXNCAT), categ(15) ,&
    tclass(MXNCAT)
  INTEGER iptr(MXNCAT), jptr(MXNCAT), kptr(MXNCAT)
  CHARACTER(80) :: class(MXNCL+1), stmts(MXNCL)
  INTEGER iptrl(7*MXNRN), iptrr(7*MXNRN)
  !
  !     External functions.
  !
  INTEGER FIND, LENSTR
  LOGICAL IFDECK, IFIF, IFSID
  CHARACTER(10) :: CVTCAT
  EXTERNAL CVTCAT, FIND, IFDECK, IFIF, IFSID, LENSTR
  !
  !     Intrinsic functions.
  !
  INTRINSIC ABS, INDEX, MAX, MOD
  !***FIRST EXECUTABLE STATEMENT  SLPREP
  !     OPEN (UNIT=LU6, FILE=FOUT, STATUS='UNKNOWN', FORM='FORMATTED',
  !    +      IOSTAT = INFO)
  !     IF (INFO .NE. 0) THEN
  !       MSG = 'Failure in attempting to open ' // FOUT
  !       NERR = 1
  !       GO TO 240
  !     ENDIF
  !
  !     OPEN (UNIT=LU5, FILE=FINPUT, STATUS='UNKNOWN', FORM='FORMATTED',
  !    +      IOSTAT = INFO)
  !     IF (INFO .NE. 0) THEN
  !       MSG = 'Failure in attempting to open ' // FINPUT
  !       NERR = 1
  !       GO TO 240
  !     ENDIF
  !
  finp = ' '
  WRITE (UNIT=LU6,FMT=99015) DFINP
  READ (UNIT=LU5,FMT=99001,END=100) finp
  100 CONTINUE
  IF ( LENSTR(finp)==0 ) finp = DFINP
  OPEN (UNIT=LU15,FILE=finp,STATUS='OLD',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//finp
    nerr = 1
    GOTO 800
  ENDIF
  !
  fclass = ' '
  WRITE (UNIT=LU6,FMT=99013) DFCLAS
  READ (UNIT=LU5,FMT=99001,END=200) fclass
  200 CONTINUE
  IF ( LENSTR(fclass)==0 ) fclass = DFCLAS
  OPEN (UNIT=LU13,FILE=fclass,STATUS='OLD',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//fclass
    nerr = 1
    GOTO 800
  ENDIF
  !
  fcat = ' '
  WRITE (UNIT=LU6,FMT=99014) DFCAT
  READ (UNIT=LU5,FMT=99001,END=300) fcat
  300 CONTINUE
  IF ( LENSTR(fcat)==0 ) fcat = DFCAT
  OPEN (UNIT=LU14,FILE=fcat,STATUS='NEW',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//fcat
    nerr = 1
    GOTO 800
  ENDIF
  !
  fdaf = ' '
  WRITE (UNIT=LU6,FMT=99016) DFDAF
  READ (UNIT=LU5,FMT=99001,END=400) fdaf
  400 CONTINUE
  IF ( LENSTR(fdaf)==0 ) fdaf = DFDAF
  OPEN (UNIT=LU17,FILE=fdaf,STATUS='NEW',ACCESS='DIRECT',FORM='FORMATTED',&
    RECL=LLN,IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//fdaf
    nerr = 1
    GOTO 800
  ENDIF
  !
  ftbl = ' '
  WRITE (UNIT=LU6,FMT=99017) DFTBL
  READ (UNIT=LU5,FMT=99001,END=500) ftbl
  500 CONTINUE
  IF ( LENSTR(ftbl)==0 ) ftbl = DFTBL
  OPEN (UNIT=LU18,FILE=ftbl,STATUS='NEW',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//ftbl
    nerr = 1
    GOTO 800
  ENDIF
  !
  fkwd = ' '
  WRITE (UNIT=LU6,FMT=99018) DFKWD
  READ (UNIT=LU5,FMT=99001,END=600) fkwd
  600 CONTINUE
  IF ( LENSTR(fkwd)==0 ) fkwd = DFKWD
  OPEN (UNIT=LU19,FILE=fkwd,STATUS='NEW',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//fkwd
    nerr = 1
    GOTO 800
  ELSE
    ntkwd = 0
  ENDIF
  !
  flog = ' '
  WRITE (UNIT=LU6,FMT=99012) DFLOG
  READ (UNIT=LU5,FMT=99001,END=700) flog
  700 CONTINUE
  IF ( LENSTR(flog)==0 ) flog = DFLOG
  OPEN (UNIT=LU12,FILE=flog,STATUS='NEW',FORM='FORMATTED',IOSTAT=info)
  IF ( info/=0 ) THEN
    msg = 'Failure in attempting to open '//flog
    nerr = 1
    GOTO 800
  ENDIF
  !
  WRITE (UNIT=LU6,FMT=99011) flog
  !
  !     Write the names of all files to the transaction log file.
  !
  WRITE (UNIT=LU12,FMT=99010) finp, FOUT, fclass, fcat, fdaf, ftbl ,&
    fkwd
  !
  !     IRD is the "READ" line number and IWR is the "WRITE" line number.
  !
  ird = 0
  iwr = 0
  !
  !     NUMR is the number of routines in the library and MXLR is the
  !     maximum length of all these routines.
  !
  numr = 0
  numrr = 0
  mxlr = 0
  mxnca = 0
  ntcat = 0
  inext = MXNKWD + 1
  mxnkw = 0
  mxlkw = 0
  DO
    !
    READ (UNIT=LU15,FMT=99001,END=900) line
    ird = ird + 1
    IF ( .NOT.(IFDECK(line).OR.IFIF(line)) ) THEN
      iwr = iwr + 1
      WRITE (UNIT=LU17,FMT=99001,REC=iwr) line
      !
      IF ( line(1:6)=='      '.OR.IFSID(line) ) THEN
        !
        !       Subprogram statement.  Save starting record number.
        !
        ib = iwr
        ncat = 1
        categ(1) = ' '
        nkwd = 0
      ENDIF
      !
      IF ( line(2:6)=='***BE' ) THEN
        !
        !       BEGIN PROLOGUE record.  Extract routine name.
        !
        rtname = line(21:20+MXLRN)
        !
        !       Check case of routine name;  if not UPPER, terminate.
        !
        CALL UPCASE(rtname,rtname)
        IF ( rtname/=line(21:20+MXLRN) ) THEN
          msg = 'Routine name not in UPPER case'
          nerr = 2
          EXIT
        ENDIF
        numr = numr + 1
        IF ( numr>MXNRN ) THEN
          msg = 'Too many routine names.  Recompile code with larger MXNRN.'
          nerr = 2
          EXIT
        ENDIF
        !
        !       Write every 10th routine name to standard output.
        !
        IF ( MOD(numr,10)==0 ) WRITE (UNIT=LU6,FMT=99006) rtname, numr
      ENDIF
      !
      !     Effective with version 4.0, every routine MUST have a PURPOSE.
      !
      IF ( line(2:6)=='***PU' ) THEN
        !
        !       Record starting line number of Purpose section.
        !
        ips = iwr
        DO
          READ (UNIT=LU15,FMT=99001,END=900) line
          ird = ird + 1
          iwr = iwr + 1
          WRITE (UNIT=LU17,FMT=99001,REC=iwr) line
          IF ( line(2:4)=='***' ) THEN
            !
            !       Record ending line number of Purpose section.
            !
            ipe = iwr - 1
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !
      IF ( line(2:6)=='***CA' ) THEN
        !
        !       "CATEGORY" record.
        !
        ncat = 0
        !
        !       Initialize the pointer to point to the first category.
        !
        710        ic = 15
        ilen = LENSTR(line)
        IF ( ic>ilen ) THEN
          msg = 'No category on CATEGORY record'
          nerr = 3
          EXIT
        ENDIF
        DO
          !
          !       Get the offset location of the next delimiter.
          !
          icom = INDEX(line(ic:ilen),',')
          IF ( icom==0 ) THEN
            id = ilen + 1
          ELSE
            id = ic + icom - 1
          ENDIF
          ncat = ncat + 1
          IF ( ncat>15 ) THEN
            msg = 'Too many categories in a routine.  Recompile code with larger MXNCAT.'
            nerr = 3
            GOTO 800
          ENDIF
          !
          !       Put category into table in UPPER case.
          !
          categ(ncat) = ' '
          CALL UPCASE(line(ic:id-1),categ(ncat)(1:id-ic))
          IF ( line(ic:id-1)/=categ(ncat)(1:id-ic) ) THEN
            msg = 'Category not in UPPER case'
            nerr = 3
            GOTO 800
          ENDIF
          !
          !       Check to see if this category is already in the table of
          !       categories.
          !
          ifind = FIND(tcat,ntcat,categ(ncat))
          IF ( ifind<=0 ) THEN
            !
            !         Check to see if there is room in the table.
            !
            ntcat = ntcat + 1
            IF ( ntcat>MXNCAT ) THEN
              msg = 'Too many categories in the library.  Recompile code with larger MXNCAT.'
              nerr = 3
              GOTO 800
            ENDIF
            !
            !         Put category in table and check case.  If not UPPER case,
            !         terminate.
            !
            DO i = ntcat - 1, ABS(ifind) + 1, -1
              tcat(i+1) = tcat(i)
            ENDDO
            tcat(ABS(ifind)+1) = categ(ncat)
          ENDIF
          !
          !       Is this the last category for this subprogram?
          !
          IF ( icom==0 ) THEN
            mxnca = MAX(mxnca,ncat)
            EXIT
          ELSE
            !
            !       Is this the last category on the current line?
            !
            IF ( id>=ilen ) THEN
              !
              !         We are at the end of the line and need to read again.
              !
              READ (UNIT=LU15,FMT=99001,END=900) line
              ird = ird + 1
              IF ( line(1:5)/='C   '.AND.line(1:5)/='*   ' ) THEN
                msg = 'CATEGORY section not in correct form'
                nerr = 3
                GOTO 800
              ENDIF
              iwr = iwr + 1
              WRITE (UNIT=LU17,FMT=99001,REC=iwr) line
              GOTO 710
            ENDIF
            !
            !       Check to see that there is a blank character following the
            !       delimiter after the category just processed.
            !
            IF ( line(id+1:id+1)/=' ' ) THEN
              msg = 'CATEGORY section not in correct form'
              nerr = 3
              GOTO 800
            ENDIF
            !
            !       Set the pointer IC to the next category.
            !
            ic = id + 2
          ENDIF
        ENDDO
      ENDIF
      !
      IF ( line(2:6)=='***KE' ) THEN
        !
        !       "KEYWORD" section.
        !
        !
        !       Initialize the pointer to point to the first keyword phrase.
        !
        720        ic = 15
        ilen = LENSTR(line)
        IF ( ic>ilen ) THEN
          msg = 'No keyword phrase on KEYWORD record'
          nerr = 4
          EXIT
        ENDIF
        DO
          !
          !       Get the offset location of the next delimiter.
          !
          icom = INDEX(line(ic:ilen),',')
          IF ( icom==0 ) THEN
            id = ilen + 1
          ELSE
            id = ic + icom - 1
          ENDIF
          nkwd = nkwd + 1
          IF ( nkwd>KMAXJ ) THEN
            msg = 'Too many keyword phrases in a routine.  Recompile code with larger KMAXJ.'
            nerr = 4
            GOTO 800
          ENDIF
          kwrds(nkwd) = ' '
          IF ( id-ic>KMAXI ) THEN
            msg = 'Keyword phrase is too long'
            nerr = 4
            GOTO 800
          ENDIF
          !
          !       Put keyword phrase in table and check case.  If not UPPER case,
          !       terminate.
          !
          CALL UPCASE(line(ic:id-1),kwrds(nkwd)(1:id-ic))
          IF ( line(ic:id-1)/=kwrds(nkwd)(1:id-ic) ) THEN
            msg = 'Keyword phrase not in UPPER case'
            nerr = 4
            GOTO 800
          ENDIF
          mxlkw = MAX(mxlkw,id-ic)
          !
          !       Is this the last keyword for this subprogram?
          !
          IF ( icom==0 ) THEN
            mxnkw = MAX(mxnkw,nkwd)
            EXIT
          ELSE
            IF ( id>=ilen ) THEN
              !
              !         We are at the end of the line and need to read again.
              !
              READ (UNIT=LU15,FMT=99001,END=900) line
              ird = ird + 1
              IF ( line(1:5)/='C   '.AND.line(1:5)/='*   ' ) THEN
                msg = 'KEYWORD section not in correct form'
                nerr = 4
                GOTO 800
              ENDIF
              iwr = iwr + 1
              WRITE (UNIT=LU17,FMT=99001,REC=iwr) line
              GOTO 720
            ENDIF
            !
            !       Check to see that there is a blank character following the
            !       delimiter after the category just processed.
            !
            IF ( line(id+1:id+1)/=' ' ) THEN
              msg = 'KEYWORD section not in correct form'
              nerr = 4
              GOTO 800
            ENDIF
            ic = id + 2
          ENDIF
        ENDDO
      ENDIF
      !
      IF ( line(2:6)=='***EN' ) THEN
        !
        !       "END PROLOGUE" statement.
        !
        IF ( line(19:18+MXLRN)==rtname(1:MXLRN) ) THEN
          mxlr = MAX(mxlr,iwr-ib+1)
          !
          !         Write the categories for this routine, the routine name, etc.
          !         to the file FTBL.
          !
          WRITE (UNIT=LU18,FMT=99005) (categ(j),rtname,ib,iwr,ips,ipe,j=1,&
            ncat)
          numrr = numrr + ncat
          !
          IF ( nkwd>0 ) THEN
            mxnkw = MAX(mxnkw,nkwd)
            DO j = nkwd, 1, -1
              ifind = FIND(tkwd(1),ntkwd,kwrds(j))
              IF ( ifind<=0 ) THEN
                ntkwd = ntkwd + 1
                IF ( ntkwd>MXNKWD ) THEN
                  msg = 'Too many keyword phrases.  Recompile code with larger MXNKWD.'
                  nerr = 2
                  GOTO 800
                ENDIF
                DO i = ntkwd - 1, ABS(ifind) + 1, -1
                  tkwd(i+1) = tkwd(i)
                  iptrr(i+1) = iptrr(i)
                  iptrl(i+1) = iptrl(i)
                ENDDO
                tkwd(ABS(ifind)+1) = kwrds(j)
                iptrr(ABS(ifind)+1) = numrr
                iptrl(ABS(ifind)+1) = 0
              ELSE
                ientry = iptrl(ifind)
                IF ( ientry==0 ) THEN
                  !
                  !                 No entries in the chain.
                  !
                  iptrl(ifind) = inext
                ELSE
                  DO
                    nextl = iptrl(ientry)
                    IF ( nextl/=0 ) THEN
                      !
                      !                   Not at end of chain.
                      !
                      ientry = nextl
                      CYCLE
                    ENDIF
                    iptrl(ientry) = inext
                    EXIT
                  ENDDO
                ENDIF
                iptrl(inext) = 0
                iptrr(inext) = numrr
                inext = inext + 1
              ENDIF
            ENDDO
          ENDIF
        ELSE
          msg = 'BEGIN name and END name differ'
          nerr = 5
          EXIT
        ENDIF
        DO
          !
          !       Read to the next *DECK record.
          !
          READ (UNIT=LU15,FMT=99001,END=900) line
          ird = ird + 1
          IF ( IFDECK(line) ) EXIT
        ENDDO
        !
        !     GOTO read.
        !
      ENDIF
    ENDIF
  ENDDO
  !
  !     Abnormal termination.
  !
  800  CALL XERMSG(' ','SLPREP',msg,nerr,1)
  WRITE (UNIT=LU6,FMT=99009) rtname, ird, line
  !
  !     Write the job abort message to the transaction log file also.
  !
  WRITE (UNIT=LU12,FMT=99009) rtname, ird, line
  CLOSE (UNIT=LU14,STATUS='DELETE')
  CLOSE (UNIT=LU17,STATUS='DELETE')
  CLOSE (UNIT=LU18,STATUS='DELETE')
  CLOSE (UNIT=LU19,STATUS='DELETE')
  GOTO 1100
  !
  900 CONTINUE
  DO j = 1, ntcat
    etcat(j) = CVTCAT(tcat(j))
  ENDDO
  CALL SORT(etcat,ntcat,1,tcat)
  !
  !     Read in the GAMS classification file.
  !
  i = 0
  DO
    READ (UNIT=LU13,FMT=99001,END=1000) clline
    IF ( i>=MXNCL ) THEN
      msg = 'Too many lines in file '//fclass
      nerr = 6
      GOTO 800
    ENDIF
    i = i + 1
    class(i) = clline
  ENDDO
  1000 nclass = i
  class(i+1) = '@'
  mncl = MXNCL
  CALL PSCAT(etcat,ntcat,class,mncl,ncc,tclass,iptr,jptr,kptr,nstmts,stmts,&
    nerr)
  IF ( nerr/=0 ) THEN
    msg = 'Could not locate '//tclass(nerr)//' in GAMS file.'
    line = 'Input file completely read.'
    nerr = 7
    GOTO 800
  ENDIF
  !
  !     The reason NCC+1 is being written is to call attention to
  !     the fact that an extra line is included in the KPTR array
  !     which gives the line beyond the last line of the STMTS array.
  !
  WRITE (UNIT=LU14,FMT=99002) ncc + 1
  WRITE (UNIT=LU14,FMT=99003) (iptr(i),jptr(i),kptr(i),tclass(i)(1:LENSTR(&
    tclass(i))),i=1,ncc)
  WRITE (UNIT=LU14,FMT=99004) kptr(ncc+1)
  WRITE (UNIT=LU14,FMT=99001) (stmts(i)(1:LENSTR(stmts(i))),i=1,nstmts)
  !
  WRITE (UNIT=LU19,FMT=99002) ntkwd
  DO j = 1, ntkwd
    ilen = LENSTR(tkwd(j))
    WRITE (UNIT=LU19,FMT=99001) tkwd(j)(1:ilen)
  ENDDO
  !
  !     We now need to compress the IPTRL and IPTRR tables and remove
  !     any unneeded cells between the original allocated cells (MXNKWD)
  !     and the final number (NTKWD) of pointer cells.
  !
  ichng = MXNKWD - ntkwd
  IF ( ichng>0 ) THEN
    DO j = 1, ntkwd
      IF ( iptrl(j)/=0 ) iptrl(j) = iptrl(j) - ichng
    ENDDO
    jj = ntkwd + 1
    DO j = MXNKWD + 1, inext - 1
      IF ( iptrl(j)/=0 ) THEN
        iptrl(jj) = iptrl(j) - ichng
      ELSE
        !           IPTRL(JJ)=IPTRL(J)
        iptrl(jj) = 0
      ENDIF
      iptrr(jj) = iptrr(j)
      jj = jj + 1
    ENDDO
  ENDIF
  !
  !     Set last set of pointers to 0.
  !
  iptrl(inext-ichng) = 0
  iptrr(inext-ichng) = 0
  DO j = 1, inext - ichng
    WRITE (UNIT=LU19,FMT=99002) iptrl(j), iptrr(j)
  ENDDO
  !
  !     Normal termination.
  !
  WRITE (UNIT=LU6,FMT=99007) numr, ird, iwr, ntcat, ntkwd, mxlr ,&
    mxnca, mxnkw, mxlkw, numrr, ntkwd ,&
    inext - ichng
  WRITE (UNIT=LU6,FMT=99008) nclass, ncc, nstmts
  !
  !     Write summary information to the transaction log file also.
  !
  WRITE (UNIT=LU12,FMT=99007) numr, ird, iwr, ntcat, ntkwd, mxlr ,&
    mxnca, mxnkw, mxlkw, numrr, ntkwd ,&
    inext - ichng
  WRITE (UNIT=LU12,FMT=99008) nclass, ncc, nstmts
  CLOSE (UNIT=LU14)
  !
  !     Call to FTRUN for CTSS/CFTLIB systems only.
  !
  CLOSE (UNIT=LU17)
  CLOSE (UNIT=LU18)
  CLOSE (UNIT=LU19)
  !
  1100 CLOSE (UNIT=LU12)
  CLOSE (UNIT=LU13)
  CLOSE (UNIT=LU15)
  !     CLOSE (UNIT=LU5)
  !     CLOSE (UNIT=LU6)
  STOP
  !
  99001 FORMAT (A)
  99002 FORMAT (I5,2X,I5)
  99003 FORMAT (3I5,3X,A)
  99004 FORMAT (I15)
  99005 FORMAT (1X,2A,4I8)
  99006 FORMAT (' Processing routine ',A,',  number ',I4)
  99007 FORMAT (//' S U M M A R Y There were ',I4,&
    ' documentation modules and a total ','of ',I6,&
    ' lines read '/' from the sequential input file and a total of ',&
    I6,' lines written '/' to the direct access documentation file.'/&
    ' The library has a total of ',I3,' distinct ','categories and ',&
    I3,' distinct keywords.'/' At least one module has ',I4,&
    ' lines, one module ','has ',I2,' categories,'/' one module has ',&
    I2,' keyword phrases and ','the longest keyword phrase is '/1X,I2,&
    ' characters.'/' There were ',I4,&
    ' category/routine lines written to ',&
    'the category file.'/' There were ',I3,&
    ' distinct keyword phrases written ',&
    'to the keyword file, and '/1X,I4,&
    ' keyword pointers written to the file.')
  99008 FORMAT (' There were ',I3,' lines read from the GAMS ',&
    'classification scheme file.'/' There were ',I3,&
    ' classification categories written ',&
    'to the classification file'/' and a total of ',I3,&
    ' classification descriptors ','written to the file.'//)
  99009 FORMAT (//' J O B   A B O R T The job aborted in routine ',A,&
    ' after processing ','approximately ',I6,&
    ' lines'/' of the input file.'/&
    ' The last line read from the input file was '/1X,A//)
  99010 FORMAT (//5X,'T R A N S A C T I O N   L O G'//2X,'Input source file: ',&
    A/2X,'Output source file: ',A/2X,'GAMS classification file: ',&
    A/2X,'Sequential file containing linked list of all ',&
    'classifications used: ',A/2X,&
    'Direct access file containing the documentation: ',A/2X,&
    'Sequential file containing the deck names and ',&
    'prologue line numbers: ',A/2X,&
    'Sequential file of keyword phrases and pointer ','arrays: ',A/)
  99011 FORMAT (2X,'Transaction log file: ',A/)
  99012 FORMAT (//,' Enter the name of your transaction log file or',' <cr>',/,&
    ' (The default is ''',A,''')')
  99013 FORMAT (//,' Enter the name of the GAMS (SLATEC) classification',&
    ' file or <cr>',/,' (The default is ''',A,''')')
  99014 FORMAT (//,' Enter the name of the output file to contain the',&
    ' classifications used ',/,' or <cr>',/,' (The default is ''',A,''')')
  99015 FORMAT (//,' Enter the name of your input source file or <cr>',/,&
    ' (The default is ''',A,''')')
  99016 FORMAT (//,' Enter the name of the output direct access file',&
    ' which will contain the ',/,' documentation or',' <cr>',/,&
    ' (The default is ''',A,''')')
  99017 FORMAT (//,' Enter the name of the output file which will',&
    ' contain the line numbers',/,' of data in the',&
    ' documentation file or <cr>',/,' (The default is ''',A,''')')
  99018 FORMAT (//,' Enter the name of the output file which will',&
    ' contain the keyword phrases ',/,' and their',&
    ' pointers or <cr>',/,' (The default is ''',A,''')')
END PROGRAM SLPREP
!DECK IFDECK
LOGICAL FUNCTION IFDECK(Line)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  IFDECK
  !***PURPOSE  Determine if a line is an *DECK or *DK statement.
  !***LIBRARY   SLATEC
  !***CATEGORY  R3C
  !***TYPE      LOGICAL (IFDECK-L)
  !***KEYWORDS  DECK
  !***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   This FUNCTION determines if a record is an  *DECK or *DK  statement.
  !   If it is, the value  .TRUE.  is returned; otherwise, the value
  !   .FALSE. is returned.
  !
  !   INPUT
  !
  !      LINE   - a character variable containing the string to be
  !               examined.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  UPCASE
  !***REVISION HISTORY  (YYMMDD)
  !   840427  DATE WRITTEN
  !   840430  REVISION DATE from the pre-1990 prologue.
  !   891215  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  IFDECK
  !     .. Scalar Arguments ..
  CHARACTER*(*) Line
  !     .. Local Scalars ..
  CHARACTER(6) :: temp
  !     .. External Subroutines ..
  EXTERNAL UPCASE
  !***FIRST EXECUTABLE STATEMENT  IFDECK
  CALL UPCASE(Line(1:6),temp(1:6))
  IFDECK = .TRUE.
  IF ( temp(1:6)/='*DECK '.AND.temp(1:4)/='*DK ' ) IFDECK = .FALSE.
END FUNCTION IFDECK
!DECK IFIF
LOGICAL FUNCTION IFIF(Line)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  IFIF
  !***PURPOSE  Determine if a line is an  *IF  statement.
  !***LIBRARY   (NONE)
  !***CATEGORY  (NONE)
  !***KEYWORDS  (NONE)
  !***AUTHOR  Boland, W. Robert, C-8, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   This FUNCTION determines if a record is an  *IF  statement.
  !   If it is, the value  .TRUE.  is returned; otherwise, the value
  !   .FALSE. is returned.
  !
  !   INPUT
  !
  !      LINE   - a character variable containing the string to be
  !               examined.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  UPCASE
  !***REVISION HISTORY  (YYMMDD)
  !   910208  DATE WRITTEN
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  IFIF
  !     .. Scalar Arguments ..
  CHARACTER*(*) Line
  !     .. Local Scalars ..
  CHARACTER(3) :: temp
  !     .. External Subroutines ..
  EXTERNAL UPCASE
  !***FIRST EXECUTABLE STATEMENT  IFIF
  CALL UPCASE(Line(1:3),temp(1:3))
  IFIF = .TRUE.
  IF ( temp(1:3)/='*IF' ) IFIF = .FALSE.
END FUNCTION IFIF
!DECK IFSID
LOGICAL FUNCTION IFSID(Line)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  IFSID
  !***PURPOSE  Determine if a line is a special *...IDENT statement.
  !***LIBRARY   (NONE)
  !***CATEGORY  (NONE)
  !***KEYWORDS  (NONE)
  !***AUTHOR  Bacon, Barbara, C-10, (LANL)
  !***DESCRIPTION
  !
  !   This FUNCTION determines if a record is a  *...IDENT  statement.
  !   If it is, the value  .TRUE.  is returned; otherwise, the value
  !   .FALSE. is returned.
  !
  !   INPUT
  !
  !     LINE   - a character variable containing the string to be
  !              examined.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  UPCASE
  !***REVISION HISTORY  (YYMMDD)
  !   910211  DATE WRITTEN
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  IFSID
  !     .. Scalar Arguments ..
  CHARACTER*(*) Line
  !     .. Local Scalars ..
  INTEGER i, id
  CHARACTER(20) :: temp
  !     .. External Subroutines ..
  EXTERNAL UPCASE
  !     .. Intrinsic Functions ..
  INTRINSIC INDEX
  !***FIRST EXECUTABLE STATEMENT  IFSID
  CALL UPCASE(Line(1:20),temp(1:20))
  IFSID = .FALSE.
  IF ( Line(1:1)=='*' ) THEN
    id = INDEX(Line,'IDENT')
    IF ( id/=0 ) THEN
      !         IF (LINE(2:ID-1) .EQ. ' ') IFSID = .TRUE.
      DO i = 2, id - 1
        IF ( Line(i:i)/=' ' ) RETURN
      ENDDO
      !
      !         We have a line with *^^...^^IDENT
      !
      IFSID = .TRUE.
    ENDIF
  ENDIF
  RETURN
END FUNCTION IFSID
!DECK PSCAT
SUBROUTINE PSCAT(Ecat,Ncat,Class,Mncl,Ncc,Tclass,Iptr,Jptr,Kptr,Istmt,&
    Stmts,Nerr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PSCAT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SLPREP
  !***LIBRARY   (NONE)
  !***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !     Creates the data for the FCAT file consisting of two parts:
  !     1) linked lists enabling one to obtain the various categories
  !        and subcategories present in the library being examined
  !     2) statements from the classification file associated with
  !        the categories.
  !
  !   INPUT
  !
  !     ECAT   - an array of "extended" categories found in the library
  !
  !     NCAT   - the number of categories in the above arrays
  !
  !   OUTPUT
  !
  !     NCC    - the number of elements in TCLASS, IPTR, etc.
  !
  !     TCLASS - the expanded table of categories for this library;
  !              e.g., if A1C were a category in the library, then
  !              TCLASS would contain A, A1, and A1C.
  !
  !     IPTR   - the linked list allowing one to travel down categories
  !              of the same length (i.e., A to C to ... to Z;
  !              or  A1 to A3 to ... to An).
  !
  !     JPTR   - the linked list allowing one to get the next larger
  !              element in a category (i.e., A to A1)
  !
  !     KPTR   - the location in the STMTS array of the statement
  !              associated with the element in TCLASS.
  !
  !     ISTMT  - the number of statements in STMTS
  !
  !     STMTS  - the statements from the classification file associated
  !              with TCLASS.
  !
  !     NERR   - the number of the element in TCLASS which could not
  !              be located in the GAMS file.  This value is set to
  !              zero upon entry to PSCAT and will indicate an error
  !              to the calling routine if non-zero upon exit.
  !
  !     The easiest way to describe how this program works is to
  !     use examples.
  !     1) A02B05
  !        At the beginning of the DO loop, the following values
  !        will be set:
  !        OPART1=OPART2=OPART3=OPART4=OPART5=OPART6=OPART7=' '
  !        PART1='A'   PART2='02'   PART3='B'   PART4='05'
  !        PART5=' '   PART6='  '   PART7=' '
  !
  !        Since PART1 .NE. OPART1, we will go to 100, clear all
  !        the OPART's.  Then since OPART1 is not equal to PART1,
  !        OPART1 will be set to PART1 and the first member of
  !        TCLASS will be set to 'A'.  Then we will loop back to
  !        statement 160 to build the rest of A02B05.
  !        Now OPART1 is equal to PART1, so go to the next IF block.
  !        Since OPART2=' ' and PART2='02', we will set OPART2=PART2,
  !        and build the next member of TCLASS to be "A02" by
  !        concatenating OPART1 and OPART2.
  !        Then loop back to statement 160.
  !        This process continues down through PART4, so that we will
  !        build A02B and A02B05 in TCLASS.  So on the fifth time
  !        down we will have OPART1='A', OPART2='02', OPART3='B',
  !        OPART4='05', OPART5=OPART6=OPART7=' '.
  !        Then since PART5=' ', we will transfer to the end of the
  !        do 180 loop and return to the loop's beginning to pick up
  !        a new member of ECAT.
  !     2) A02B05C
  !        At this point, the following values will be set when
  !        we reach the first IF statement:
  !        OPART1='A'  OPART2='02'  OPART3='B'  OPART4='05'
  !        OPART5=' '  OPART6='  '  OPART7=' '
  !        PART1 ='A'  PART2 ='02'  PART3 ='B'  PART4 ='05'
  !        PART5 ='C'  PART6 =' '   PART7 =' '
  !        So we will transfer to 100 and be sure that OPART5, OPART6,
  !        and OPART7 are ' ' (it is possible for those to be non-blank
  !        if the previous value of ECAT had been longer).
  !        Then we will bypass the first four IF blocks and do number
  !        five.  Since PART5 is not blank and OPART5.NE.PART5, we
  !        will set OPART5=PART5 which is 'C' and build the next
  !        member of TCLASS as "A02B05C' by concatenating OPART1,
  !        OPART2, OPART3, OPART4, and OPART5.
  !        Loop back to 160 where we will completely fall through all
  !        the IF blocks to the end of the loop.
  !     3) A03
  !        Now PART1='A'  PART2='03' and PART3 through PART7 will be ' '.
  !        But OPART1='A', OPART2='02', OPART3='B', OPART4='05',
  !        OPART5='C' and OPART6=OPART7=' '.
  !        This time we will branch to statement 40 since OPART1=PART1,
  !        but OPART2 is not equal to PART2.
  !        OPART2 through OPART7 will be set to ' '.
  !        We will bypass the first IF block and do the second one,
  !        since PART2 is not blank and OPART2.NE.PART2.
  !        OPART2 will be set to PART2 (which is '03') and we will
  !        build the next member of TCLASS by concatenating OPART1 and
  !        OPART2 ('A03').  Then go back to statement 160.  All the IF
  !        blocks will fail and we will drop through to the end of the
  !        loop and go back to the beginning of the loop to get a new
  !        value of ECAT.
  !     4) This process continues until all the members of ECAT have been
  !        processed.  The complete table is built in TCLASS and the
  !        number of elements in it is NCC.
  !
  !
  !     The next two DO blocks (the DO 280 and the DO 320) build the
  !     IPTR and JPTR arrays, such that the members of IPTR point to
  !     the next member of TCLASS which has the same length of "CLASS"
  !     and the (n-1)th part of the member is the same but the n-th part
  !     is different.  For example, the IPTR pointer associated with 'A'
  !     would point to where 'C' is in TCLASS; then the IPTR pointer
  !     associated with 'C' would point to where 'D' is, etc.  The final
  !     IPTR pointer associated with a single letter in TCLASS will be
  !     zero.  The IPTR pointer for A04 would point to where A06 is in
  !     TCLASS and the IPTR pointer for A06 would be zero because there
  !     are no more An.  The IPTR pointer for C01 points to where C02
  !     is in TCLASS.  The IPTR pointer for C02 points to where C03 is.
  !     ...  C14 is the last member of this type and its IPTR value is
  !     zero.
  !
  !     The JPTR pointers are similar to the IPTR pointers, but point to
  !     where the next sized member of TCLASS will be located.  For
  !     example, the JPTR pointer for A points to where A04 is in TCLASS,
  !     the JPTR pointer for A04 points to where A04A is in TCLASS, and
  !     the JPTR pointer for A04A is zero since there is nothing more in
  !     TCLASS which starts with A04A.
  !
  !     The final section puts into the array STMTS the statements
  !     associated with each member of TCLASS.  Since a statement may take
  !     up more than one line, we need another pointer table (KPTR) which
  !     will point to the member of STMTS associated with each member of
  !     TCLASS.  E.g. The statement for A is in STMTS(1), ..., the
  !     statement for C01 is in STMTS(7), but the statement for C02 is
  !     STMTS(9) because the one for C01 is two lines long.
  !     By setting up KPTR this way (and putting a dummy value in the
  !     cell beyond the last entry corresponding to a member of TCLASS),
  !     we can find out how long a statement is by taking the KPTR pointer
  !     for the member of TCLASS in which we are interested and subtract
  !     it from the next member of KPTR.
  !
  !     The following is an excerpt from a file produced by this
  !     code to show how it works.  Column 1 are the IPTR pointers,
  !     column 2 are the JPTR pointers, column 3 are the KPTR pointers,
  !     and column 4 are the members of TCLASS.
  !
  !       |<-------6       <----2       1   A
  !       | |<-----4       |<-->3       2   A04
  !       | |      0       |--->0       3   A04A
  !       | |----->0       <----5       4   A06
  !       |        0       |--->0       5   A06B
  !       |<----->37       <----7       6   C
  !       | |<-----8       |--->0       7   C01
  !       | |<---->9            0       9   C02
  !       | |<--->12       <---10      10   C03
  !       | |      0       |<->11      11   C03A
  !       | |      0       |--->0      12   C03A02
  !       | |<--->16       <---13      13   C04
  !       | | |--<14       |--->0      14   C04A
  !       | | |<->15            0      15   C04B
  !       | | |--->0            0      16   C04C
  !       | |<--->17            0      17   C05
  !       | |<--->23       <---18      18   C07
  !       | | |<--19       |--->0      19   C07A
  !       | | |<->20            0      20   C07B
  !       | | |<->21            0      21   C07C
  !       | | |<->22            0      22   C07E
  !       | | |--->0            0      23   C07F
  !       | |<--->26       <---24      24   C08
  !       | | |--<25       |--->0      25   C08A
  !       | | |--->0            0      27   C08C
  !       | |<--->35       <---27      28   C10
  !       | | |--<30       |<->28      29   C10A
  !       | | | |<29       |--->0      30   C10A01
  !       | | | |->0            0      31   C10A03
  !       | | |<->33       <---31      32   C10B
  !       | | | |<32       |--->0      33   C10B01
  !       | | | |->0            0      34   C10B03
  !       | | |<->34            0      35   C10D
  !       | | |--->0            0      36   C10F
  !       | |<--->36            0      37   C11
  !       | |----->0            0      38   C14
  !       |<---->123       <---38      39   D
  !
  !***ROUTINES CALLED  CVTCAT, LENSTR
  !***REVISION HISTORY  (YYMMDD)
  !   891215  DATE WRITTEN
  !   891215  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  PSCAT
  !     .. Scalar Arguments ..
  INTEGER Istmt, Mncl, Ncat, Ncc, Nerr
  !     .. Array Arguments ..
  INTEGER Iptr(*), Jptr(*), Kptr(*)
  CHARACTER*(*) Class(*), Ecat(*), Stmts(*), Tclass(*)
  !     .. Local Scalars ..
  INTEGER i, iclass, ilen, iper, istart, j, k, nlen
  CHARACTER :: opart1, opart3, opart5, opart7, part1, part3, part5 ,&
    part7
  CHARACTER(2) :: opart2, opart4, opart6, part2, part4, part6
  !     .. Local Arrays ..
  INTEGER size(7)
  !     .. External Functions ..
  INTEGER LENSTR
  CHARACTER(10) :: CVTCAT
  EXTERNAL LENSTR, CVTCAT
  !     .. Intrinsic Functions ..
  INTRINSIC INDEX
  !     .. Data statements ..
  DATA size/1, 3, 4, 6, 7, 9, 10/
  !***FIRST EXECUTABLE STATEMENT  PSCAT
  Nerr = 0
  Ncc = 0
  opart1 = ' '
  opart2 = ' '
  opart3 = ' '
  opart4 = ' '
  opart5 = ' '
  opart6 = ' '
  opart7 = ' '
  DO i = 1, Ncat
    part1 = Ecat(i)(1:1)
    part2 = Ecat(i)(2:3)
    part3 = Ecat(i)(4:4)
    part4 = Ecat(i)(5:6)
    part5 = Ecat(i)(7:7)
    part6 = Ecat(i)(8:9)
    part7 = Ecat(i)(10:10)
    IF ( part1/=opart1 ) THEN
      opart1 = ' '
      opart2 = ' '
      opart3 = ' '
      opart4 = ' '
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part2/=opart2 ) THEN
      opart2 = ' '
      opart3 = ' '
      opart4 = ' '
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part3/=opart3 ) THEN
      opart3 = ' '
      opart4 = ' '
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part4/=opart4 ) THEN
      opart4 = ' '
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part5/=opart5 ) THEN
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part6/=opart6 ) THEN
      opart6 = ' '
      opart7 = ' '
    ELSEIF ( part7/=opart7 ) THEN
      opart7 = ' '
    ELSE
      opart1 = ' '
      opart2 = ' '
      opart3 = ' '
      opart4 = ' '
      opart5 = ' '
      opart6 = ' '
      opart7 = ' '
    ENDIF
    DO
      IF ( opart1/=part1 ) THEN
        opart1 = part1
        Ncc = Ncc + 1
        Tclass(Ncc) = opart1
      ENDIF
      IF ( part2/=' ' ) THEN
        IF ( opart2/=part2 ) THEN
          opart2 = part2
          Ncc = Ncc + 1
          Tclass(Ncc) = opart1//opart2
          CYCLE
        ENDIF
        IF ( part3/=' ' ) THEN
          IF ( opart3/=part3 ) THEN
            opart3 = part3
            Ncc = Ncc + 1
            Tclass(Ncc) = opart1//opart2//opart3
            CYCLE
          ENDIF
          IF ( part4/=' ' ) THEN
            IF ( opart4/=part4 ) THEN
              opart4 = part4
              Ncc = Ncc + 1
              Tclass(Ncc) = opart1//opart2//opart3//opart4
              CYCLE
            ENDIF
            IF ( part5/=' ' ) THEN
              IF ( opart5/=part5 ) THEN
                opart5 = part5
                Ncc = Ncc + 1
                Tclass(Ncc) = opart1//opart2//opart3//opart4//opart5
                CYCLE
              ENDIF
              IF ( part6/=' ' ) THEN
                IF ( opart6/=part6 ) THEN
                  opart6 = part6
                  Ncc = Ncc + 1
                  Tclass(Ncc) = opart1//opart2//opart3//opart4//opart5//&
                    opart6
                  CYCLE
                ENDIF
                IF ( part7/=' ' ) THEN
                  IF ( opart7/=part7 ) THEN
                    opart7 = part7
                    Ncc = Ncc + 1
                    Tclass(Ncc) = opart1//opart2//opart3//opart4//opart5//&
                      opart6//opart7
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      EXIT
    ENDDO
  ENDDO
  DO i = 1, Ncc
    Iptr(i) = 0
    Jptr(i) = 0
  ENDDO
  DO j = 1, 7
    istart = 1
    DO
      DO i = istart, Ncc
        ilen = LENSTR(Tclass(i))
        IF ( ilen==size(j) ) THEN
          DO k = i + 1, Ncc
            nlen = LENSTR(Tclass(k))
            IF ( ilen==nlen ) THEN
              IF ( j==1 ) THEN
                Iptr(i) = k
              ELSEIF ( Tclass(i)(1:size(j-1))==Tclass(k)(1:size(j-1)) ) THEN
                Iptr(i) = k
              ELSE
                Iptr(i) = 0
              ENDIF
              istart = k
              GOTO 50
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      EXIT
      50 CONTINUE
    ENDDO
  ENDDO
  DO j = 1, 7
    DO i = 1, Ncc
      ilen = LENSTR(Tclass(i))
      IF ( ilen==size(j) ) THEN
        nlen = LENSTR(Tclass(i+1))
        IF ( ilen<nlen ) THEN
          IF ( Tclass(i)(1:size(j))==Tclass(i+1)(1:size(j)) ) Jptr(i)&
            = i + 1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO i = 1, Ncc
    Kptr(i) = 0
  ENDDO
  !
  iclass = 1
  Istmt = 1
  DO i = 1, Ncc
    DO
      !
      iper = INDEX(Class(iclass)(1:8),'.')
      IF ( iper==0 ) THEN
        iclass = iclass + 1
        IF ( iclass>Mncl ) THEN
          Nerr = i
          RETURN
        ENDIF
        CYCLE
      ENDIF
      !
      IF ( Tclass(i)==CVTCAT(Class(iclass)(1:iper-1)) ) THEN
        Kptr(i) = Istmt
        DO
          Stmts(Istmt) = Class(iclass)(iper+2:80)
          Istmt = Istmt + 1
          iclass = iclass + 1
          IF ( Class(iclass)(1:1)/=' ' ) GOTO 100
        ENDDO
      ELSE
        iclass = iclass + 1
        IF ( iclass>Mncl ) THEN
          Nerr = i
          RETURN
        ENDIF
      ENDIF
    ENDDO
    100 CONTINUE
  ENDDO
  Kptr(Ncc+1) = Istmt
END SUBROUTINE PSCAT
!DECK SORT
SUBROUTINE SORT(R,N,Nr,Cr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SORT
  !***PURPOSE  Alphabetise a character array and re-order one dependent
  !            character array.
  !***LIBRARY   (NONE)
  !***CATEGORY  (NONE)
  !***KEYWORDS  (NONE)
  !***AUTHOR  Bacon, Barbara A., C-10, Los Alamos National Laboratory
  !***DESCRIPTION
  !
  !   This subroutine alphabetises a character array.  It then accordingly
  !   reorders up to one dependent character array.
  !
  !   INPUT
  !
  !      R      - the character array to be alphabetised.
  !      N      - the number of elements in the array R.
  !      NR     - the number (0 or 1) of dependent arrays.
  !      CR     - the dependent array of length at most 15 characters.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891101  DATE WRITTEN
  !   891215  Prologue converted to Version 4.0 format.  (BAB)
  !   920911  Declarations section restructured.  (WRB)
  !***END PROLOGUE  SORT
  !     .. Scalar Arguments ..
  INTEGER N, Nr
  !     .. Array Arguments ..
  CHARACTER*(*) Cr(*), R(*)
  !     .. Local Scalars ..
  INTEGER i, j, k, m
  CHARACTER(15) :: it
  !     .. Intrinsic Functions ..
  INTRINSIC MOD
  !***FIRST EXECUTABLE STATEMENT  SORT
  m = N
  !     DO 5 I = 1,N
  !       IR(I) = I
  !   5 CONTINUE
  100 CONTINUE
  IF ( m>15 ) THEN
    m = m/4
  ELSE
    m = m/2
    IF ( m==0 ) RETURN
  ENDIF
  IF ( MOD(m,2)==0 ) m = m + 1
  k = N - m
  j = 1
  i = j
  200 CONTINUE
  DO WHILE ( R(i)>R(i+m) )
    it = R(i)
    R(i) = R(i+m)
    R(i+m) = it
    IF ( Nr>0 ) THEN
      it = Cr(i)
      Cr(i) = Cr(i+m)
      Cr(i+m) = it
    ENDIF
    i = i - m
    IF ( i<1 ) EXIT
  ENDDO
  j = j + 1
  IF ( j>k ) GOTO 100
  i = j
  GOTO 200
END SUBROUTINE SORT
