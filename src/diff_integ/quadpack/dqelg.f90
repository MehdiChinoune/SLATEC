!** DQELG
SUBROUTINE DQELG(N,Epstab,Result,Abserr,Res3la,Nres)
  !> The routine determines the limit of a given sequence of
  !            approximations, by means of the Epsilon algorithm of
  !            P.Wynn. An estimate of the absolute error is also given.
  !            The condensed Epsilon table is computed. Only those
  !            elements needed for the computation of the next diagonal
  !            are preserved.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (QELG-S, DQELG-D)
  !***
  ! **Keywords:**  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **Description:**
  !
  !           Epsilon algorithm
  !           Standard fortran subroutine
  !           Double precision version
  !
  !           PARAMETERS
  !              N      - Integer
  !                       EPSTAB(N) contains the new element in the
  !                       first column of the epsilon table.
  !
  !              EPSTAB - Double precision
  !                       Vector of dimension 52 containing the elements
  !                       of the two lower diagonals of the triangular
  !                       epsilon table. The elements are numbered
  !                       starting at the right-hand corner of the
  !                       triangle.
  !
  !              RESULT - Double precision
  !                       Resulting approximation to the integral
  !
  !              ABSERR - Double precision
  !                       Estimate of the absolute error computed from
  !                       RESULT and the 3 previous results
  !
  !              RES3LA - Double precision
  !                       Vector of dimension 3 containing the last 3
  !                       results
  !
  !              NRES   - Integer
  !                       Number of calls to the routine
  !                       (should be zero at first call)
  !
  !***
  ! **See also:**  DQAGIE, DQAGOE, DQAGPE, DQAGSE
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : D1MACH
  !
  REAL(DP) :: Abserr, delta1, delta2, delta3, epmach, &
    epsinf, Epstab(52), error, err1, err2, err3, e0, e1, &
    e1abs, e2, e3, oflow, res, Result, Res3la(3), ss, tol1, tol2, tol3
  INTEGER :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, N, newelm, Nres, num
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
  !           E1       ELEMENT IN THE EPSILON TABLE IS BASED
  !           E2
  !           E3                 E0
  !                        E3    E1    NEW
  !                              E2
  !           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
  !                    DIAGONAL
  !           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
  !           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
  !                    OF ERROR
  !
  !           MACHINE DEPENDENT CONSTANTS
  !           ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
  !           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
  !           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
  !           DIAGONAL OF THE EPSILON TABLE IS DELETED.
  !
  !* FIRST EXECUTABLE STATEMENT  DQELG
  epmach = D1MACH(4)
  oflow = D1MACH(2)
  Nres = Nres + 1
  Abserr = oflow
  Result = Epstab(N)
  IF( N>=3 ) THEN
    limexp = 50
    Epstab(N+2) = Epstab(N)
    newelm = (N-1)/2
    Epstab(N) = oflow
    num = N
    k1 = N
    DO i = 1, newelm
      k2 = k1 - 1
      k3 = k1 - 2
      res = Epstab(k1+2)
      e0 = Epstab(k3)
      e1 = Epstab(k2)
      e2 = res
      e1abs = ABS(e1)
      delta2 = e2 - e1
      err2 = ABS(delta2)
      tol2 = MAX(ABS(e2),e1abs)*epmach
      delta3 = e1 - e0
      err3 = ABS(delta3)
      tol3 = MAX(e1abs,ABS(e0))*epmach
      IF( err2>tol2 .OR. err3>tol3 ) THEN
        e3 = Epstab(k1)
        Epstab(k1) = e1
        delta1 = e1 - e3
        err1 = ABS(delta1)
        tol1 = MAX(e1abs,ABS(e3))*epmach
        !
        !           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
        !           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
        !
        IF( err1>tol1 .AND. err2>tol2 .AND. err3>tol3 ) THEN
          ss = 1._DP/delta1 + 1._DP/delta2 - 1._DP/delta3
          epsinf = ABS(ss*e1)
          !
          !           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
          !           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
          !           OF N.
          !
          IF( epsinf>0.1D-03 ) THEN
            !
            !           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
            !           THE VALUE OF RESULT.
            !
            res = e1 + 1._DP/ss
            Epstab(k1) = res
            k1 = k1 - 2
            error = err2 + ABS(res-e2) + err3
            IF( error<=Abserr ) THEN
              Abserr = error
              Result = res
            END IF
            CYCLE
          END IF
        END IF
        N = i + i - 1
        !- **JUMP OUT OF DO-LOOP
        EXIT
      ELSE
        !
        !           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
        !           ACCURACY, CONVERGENCE IS ASSUMED.
        !           RESULT = E2
        !           ABSERR = ABS(E1-E0)+ABS(E2-E1)
        !
        Result = res
        Abserr = err2 + err3
        !- **JUMP OUT OF DO-LOOP
        GOTO 100
      END IF
    END DO
    !
    !           SHIFT THE TABLE.
    !
    IF( N==limexp ) N = 2*(limexp/2) - 1
    ib = 1
    IF( (num/2)*2==num ) ib = 2
    ie = newelm + 1
    DO i = 1, ie
      ib2 = ib + 2
      Epstab(ib) = Epstab(ib2)
      ib = ib2
    END DO
    IF( num/=N ) THEN
      indx = num - N + 1
      DO i = 1, N
        Epstab(i) = Epstab(indx)
        indx = indx + 1
      END DO
    END IF
    IF( Nres>=4 ) THEN
      !
      !           COMPUTE ERROR ESTIMATE
      !
      Abserr = ABS(Result-Res3la(3)) + ABS(Result-Res3la(2)) + ABS(Result-Res3la(1))
      Res3la(1) = Res3la(2)
      Res3la(2) = Res3la(3)
      Res3la(3) = Result
    ELSE
      Res3la(Nres) = Result
      Abserr = oflow
    END IF
  END IF
  100  Abserr = MAX(Abserr,5._DP*epmach*ABS(Result))
END SUBROUTINE DQELG
