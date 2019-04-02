MODULE DDEBD1
  IMPLICIT NONE
  REAL(8) :: TOLd, CONit, CRAte, EL(13), ELCo(13,12), HOLd, RC, RMAx, TESco(3,12), &
    EL0, H, HMIn, HMXi, HU, TN, UROund
  INTEGER :: IQUit, INIt, IYH, IEWt, IACor, ISAvf, IWM, KSTeps, IALth, IPUp, &
    LMAx, MEO, NQNyh, NSTepj, IBEgin, ITOl, IINteg, ITStop, IJAc, IBAnd, IER, &
    JSTart, KFLag, L, METh, MITer, MAXord, N, NQ, NST, NFE, NJE, NQU
END MODULE DDEBD1