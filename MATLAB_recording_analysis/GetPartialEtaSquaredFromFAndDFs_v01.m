function pE2 = GetPartialEtaSquaredFromFAndDFs_v01(Fstat,DFeffect,DFerror)


pE2 = (Fstat*DFeffect) / ( (Fstat*DFeffect) + DFerror );

end % of function