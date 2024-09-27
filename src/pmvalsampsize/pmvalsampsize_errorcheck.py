### error check function 
import pandas as pd
 
def pmvalsampsize_errorcheck(type,prevalence,cstatistic,oe,oeciwidth,cslope,
                             csciwidth,cstatciwidth,simobs,lpnormal,lpbeta,
                             lpcstat,tolerance,increment,oeseincrement,seed,
                             graph,trace,sensitivity,specificity,threshold,
                             nbciwidth,nbseincrement): 

   
    if type not in ["c", "b", "s"]:
        raise ValueError('type must be "c" for continuous, "b" for binary, or "s" for survival')
    if not isinstance(simobs, int):
        raise ValueError('simobs must be an integer')
    if simobs != round(simobs):
        raise ValueError('simobs must be an integer')
 
# parameters for binary
    if type == "b":
  # parameters needed
        if pd.isna(prevalence):
            raise ValueError('prevalence must be specified for binary sample size')
        if pd.isna(cstatistic):
            raise ValueError('cstatistic must be specified for binary outcome models')

        if not pd.isna(lpnormal):
            if not pd.isna(lpbeta) or not pd.isna(lpcstat):
                raise ValueError('Only one LP distribution option can be specified')

        elif not pd.isna(lpbeta):
            if not pd.isna(lpcstat):
                raise ValueError('Only one LP distribution option can be specified')
    
        elif pd.isna(lpnormal) and pd.isna(lpbeta) and pd.isna(lpcstat):
            raise ValueError('An LP distribution must be specified')

# parameter conditions
        if not isinstance(prevalence, (int, float)):
            raise ValueError('prevalence must be numeric')
        if cstatistic < 0 or cstatistic > 1:
            raise ValueError('cstatistic must be between 0 and 1')
        if not isinstance(cstatistic, (int, float)):
            raise ValueError('cstatistic must be numeric')
        if not isinstance(cslope, (int, float)):
            raise ValueError('cslope must be numeric')
        if prevalence <= 0 or prevalence >= 1:
            raise ValueError('prevalence must be between 0 and 1')