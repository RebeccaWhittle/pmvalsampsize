def pmvalsampsize(type, prevalence=None, cstatistic=None, oe=1, oeciwidth=0.2, 
                  cslope=1, csciwidth=0.2, cstatciwidth=0.1, simobs=1000000,
                  lpnormal=None, lpbeta=None, lpcstat=None, tolerance=5e-04, 
                  increment=0.1, oeseincrement=1e-04, seed=123456, 
                  graph=False,trace=False,sensitivity=None, specificity=None, 
                  threshold=None, nbciwidth=0.2, nbseincrement=1e-04, noprint=None): 
    """Computes the minimum sample size required for the external validation of an existing multivariable prediction model
 
    Parameters
    ----------
    type: str
        specifies the type of analysis for which sample size is being calculated
        "b" specifies sample size calculation for a prediction model with a binary outcome
    cslope: float, default=1
        specifies the anticipated c-slope performance in the validation sample. 
        The value could alternatively be based on a previous validation study for example. 
        For binary outcomes the c-slope calculation requires the user to specify a distribution 
        for the assumed LP in the validation sample 
        (or alternatively the distribution of predicted probabilities in the validation sample). 
        See lp*() options below.
    csciwidth: float, default=0.2
        specifies the target CI width (acceptable precision) for the c-slope performance. 
    oe: float, default=1
        specifies the anticipated O/E performance in the validation sample. 
    oeciwidth: float, default=0.2
        specifies the target CI width (acceptable precision) for the E/O performance. 
        The choice of CI width is context specific, 
        and depends on the event probability in the population. 
        See Riley et al. for further details.
    cstatistic: float
        specifies the anticipated c-statistic performance in the validation sample. 
        This is a required input. 
        May be based on the optimism-adjusted c-statistic reported in the development study 
        for the existing prediction model. Ideally, this should be an optimism-adjusted c-statistic. 
        NB: This input is also used when using the lpcstat() option.
    cstatciwidth: float, default=0.1
        specifies the target CI width (acceptable precision) for the c-statistic performance. 
    simobs: int, default=1000000
        specifies the number of observations to use when simulating the 
        LP distribution for c-slope calculation in criteria 2. 
        Higher simobs() values will reduce random variation further. 
    lpnormal: float, optional
        defines parameters to simulate the LP distribution for criteria 2 from a normal distribution. 
        The user must specify the mean and standard deviation (in this order) of the LP distribution.
    lpbeta: float, optional
        defines parameters to simulate the distribution of predicted probabilities
        for criteria 2 from a beta distribution. 
        The user must specify the alpha and beta parameters (in this order) of the
        probability distribution. 
        The LP distribution is then generated internally using this probability distribution.
    lpcstat: float, optional
        defines parameters to simulate the LP distribution for criteria 2 assuming
        that the distribution of events and non-events are normal with a common variance. 
        The user specifies a single input value- the expected mean for the non-events distribution.
        This could be informed by clinical guidance. 
        However, this input is taken as a starting value and an iterative process 
        is used to identify the most appropriate values for the event and non-event
        distributions so as to closely match the anticipated prevalence in the 
        validation sample. 
        NB: this approach makes strong assumptions of normality and equal variances
        in each outcome group, which may be unrealistic in most situations.
    tolerance: float, default=0.0005
        for use with lpcstat option. Sets the tolerance for agreement between the 
        simulated and expected event proportion during the iterative procedure 
        for calculating the mean for the non-events distribution.
    increment: float, default=0.1
        for use with lpcstat option. Sets increment by which to iterate the value 
        of the mean for the non-events distribution. 
        Trial and error may be necessary as it is dependent on how close the 
        initial input for the non-event mean in lpcstat is to the required value. 
        If the process takes a particularly long time then the user could try an 
        alternative increment value, or an alternative non-event mean value in lpcstat. 
        The trace option may be useful in such circumstances.
    oeseincrement: float, default=0.0001
        sets the increment by which to iterate when identifying the SE(ln(OE)) to 
        meet the target CI width specified for OE. 
        In the majority of cases this will be suitably small to ensure a precise 
        SE is identified. 
        The user should check the output table to ensure that the target CI width 
        has been attained and adjust the increment if necessary.
    graph: bool, default=False
        specifies that a histogram of the simulated LP distribution for criteria 2 is produced. 
        The graph also details summary statistics for the simulated distribution. 
        Useful option for checking the simulated LP distribution against the source of input parameters. Also useful for reporting at publication.
    trace: bool, default=False
        for use with `lpcstat` option. 
        Specifies that a trace of the values obtained in each iteration when 
        identifying the non-event mean is output. 
        Useful when finding the appropriate values for `lpcstat` & `increment` is proving difficult!
    prevalence: float, optional
        specifies the overall outcome proportion (for a prognostic model) or 
        overall prevalence (for a diagnostic model) expected 
        within the model validation sample. 
        This is a required input. 
        This should be derived based on previous studies in the same population 
        or directly from the validation sample if to hand.
    seed: int, default=123456
        specifies the initial value of the random-number seed used by the 
        random-number functions when simulating data to approximate the LP 
        distribution for criteria 2.
    sensitivity: float, optional
        specifies the anticipated sensitivity performance in the validation sample 
        at the chosen risk threshold (specified using `threshold`). 
        If sensitivity and specificity are not provided then `pmvalsampsize` 
        uses the simulated LP distribution from criteria 2 and the user-specified 
        risk threshold to estimate the anticipated sensitivity and specificity 
        to be used in calculation of net benefit. 
        NB: net benefit criteria is not calculated if either
        i) `sensitivity`, `specificity` and `threshold` 
        or ii) `threshold` option are not specified.
    specificity: float, optional
        specifies the anticipated specificity performance in the validation sample 
        at the chosen risk threshold (specified using `threshold`). 
        If sensitivity and specificity are not provided then `pmvalsampsize` 
        uses the simulated LP distribution from criteria 2 and the user-specified 
        risk threshold to estimate the anticipated sensitivity and specificity 
        to be used in calculation of net benefit. 
        NB: net benefit criteria is not calculated if either
        i) `sensitivity`, `specificity` and `threshold` 
        or ii) `threshold` option are not specified.
    threshold: float, optional
        specifies the risk threshold to be used for calculation of net benefit 
        performance of the model in the validation sample. 
        If `sensitivity` and `specificity` are not provided then `threshold` must 
        be given in order for `pmvalsampsize` to assess sample size 
        requirements for net benefit. 
        NB: net benefit criteria is not calculated if either
        i) `sensitivity`, `specificity` and `threshold` 
        or ii) `threshold` option are not specified.
    nbciwidth: float, default=0.2
        specifies the target CI width (acceptable precision) for the standardised 
        net benefit performance. 
        The choice of CI width is context specific. See Riley et al. for further details.   
    nbseincrement: float, default=0.0001
        sets the increment by which to iterate when identifying the 
        SE(standardised net benefit) to meet the target CI width specified 
        for standardised net benefit. 
        In the majority of cases this will be suitably small to ensure a 
        precise SE is identified. The user should check the output table to ensure 
        that the target CI width has been attained and adjust the increment if necessary.
    noprint: bool, default=False
        supresses output being printed    
    """
    
#error checking
    pmvalsampsize_errorcheck(type=type,prevalence=prevalence,
                             cstatistic=cstatistic,oe=oe,oeciwidth=oeciwidth,
                             cslope=cslope,csciwidth=csciwidth,
                             cstatciwidth=cstatciwidth,simobs=simobs,
                             lpnormal=lpnormal,lpbeta=lpbeta,lpcstat=lpcstat,
                             tolerance=tolerance,seed=seed,graph=graph,
                             oeseincrement=oeseincrement,increment=increment,
                             trace=trace,sensitivity=sensitivity,
                             specificity=specificity,threshold=threshold,
                             nbciwidth=nbciwidth,nbseincrement=nbseincrement, noprint=noprint)

    if type == "b":
        out = pmvalsampsize_bin(prevalence=prevalence,cstatistic=cstatistic,
                                oe=oe,oeciwidth=oeciwidth,cslope=cslope,
                                csciwidth=csciwidth,cstatciwidth=cstatciwidth,
                                simobs=simobs,lpnormal=lpnormal,lpbeta=lpbeta,
                                lpcstat=lpcstat,tolerance=tolerance,seed=seed,
                                increment=increment,graph=graph,trace=trace,
                                oeseincrement=oeseincrement,
                                sensitivity=sensitivity,nbciwidth=nbciwidth,
                                specificity=specificity,threshold=threshold,
                                nbseincrement=nbseincrement, noprint=noprint)
    return out


# error check function
import pandas as pd
def pmvalsampsize_errorcheck(type,prevalence,cstatistic,oe,oeciwidth,cslope,
                             csciwidth,cstatciwidth,simobs,lpnormal,lpbeta,
                             lpcstat,tolerance,increment,oeseincrement,seed,
                             graph,trace,sensitivity,specificity,threshold,
                             nbciwidth,nbseincrement,noprint): 

   
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
        


# Binary option
import itertools as it
import math 
import random 
import matplotlib
import matplotlib.pyplot as plt
from tabulate import tabulate, SEPARATING_LINE

import numpy as np
import pandas as pd
from scipy.stats import norm, beta
from scipy.special import logit
from sklearn.metrics import (accuracy_score, confusion_matrix, 
                                 roc_auc_score, roc_curve)


### binary option

def pmvalsampsize_bin(prevalence,cstatistic,oe,oeciwidth,cslope,csciwidth,
                      cstatciwidth,simobs,lpnormal,lpbeta,lpcstat,tolerance,
                      increment,oeseincrement,seed,graph,trace,sensitivity,
                      specificity,threshold,nbciwidth,nbseincrement,noprint): 



    np.random.seed(seed)
 
#############################
  # criteria 1 - o/e

    width_oe = 0
    se_oe = 0
    while (width_oe < oeciwidth):
        se_oe += oeseincrement
        ub_oe = np.exp(np.log(oe) + (1.96*se_oe))
        lb_oe = np.exp(np.log(oe) - (1.96*se_oe))
        width_oe = ub_oe - lb_oe

    n1 = math.ceil((1-prevalence)/(prevalence*(se_oe**2)))
    E1 = n1*prevalence
 
#############################
# criteria 2 - c-slope

    if not pd.isna(lpnormal):
        lpdist = "normal"
        if noprint==None:
            print("\n", "Normal LP distribution with parameters - mean =",
                  lpnormal[0],", sd =",lpnormal[1],"\n")
        LP = norm.rvs(loc=lpnormal[0], scale=lpnormal[1], size=simobs)

    elif not pd.isna(lpbeta):
        lpdist = "beta"
        if noprint==None:
            print("\n", "Beta P distribution with parameters - alpha =", lpbeta[0],
                  ", beta =", lpbeta[1], "\n")
        P = beta.rvs(a=lpbeta[0], b=lpbeta[1], size=simobs)
        LP = np.log(P/(1-P))
   

    elif not pd.isna(lpcstat):
        lpdist = "cstat"
        m2 = lpcstat
        var = 2*(norm.ppf(cstatistic)**2)

        outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
    
        LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
        LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
    
        LP[outcome==1] = LP2[outcome==1]

        P = np.exp(LP)/(1+np.exp(LP))
        outcome_test = np.random.binomial(n=1, p=P, size=simobs)
        diff = abs(prevalence - outcome_test.mean())
   
        if trace==True:
          
            if (diff>tolerance):
                if noprint==None:
                    print("\n",
                          "Proportion of observed outcome events does not match",
                          " input prevalence","\n",
                          "Beginning iterative approach ...")
                    print("-------------- TRACE ON --------------","\n")
        
                n = 1
                diff_new = diff
                m2 += increment
            
                outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
            
                LP[outcome==1] = LP2[outcome==1]
              
                P = np.exp(LP)/(1+np.exp(LP))
                outcome_test = np.random.binomial(n=1, p=P, size=simobs)
    
                diff_new = abs(prevalence - outcome_test.mean())
        
                if (diff_new < diff):
                    while (diff_new > tolerance):
                        m2 += increment
            
                        outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                        LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                        LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
            
                        LP[outcome==1] = LP2[outcome==1]
            
                        P = np.exp(LP)/(1+np.exp(LP))
                        outcome_test = np.random.binomial(n=1, p=P, size=simobs)
            
                        diff_new = abs(prevalence - outcome_test.mean())
            
                        if noprint==None:
                            print("Proportion of outcome events under simulation =",
                                  outcome_test.mean(),"\n","Target prevalence =",
                                  prevalence,"\n","Mean in non-event group = ",m2)
            
                else:
                    while (diff_new > tolerance):
                        m2 -= increment
            
                        outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                        LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                        LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
            
                        LP[outcome==1] = LP2[outcome==1]
              
                        P = np.exp(LP)/(1+np.exp(LP))
                        outcome_test = np.random.binomial(n=1, p=P, size=simobs)
            
                        diff_new = abs(prevalence - outcome_test.mean())
            
                        if noprint==None:
                            print("Proportion of outcome events under simulation =",
                                  outcome_test.mean(),"\n","Target prevalence =",
                                  prevalence,"\n","Mean in non-event group = ",
                                  m2,"\n")
                if noprint==None:
                    print("-------------- TRACE OFF --------------","\n")
                    print("Proportion of observed outcome events is within tolerance",
                      "\n","Proportion of outcome events under simulation =",
                      outcome_test.mean(),"\n","Target prevalence =",
                      prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
        
            else:
                if noprint==None:
                    print("\n", "Proportion of observed outcome events is within tolerance",
                          "\n","Proportion of outcome events under simulation =",
                          outcome_test.mean(),"\n","Target prevalence =",
                          prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
    
        else:
            if (diff > tolerance):
                if noprint==None:
                    print("\n","Proportion of observed outcome events does not ",
                          "match input prevalence","\n",
                          "Beginning iterative approach ...","\n")
    
                n = 1
                diff_new = diff
                m2 += increment
    
                outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
    
                LP[outcome==1] = LP2[outcome==1]
    
                P = np.exp(LP)/(1+np.exp(LP))
                outcome_test = np.random.binomial(n=1, p=P, size=simobs)
        
                diff_new = abs(prevalence - outcome_test.mean())
    
                if (diff_new < diff):
                    while (diff_new > tolerance):
                        m2 += increment
    
                        outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                        LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                        LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
    
                        LP[outcome==1] = LP2[outcome==1]
    
                        P = np.exp(LP)/(1+np.exp(LP))
                        outcome_test = np.random.binomial(n=1, p=P, size=simobs)
    
                        diff_new = abs(prevalence - outcome_test.mean())
    
    
                else:
                    while (diff_new > tolerance):
                        m2 -= increment
    
                        outcome = np.random.binomial(n=1, p=prevalence, size=simobs)
                        LP = norm.rvs(loc=m2, scale=math.sqrt(var), size=simobs)
                        LP2 = norm.rvs(loc=m2+var, scale=math.sqrt(var), size=simobs)
    
                        LP[outcome==1] = LP2[outcome==1]
    
                        P = np.exp(LP)/(1+np.exp(LP))
                        outcome_test = np.random.binomial(n=1, p=P, size=simobs)
    
                        diff_new = abs(prevalence - outcome_test.mean())
    
    
                if noprint==None:
                    print("\n","Proportion of observed outcome events is within tolerance",
                          "\n","Proportion of outcome events under simulation =",
                          outcome_test.mean(),"\n","Target prevalence =",
                          prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
    
            else:
                if noprint==None:
                    print("\n","Proportion of observed outcome events is within tolerance",
                          "\n","Proportion of outcome events under simulation =",
                          outcome_test.mean(),"\n","Target prevalence =",
                          prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
   
# check c-statistic
        simulated_data_cstat = roc_auc_score(outcome_test, P)
   
### histogram
# graph = True
    if (graph==True):
        plt.hist(LP, density=1, bins=50)
	 
# input assumed parameters of calibration model (in future vr these could be options)
    beta0 = 0
    beta1 = 1

# calculate elements of I matrix
    Borenstein_00 = np.exp(beta0 + (beta1*LP)) / ((1 + np.exp(beta0 + (beta1*LP)))**2)
    Borenstein_01 = LP * np.exp(beta0 + (beta1*LP)) / ((1+ np.exp(beta0 + (beta1*LP)))**2)
    Borenstein_11 = (LP * LP * np.exp(beta0 + (beta1*LP)) 
                     / ((1 + np.exp(beta0 + (beta1*LP)))**2))

    I_00 = Borenstein_00.mean()  
    I_01 = Borenstein_01.mean()
    I_11 = Borenstein_11.mean()

# calculate SE from input target CI width
    se_cslope = csciwidth/(2*1.96)

# calculate sample size
    n2 = math.ceil(I_00 / (se_cslope*se_cslope * ((I_00*I_11) - (I_01*I_01))))
    E2 = n2*prevalence
   
   
#############################
# criteria 3 - c-statistic


    cstat_df = pd.DataFrame(index=np.arange(1000000))
    cstat_df['size'] = cstat_df.index+1
    cstat_df['se_cstatsq'] = (cstatistic
                              *(1-cstatistic)*(1+(((cstat_df['size']/2)-1)
                                                  *((1-cstatistic)/(2-cstatistic))) 
                                                  +((((cstat_df['size']/2)-1)
                                                     *cstatistic)/(1+cstatistic)))
                                                     /(cstat_df['size']
                                                       *cstat_df['size']
                                                       *prevalence*(1-prevalence)))
    cstat_df['se_cstat'] = cstat_df['se_cstatsq']**0.5
    cstat_df['CIwidth'] = 2*1.96*cstat_df['se_cstat']

    cstat_df2 = cstat_df.loc[cstat_df.CIwidth<=cstatciwidth]
    n3 = cstat_df2['size'].min()

    se_cstat = math.sqrt(cstatistic*
                         (1-cstatistic)*
                         (1+(((n3/2)-1)*((1-cstatistic)
                                         /(2-cstatistic))) 
                          +((((n3/2)-1)*cstatistic)/(1+cstatistic)))
                          /(n3*n3*prevalence*(1-prevalence)))


#############################
# criteria 4 - net benefit


    if sensitivity is not None and specificity is not None:
        nb = (sensitivity*prevalence) 
        - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
        standardised_nb = nb/prevalence

        w = ((1-prevalence)/prevalence)*(threshold/(1-threshold))

        width_nb = 0
        se_nb = 0
        while width_nb<nbciwidth:
            se_nb += nbseincrement
            ub_nb = (standardised_nb) + (1.96*se_nb)
            lb_nb = (standardised_nb) - (1.96*se_nb)
            width_nb = ub_nb-lb_nb
        
        se_nb = round(se_nb,3)

        n4 = math.ceil((1/(se_nb**2))
                       *(((sensitivity*(1-sensitivity))/prevalence)
                         +(w*w*specificity*(1-specificity)/(1-prevalence))
                         +(w*w*(1-specificity)*(1-specificity)
                           /(prevalence*(1-prevalence)))))

        no_nb = False

    elif not pd.isna(threshold): 
        nb_p = np.exp(LP)/(1+np.exp(LP))
        nb_outcome = np.random.binomial(n=1, p=nb_p, size=simobs)

        nb_df = pd.DataFrame(data = {'nb_outcome': nb_outcome, 'nb_p': nb_p})
        nb_df['classification'] = 1
        nb_df.loc[nb_df['nb_p'] < threshold, 'classification'] = 0

        sensitivity = round(nb_df.classification[nb_df.nb_outcome==1].sum() 
                            / (nb_df['nb_outcome']==1).sum(), 3)
        specificity = round(((nb_df['nb_outcome']==0).sum() 
                             - nb_df.classification[nb_df.nb_outcome==0].sum()) 
                             / (nb_df['nb_outcome']==0).sum(), 3)
 
        nb = (sensitivity*prevalence) 
        - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
        standardised_nb = nb/prevalence

        w = ((1-prevalence)/prevalence)*(threshold/(1-threshold))

   # calc se for target ci width
        width_nb = 0
        se_nb = 0
        while width_nb<nbciwidth:
            se_nb += nbseincrement
            ub_nb = (standardised_nb) + (1.96*se_nb)
            lb_nb = (standardised_nb) - (1.96*se_nb)
            width_nb = ub_nb-lb_nb

        se_nb = round(se_nb, 3)

   # calculate sample size
        n4 = math.ceil((1/(se_nb**2))
                       *(((sensitivity*(1-sensitivity))/prevalence)
                         +(w*w*specificity*(1-specificity)/(1-prevalence))
                         +(w*w*(1-specificity)*(1-specificity)
                           /(prevalence*(1-prevalence)))))

        no_nb = False

    else:  
        no_nb = True
   
   
   
####### summary
    if no_nb==True:
 # minimum n
        nfinal = max(n1,n2,n3)
        efinal = nfinal*prevalence       
   
# create output table
        res = [["Criteria 1 - O/E", n1, round(oe, 3), round(se_oe, 3), round(width_oe, 3)], 
              ["Criteria 2 - C-slope", n2, round(cslope, 3), round(se_cslope, 3), round(csciwidth, 3)], 
              ["Criteria 3 - C statistic", n3, round(cstatistic, 3), round(se_cstat, 3), round(cstatciwidth, 3)], 
              ["Final SS", 1, 1, 1, 1]]

        col_names = ["Criteria", "Sample size", "Perf", "SE", "CI width"]
 

        res_sort = res.sort(key = lambda res:res[1], reverse=False)
 
        res = [["Criteria 1 - O/E", n1, round(oe, 3), round(se_oe, 3), round(width_oe, 3)], 
              ["Criteria 2 - C-slope", n2, round(cslope, 3), round(se_cslope, 3), round(csciwidth, 3)], 
              ["Criteria 3 - C statistic", n3, round(cstatistic, 3), round(se_cstat, 3), round(cstatciwidth, 3)], 
              SEPARATING_LINE,
              ["Final SS", nfinal, res[3][2], res[3][3], res[3][4]]]
        
        if noprint==None:
            print(tabulate(res, headers=col_names, numalign="right"))
            print("\n","Minimum sample size required for model validation based on user inputs = ",nfinal,",","\n","with ",math.ceil(efinal)," events (assuming an outcome prevalence = ",prevalence,")", sep='')
            print("\n","Criteria 1 - precise estimation of O/E performance in the validation sample","\n",
                  "Criteria 2 - precise estimation of the calibration slope in the validation sample","\n",
                  "Criteria 3 - precise estimation of the C statistic in the validation sample","\n")


        if not pd.isna(lpcstat):
            out = {
                   "results_table": res,
                   "sample_size": nfinal,
                   "events": efinal,
                   "prevalence": prevalence,
                   "type": "binary",
                   "cstatistic": cstatistic,
                   "oe": oe,
                   "se_oe": se_oe,
                   "width_oe": width_oe,
                   "cslope": cslope,
                   "se_cslope": se_cslope,
                   "csciwidth": csciwidth,
                   "se_cstat": se_cstat,
                   "cstatciwidth": cstatciwidth,
                   "simulated_data_cstat": simulated_data_cstat,
                   }
        else:
            out = {
                   "results_table": res,
                   "sample_size": nfinal,
                   "events": efinal,
                   "prevalence": prevalence,
                   "type": "binary",
                   "cstatistic": cstatistic,
                   "oe": oe,
                   "se_oe": se_oe,
                   "width_oe": width_oe,
                   "cslope": cslope,
                   "se_cslope": se_cslope,
                   "csciwidth": csciwidth,
                   "se_cstat": se_cstat,
                   "cstatciwidth": cstatciwidth,
                   }        
    else: 
# minimum n
        nfinal = max(n1,n2,n3,n4)
        efinal = nfinal*prevalence  

# create output table
        res = [["Criteria 1 - O/E", n1, round(oe, 3), round(se_oe, 3), round(width_oe, 3)], 
              ["Criteria 2 - C-slope", n2, round(cslope, 3), round(se_cslope, 3), round(csciwidth, 3)], 
              ["Criteria 3 - C statistic", n3, round(cstatistic, 3), round(se_cstat, 3), round(cstatciwidth, 3)],
              ["Criteria 4 - St Net Benefit", n4, round(standardised_nb, 3), round(se_nb, 3), round(nbciwidth, 3)],
              ["Final SS", 1, 1, 1, 1]]

        col_names = ["Criteria", "Sample size", "Perf", "SE", "CI width"]
 

        res_sort = res.sort(key = lambda res:res[1], reverse=False)
 
        res = [["Criteria 1 - O/E", n1, round(oe, 3), round(se_oe, 3), round(width_oe, 3)], 
              ["Criteria 2 - C-slope", n2, round(cslope, 3), round(se_cslope, 3), round(csciwidth, 3)], 
              ["Criteria 3 - C statistic", n3, round(cstatistic, 3), round(se_cstat, 3), round(cstatciwidth, 3)], 
              ["Criteria 4 - St Net Benefit", n4, round(standardised_nb, 3), round(se_nb, 3), round(nbciwidth, 3)],
              SEPARATING_LINE,
              ["Final SS", nfinal, res[4][2], res[4][3], res[4][4]]]
        
        if noprint==None:

            print(tabulate(res, headers=col_names, numalign="right"))
            print("\n","Minimum sample size required for model validation based on user inputs = ",nfinal,",","\n","with ",math.ceil(efinal)," events (assuming an outcome prevalence = ",prevalence,")", sep='')
            print("\n","Criteria 1 - precise estimation of O/E performance in the validation sample","\n",
                  "Criteria 2 - precise estimation of the calibration slope in the validation sample","\n",
                  "Criteria 3 - precise estimation of the C statistic in the validation sample","\n",
                  "Criteria 4 - precise estimation of the standardised net-benefit in the validation sample","\n")
  

        if not pd.isna(lpcstat):
            out = {
                   "results_table": res,
                   "sample_size": nfinal,
                   "events": efinal,
                   "prevalence": prevalence,
                   "type": "binary",
                   "cstatistic": cstatistic,
                   "oe": oe,
                   "se_oe": se_oe,
                   "width_oe": width_oe,
                   "cslope": cslope,
                   "se_cslope": se_cslope,
                   "csciwidth": csciwidth,
                   "se_cstat": se_cstat,
                   "cstatciwidth": cstatciwidth,
                   "simulated_data_cstat": simulated_data_cstat,
                   "standardised_nb": standardised_nb,
                   "net_benefit": nb,
                   "se_st_nb": se_nb,
                   "nbciwidth": nbciwidth,
                   "sensitivity": sensitivity,
                   "specificity": specificity,
                   "threshold": threshold
                   }
        else:
            out = {
                   "results_table": res,
                   "sample_size": nfinal,
                   "events": efinal,
                   "prevalence": prevalence,
                   "type": "binary",
                   "cstatistic": cstatistic,
                   "oe": oe,
                   "se_oe": se_oe,
                   "width_oe": width_oe,
                   "cslope": cslope,
                   "se_cslope": se_cslope,
                   "csciwidth": csciwidth,
                   "se_cstat": se_cstat,
                   "cstatciwidth": cstatciwidth,
                   "standardised_nb": standardised_nb,
                   "net_benefit": nb,
                   "se_st_nb": se_nb,
                   "nbciwidth": nbciwidth,
                   "sensitivity": sensitivity,
                   "specificity": specificity,
                   "threshold": threshold
                   }
    return out


### Summary ouput
import math 
from tabulate import tabulate, SEPARATING_LINE


def summary(x, *args):


    col_names = ["Criteria", "Sample size", "Perf", "SE", "CI width"]
    
    print("\n", tabulate(x["results_table"], headers=col_names, numalign="right"))
    
    if x["type"] == "binary":
        print("\n", "Minimum sample size required for model validation based ",
              "on user inputs = ", x["sample_size"], ",", "\n", "with ", 
              math.ceil(x["events"])," events (assuming an outcome prevalence = "
              , x["prevalence"],")","\n", sep='')