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
                      specificity,threshold,nbciwidth,nbseincrement): 



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
        print("\n", "Normal LP distribution with parameters - mean =",
              lpnormal[0],", sd =",lpnormal[1],"\n")
        LP = norm.rvs(loc=lpnormal[0], scale=lpnormal[1], size=simobs)

    elif not pd.isna(lpbeta):
        lpdist = "beta"
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
            
                        print("Proportion of outcome events under simulation =",
                              outcome_test.mean(),"\n","Target prevalence =",
                              prevalence,"\n","Mean in non-event group = ",
                              m2,"\n")
                print("-------------- TRACE OFF --------------","\n")
                print("Proportion of observed outcome events is within tolerance",
                      "\n","Proportion of outcome events under simulation =",
                      outcome_test.mean(),"\n","Target prevalence =",
                      prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
        
            else:
                print("\n", "Proportion of observed outcome events is within tolerance",
                      "\n","Proportion of outcome events under simulation =",
                      outcome_test.mean(),"\n","Target prevalence =",
                      prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
    
        else:
            if (diff > tolerance):
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
    
    
                print("\n","Proportion of observed outcome events is within tolerance",
                      "\n","Proportion of outcome events under simulation =",
                      outcome_test.mean(),"\n","Target prevalence =",
                      prevalence,"\n","Mean in non-event group = ",m2,"\n","\n")
    
            else:
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
    Borenstein_11 = LP * LP * np.exp(beta0 + (beta1*LP)) / ((1 + np.exp(beta0 + (beta1*LP)))**2)

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
    cstat_df['se_cstatsq'] = cstatistic*(1-cstatistic)*(1+(((cstat_df['size']/2)-1)*((1-cstatistic)/(2-cstatistic))) +((((cstat_df['size']/2)-1)*cstatistic)/(1+cstatistic)))/(cstat_df['size']*cstat_df['size']*prevalence*(1-prevalence))
    cstat_df['se_cstat'] = cstat_df['se_cstatsq']**0.5
    cstat_df['CIwidth'] = 2*1.96*cstat_df['se_cstat']

    cstat_df2 = cstat_df.loc[cstat_df.CIwidth<=cstatciwidth]
    n3 = cstat_df2['size'].min()

    se_cstat = math.sqrt(cstatistic*(1-cstatistic)*(1+(((n3/2)-1)*((1-cstatistic)/(2-cstatistic))) +((((n3/2)-1)*cstatistic)/(1+cstatistic)))/(n3*n3*prevalence*(1-prevalence)))


#############################
# criteria 4 - net benefit


    if sensitivity is not None and specificity is not None:
        nb = (sensitivity*prevalence) - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
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

        n4 = math.ceil((1/(se_nb**2))*(((sensitivity*(1-sensitivity))/prevalence)+(w*w*specificity*(1-specificity)/(1-prevalence))+(w*w*(1-specificity)*(1-specificity)/(prevalence*(1-prevalence)))))

        no_nb = False

    elif not pd.isna(threshold): 
        nb_p = np.exp(LP)/(1+np.exp(LP))
        nb_outcome = np.random.binomial(n=1, p=nb_p, size=simobs)

        nb_df = pd.DataFrame(data = {'nb_outcome': nb_outcome, 'nb_p': nb_p})
        nb_df['classification'] = 1
        nb_df.loc[nb_df['nb_p'] < threshold, 'classification'] = 0

        sensitivity = round(nb_df.classification[nb_df.nb_outcome==1].sum() / (nb_df['nb_outcome']==1).sum(), 3)
        specificity = round(((nb_df['nb_outcome']==0).sum() - nb_df.classification[nb_df.nb_outcome==0].sum()) / (nb_df['nb_outcome']==0).sum(), 3)
 
        nb = (sensitivity*prevalence) - ((1-specificity)*(1-prevalence)*(threshold/(1-threshold)))
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
        n4 = math.ceil((1/(se_nb**2))*(((sensitivity*(1-sensitivity))/prevalence)+(w*w*specificity*(1-specificity)/(1-prevalence))+(w*w*(1-specificity)*(1-specificity)/(prevalence*(1-prevalence)))))

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



### error check function 
 
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
 



### wrapper function 
 
def pmvalsampsize(type, prevalence=None, cstatistic=None, oe=1, oeciwidth=0.2, 
                  cslope=1, csciwidth=0.2, cstatciwidth=0.1, simobs=1000000,
                  lpnormal=None, lpbeta=None, lpcstat=None, tolerance=5e-04, 
                  increment=0.1, oeseincrement=1e-04, seed=123456, 
                  graph=False,trace=False,sensitivity=None, specificity=None, 
                  threshold=None, nbciwidth=0.2, nbseincrement=1e-04): 

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
                             nbciwidth=nbciwidth,nbseincrement=nbseincrement)

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
                                nbseincrement=nbseincrement)
    return out




### Summary ouput

def summary(x, *args):


    col_names = ["Criteria", "Sample size", "Perf", "SE", "CI width"]
    
    print("\n", tabulate(x["results_table"], headers=col_names, numalign="right"))
    
    if x["type"] == "binary":
        print("\n", "Minimum sample size required for model validation based ",
              "on user inputs = ", x["sample_size"], ",", "\n", "with ", 
              math.ceil(x["events"])," events (assuming an outcome prevalence = "
              , x["prevalence"],")","\n", sep='')