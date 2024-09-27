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
