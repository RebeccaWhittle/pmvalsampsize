from pmvalsampsize.pmvalsampsize import *


def test_pmvalsampsize():

    #Test 1
    expected = [905, 4544, 4252, 4544]
    test = pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, 
                         lpnormal = (-5,2.5), oeciwidth = 1, noprint=True)
    actual = [test["results_table"][0][1], test["results_table"][1][1], 
          test["results_table"][2][1], test["results_table"][4][1]]
    assert actual == expected
    

    #Test 2
    expected = [905, 1644, 4252, 4252]
    test = pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, 
                         lpbeta = (0.5,0.5), oeciwidth = 1, noprint=True)
    actual = [test["results_table"][0][1], test["results_table"][1][1], 
          test["results_table"][2][1], test["results_table"][4][1]]
    assert actual == expected

    #Test 3
    expected = [905, 17080, 4252, 17080]
    test = pmvalsampsize(type = "b", prevalence = 0.018, cstatistic = 0.8, 
                         lpcstat = -4.7, oeciwidth = 1, noprint=True)
    actual = [test["results_table"][0][1], test["results_table"][1][1], 
          test["results_table"][2][1], test["results_table"][4][1]]
    assert actual == expected
    