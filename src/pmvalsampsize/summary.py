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