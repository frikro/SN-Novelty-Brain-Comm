import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import mixedlm
import pingouin as pg
df=pd.read_csv('out.csv')
df["diag"] = df["diag"].astype("category")
df["sub"] = df["sub"].astype("category")
model = mixedlm("SN_Volume_raw_man_ts_mask ~ diag * time", data=df, groups=df["subject"]).fit()
print(model.summary())

# Pairwise comparisons
posthoc = pg.pairwise_tests(dv="SN_Volume_raw_man_ts_mask", between="diag", subject='sub',data=df, parametric=True,padjust='sidak')
print(posthoc)


