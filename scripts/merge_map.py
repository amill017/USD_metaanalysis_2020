import pandas as pd
import numpy as np
import sys
df1 = pd.read_csv('data/current_map.txt', sep='\t')
df2 = pd.read_csv(sys.argv[1], sep='\t')
result = df1.append(df2, sort=False)
result2 = result.replace('', np.NaN)
result2.to_csv('data/current_map.txt', sep='\t', na_rep='blank')