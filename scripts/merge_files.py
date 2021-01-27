import pandas as pd
import numpy as np
import sys
df1 = pd.read_csv(sys.argv[1], sep='\t')
df2 = pd.read_csv(sys.argv[2], sep='\t')
result = pd.concat([df1,df2], sort=False, axis=1)
result2 = result.replace('', np.NaN)
result2.to_csv(sys.argv[3], sep='\t', na_rep='NA')