disease = []
dv = []
tmps = []
with open('./server_onehot_DS_ALL.txt') as f:
    for i in f:
        tmp = str(i).split('\t')
        t = eval(tmp[2])
        if t in tmps:
            # print(t)
            pass
        else:
            tmps.append(t)
            for i in t:
                if i in disease:
                    pass
                else:
                    disease.append(i)
                    dv.append(t[i])

import pandas as pd
only_dv = pd.DataFrame()
only_dv['disease'] = disease
only_dv['dv'] = dv

only_dv.to_csv('./onlydv.csv',index=False)