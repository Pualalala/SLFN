from tqdm import tqdm
disease = []
dv = []
tmps = []
with open('./rbmf_onehot_DS_UP_DOWN.txt') as f:
    for i in tqdm(f):
        tmp = str(i).split('\t')
        t = eval(tmp[2])
        if t in tmps:
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
updownsdv = pd.DataFrame()
updownsdv['disease'] = disease
updownsdv['dv'] = dv

updownsdv.to_csv('./updownsdv.csv',index=False)
