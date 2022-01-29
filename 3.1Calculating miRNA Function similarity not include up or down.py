import pandas as pd
import numpy as np
from scipy.linalg import norm
TreeAndDv = pd.read_csv(u'./server_sql_miRNA_tree_and_DV_up_down_ALL.txt',header=None,low_memory=False,sep='\t')
TreeAndDv[1] = TreeAndDv[1].apply(lambda x:eval(x))
TreeAndDv[2] = TreeAndDv[2].apply(lambda x:eval(x))

def cos_sim(vector_a, vector_b):
    """
    计算两个向量之间的余弦相似度
    :param vector_a: 向量 a
    :param vector_b: 向量 b
    :return: sim
    """
    vector_a = np.mat(vector_a)
    vector_b = np.mat(vector_b)
    num = float(vector_a * vector_b.T)
    denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
    if denom == 0:
        cos = 0
    else:
        cos = num / denom
    return cos

def get_similarity(dv1,dv2,tree1,tree2,theat):
    # 附加余弦相似性
    d1d2dict = {}
    for d1 in dv1:
        for d2 in dv2:
            # 计算两个疾病之间的相似性
            d1_set, d2_set = set([x for x in tree1[d1]]), set([x for x in tree2[d2]])
            d1_d2_set = d1_set & d2_set
            if len(d1_d2_set) == 0:
                name = sorted([d1,d2])
                name = '__'.join(name)
                inter = 0
            else:
                inter = 0
                for i in d1_d2_set:
                    inter = inter + tree1[d1][i] + tree2[d2][i]
                name = sorted([d1,d2])
                name = '__'.join(name)
            d1d2dict[name] = theat * inter / (dv1[d1] + dv2[d2])

    d1d2onehot = pd.DataFrame({'1': dict(dv1, ** d1d2dict), '2': dict(dv2,**d1d2dict)}).T.fillna(0).values
    A = np.array(d1d2onehot[0, :]).reshape(1, -1)[0]
    B = np.array(d1d2onehot[1, :]).reshape(1, -1)[0]
    return cos_sim(A, B),dict(dv1, ** d1d2dict),dict(dv2,**d1d2dict)


using = {}
result = []
mirName = []
ko = []
rbtmp = []
f2 = open('./server_onehot_DS_ALL.txt','w')
with open('./serverr_bmf_result_ALL.txt','w') as f:
    for i in TreeAndDv.values:
        mirName.append(i[0])
        tmp = []
        for j in TreeAndDv.values:
            if '{}_{}'.format(j[0],i[0]) in using:
                tmp.append(using['{}_{}'.format(j[0],i[0])])
            else:
                if i[0] == j[0]:
                    tmp.append(1.0)
                    using['{}_{}'.format(i[0],j[0])] = 1.0
                else:
                    sim_a_b,ra,rb = get_similarity(i[1],j[1],i[2],j[2],1)
                    tmp.append(sim_a_b)
                    using['{}_{}'.format(i[0], j[0])] = sim_a_b

                    if ra in rbtmp:
                        pass
                    else:
                        rbtmp.append(ra)
                        for k in ra:
                            if k in ko:
                                pass
                            else:
                                ko.append(k)
                                if ra[k] == 0:
                                    pass
                                else:
                                    f2.write(str(k) + '\t' + str(ra[k]))
                                    f2.write('\n')

                    if rb in rbtmp:
                        pass
                    else:
                        rbtmp.append(rb)
                        for k in rb:
                            if k in ko:
                                pass
                            else:
                                ko.append(k)
                                if rb[k]==0:
                                    pass
                                else:
                                    f2.write(str(k) + '\t' + str(rb[k]))
                                    f2.write('\n')
                    # f2.write(i[0] + '\t' + j[0] + '\t' + str(ra) + '\t' + str(rb) +  '\n')
        f.write('\t'.join(map(str,tmp)))
        f.write('\n')

with open('./server_miRNA_name1.txt','w') as f:
    for i in mirName:
        f.write(str(i) + '\n')
