# coding:utf-8
import pandas as pd
import numpy as np
import multiprocessing
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')
theat = 0.5

# 读取原始的结构化数据
mtrees = pd.read_csv('./mtrees2019.bin', sep=';', header=None, low_memory=False)
mtrees[0] = mtrees[0].apply(lambda x: str(x))


# 计算DV_values
# 修正版本
def DV_values(A, index, trss, usingset):

    A_tmp = A.split('.')
    A_tmp = '.'.join(A_tmp[0:index])
    usingset.append(A_tmp)
    if len(str(A).split('.')) == index:
        trss[A_tmp] = 1
        return 1
    else:
        if A_tmp in trss:
            if theat * DV_values(A, index + 1, trss,usingset) > trss[A_tmp]:
                trss[A_tmp] = theat * DV_values(A, index + 1, trss,usingset)
        else:
            trss[A_tmp] = theat * DV_values(A, index + 1, trss,usingset)
        return theat * DV_values(A, index + 1, trss,usingset)

def get_mir_dis_in_mtrees(mir_dis, mtree, mirname):
    print('miRNA {} '.format(mirname))
    res_ = {}
    tre_ = {}
    for mirdis in set(mir_dis):
        sum_ = {}
        code = mtree.values[mtree.values[:, 0] == mirdis][:, 1]
        trss = {}
        usingset = []
        DV_tree = {}
        DV = 1
        for i in code:
            tmp = str(i).split('.')
            for j in range(len(tmp), 0, -1):
                DV_values(i, j, trss, usingset)
            sum_[i] = trss
        # print(sum_)
        for s in sum_:
            dn = mtree.values[mtree.values[:, 1] == s][:, 0][0]
            DV_tree[dn] = 1
            for values in sum_[s]:
                # 根据编码找到疾病名称
                dn = mtree.values[mtree.values[:, 1] == values][:, 0][0]
                if dn in DV_tree:
                    pass
                else:
                    DV_tree[dn] = sum_[s][values]
                    DV = DV + sum_[s][values]

        res_[mirdis] = DV
        tre_[mirdis] = DV_tree
    return mirname, res_ ,tre_




# print(DV_values(['Lung Neoplasms'],'xxx',mtrees))
# a,b,c = get_mir_dis_in_mtrees(['Lung Neoplasms'],mtrees,'xxxx')
# print(a)
# print(b)
# print(c)


# exit()
# 数据格式 [miRNA {疾病集合}]
mirdis = pd.read_csv('./server_sql_DV_DATA_ALL.txt', sep='\t', header=None, low_memory=False)
mirdisValues = mirdis.groupby([0])[1].apply(np.array).reset_index().values


if __name__ == '__main__':
    f = open('./server_sql_miRNA_tree_and_DV_up_down_ALL.txt', 'w')

    pool = multiprocessing.Pool(processes=4)
    result = []
    for listParam in list(mirdisValues):
        result.append(pool.apply_async(get_mir_dis_in_mtrees, [listParam[1], mtrees, listParam[0]]))
    pool.close()
    pool.join()

    for res in result:
        t = res.get()
        a, b, c = t[0], t[1], t[2]
        f.write(str(a) + '\t' + str(b) + '\t' + str(c) + '\n')
    f.close()