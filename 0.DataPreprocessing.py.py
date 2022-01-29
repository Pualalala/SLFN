import pandas as pd

# 读取mtree数据
print('read mtree2019.bin data')
mtrees = pd.read_csv('./mtrees2019.bin',sep=';',header=None)
mtrees.columns =['DescriptorName','TreeAddress']

# 读取HMDDV3数据
# 选择具有上下调控作用的数据作为我们使用的数据集
print('read HMDDv3 date')
HMDDV3 = pd.read_excel('./HMDDv3.xlsx',sheet_name=u'miRNA-疾病关联数据')
HMDDV3[u'Up/Down\n（上下调分类，用于Comparison页面）'] = HMDDV3[u'Up/Down\n（上下调分类，用于Comparison页面）'].fillna(-1)
# HMDDV3 = HMDDV3[HMDDV3[u'Up/Down\n（上下调分类，用于Comparison页面）']!=-1]
# 只保留miRNA disease up/down
HMDDV3 = HMDDV3[[u'miRNA ID\n(请按TAM2标准检查miRNA名称是否符合规范)',u'disease',u'Up/Down\n（上下调分类，用于Comparison页面）','PMID']]
HMDDV3.columns = ['miRNA','DiseaseName','up_down','PMID']
HMDDV3 = HMDDV3.drop_duplicates()

# print(HMDDV3)

# exit()
# 读取扩展数据集
print('read date_ext')
data_ext = pd.read_excel('./HMDDv3.xlsx')
data_ext = data_ext[['TERM','MESH']]
data_ext.columns = ['DiseaseName','MESH']
data_ext = data_ext.dropna()
HMDDV3 = pd.merge(HMDDV3,data_ext,on=['DiseaseName'],how='inner',copy=False)

# 读取MeSH到疾病名称的映射数据
print('read mesh ID map')
meshToID = pd.read_csv('./meshToID.csv')

HMDDV3 = pd.merge(HMDDV3,meshToID,right_on=['DescriptorUI'],left_on=['MESH'],how='inner',copy=False)

HMDDV3 = HMDDV3[['miRNA','DescriptorName','up_down','MESH','PMID']]


HMDDV3['up_down'] = HMDDV3['up_down'].replace(-1,'N/A')
print(HMDDV3)

HMDDV3 = HMDDV3[['miRNA','MESH','PMID','up_down','DescriptorName']]
HMDDV3.columns = ['mirna','mesh','pmd','updown','disease']
HMDDV3.to_csv('./server_sql_mirnadiseasepmdupdownmesh.csv',index=False)


# 保留一份上下调控的
HMDDV3[HMDDV3['updown']!='N/A'].to_csv('./server_sql_updown.csv',index=False)

# 保留一份只有miRNAs 和疾病名称的
DV_DATA = HMDDV3[HMDDV3['updown']!='N/A'][['mirna','disease']].drop_duplicates()
DV_DATA.to_csv('./server_sql_DV_DATA.txt',index=False,header=None,sep='\t')

DV_DATA = HMDDV3[HMDDV3['updown']!='N/A'][['mirna','disease','updown']].drop_duplicates()
DV_DATA.to_csv('./server_sql_DV_DATA_up_down.txt',index=False,header=None,sep='\t')

DV_DATA = HMDDV3[['mirna','disease']].drop_duplicates()
DV_DATA.to_csv('./server_sql_DV_DATA_ALL.txt',index=False,header=None,sep='\t')

print(len(set(HMDDV3['disease'].values)))
# exit()
# HMDDV3.to_csv('./dataset.txt',index=False,header=None,sep='\t')
#
# print('disease',len(HMDDV3['DescriptorName'].unique()))
# print('miRNA',len(HMDDV3['miRNA'].unique()))
#
# DV_DATA = HMDDV3[['miRNA','DescriptorName']].drop_duplicates()
# DV_DATA.to_csv('./DV_DATA.txt',index=False,header=None,sep='\t')

