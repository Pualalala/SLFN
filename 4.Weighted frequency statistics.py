import pandas as pd

HMDD_PROCESSON = pd.read_csv('./server_sql_mirnadiseasepmdupdownmesh.csv',header=None)
HMDD_PROCESSON.columns = ['mirna','uniqueid','pmid','updown','disease']

miRNA_disease = HMDD_PROCESSON.groupby(['mirna'])['disease'].apply(list).reset_index()

miRNA_disease_updown = HMDD_PROCESSON[HMDD_PROCESSON['updown'].isin(['Up','Down'])].groupby(['mirna'])['disease'].apply(list).reset_index()

# mir1_mir2_similarity_Only_weight_DV_1.0.txt

# mir1_mir2_similarity_updown_weight_DV_1.5.txt

mir1_mir2_similarity_Only_weight_DV = pd.read_csv('./serverr_bmf_result_ALL.txt',sep='\t',header=None)
miRNA_NAME_ONLY = pd.read_csv('./server_miRNA_name1.txt',header=None)
mir1_mir2_similarity_Only_weight_DV.columns = miRNA_NAME_ONLY[0].values
mir1_mir2_similarity_Only_weight_DV.index = miRNA_NAME_ONLY[0].values
# print(miRNA_NAME_ONLY[0])
# print(mir1_mir2_similarity_Only_weight_DV)

mir1_mir2_similarity_updown_weight_DV = pd.read_csv('./rbmf_result_UP_DOWN.txt',sep='\t',header=None)
miRNA_NAME_UPDOWN = pd.read_csv('./miRNA_name.txt',header=None)
mir1_mir2_similarity_updown_weight_DV.columns = miRNA_NAME_UPDOWN[0].values
mir1_mir2_similarity_updown_weight_DV.index = miRNA_NAME_UPDOWN[0].values

# print(mir1_mir2_similarity_updown_weight_DV)
onlymirna = []
onlydisease = []
onlycount = []
from tqdm import tqdm
for i in tqdm(miRNA_NAME_ONLY.values):
    # print('---------------------------------')
    mirna = i[0]
    # print(mirna)

    # exit()
    # print(mirna)
    diseaseDict = {}
    tmp = mir1_mir2_similarity_Only_weight_DV[mirna]

    for ik,io in zip(tmp.index,tmp):
        ik_disease = miRNA_disease[miRNA_disease['mirna']==mirna]['disease'].values[0]
        if ik!=mirna and io!=0:
            # print(ik,io)
            k1tmp = miRNA_disease[miRNA_disease['mirna']==ik]['disease'].values[0]
            ktmp = []
            for x in k1tmp:
                if x not in ik_disease:
                    ktmp.append(x)
            for uk in list(set(ktmp)):
                if uk in diseaseDict.keys():
                    diseaseDict[uk] = diseaseDict.get(uk) + 1 * io
                else:
                    diseaseDict[uk] = 1 * io
        else:
            pass
            # print(ik,mirna)
    for wxz in diseaseDict:
        onlymirna.append(str(mirna).replace('\r\n', ''))
        onlydisease.append(wxz)
        onlycount.append(diseaseDict[wxz])

onlyDf = pd.DataFrame({'mirna':onlymirna,'disease':onlydisease,'count':onlycount})
onlyDf.to_csv('./onlyDf_.csv',index=False)

#####################################################################################
# print(mir1_mir2_similarity_updown_weight_DV)
onlymirna = []
onlydisease = []
onlycount = []
from tqdm import tqdm
for i in tqdm(miRNA_NAME_UPDOWN.values):
    # print('---------------------------------')
    mirna = i[0]
    # print(mirna)
    diseaseDict = {}
    tmp = mir1_mir2_similarity_updown_weight_DV[mirna]

    for ik,io in zip(tmp.index,tmp):
        ik_disease = miRNA_disease_updown[miRNA_disease_updown['mirna'] == mirna]['disease'].values[0]
        if ik!=mirna and io!=0:
            # print(ik,io)
            k1tmp = miRNA_disease_updown[miRNA_disease_updown['mirna']==ik]['disease'].values[0]
            ktmp = []
            for x in k1tmp:
                if x not in ik_disease:
                    ktmp.append(x)
            for uk in list(set(ktmp)):
                if uk in diseaseDict.keys():
                    diseaseDict[uk] = diseaseDict.get(uk) + 1 * io
                else:
                    diseaseDict[uk] = 1 * io
        else:
            pass
            # print(ik,mirna)
    for wxz in diseaseDict:
        onlymirna.append(str(mirna).replace('\r\n',''))
        onlydisease.append(wxz)
        onlycount.append(diseaseDict[wxz])

onlyDf = pd.DataFrame({'mirna':onlymirna,'disease':onlydisease,'count':onlycount})
onlyDf.to_csv('./onlyDf_updown.csv',index=False)