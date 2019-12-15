
"""
Created on Tue Sep 24 13:35:02 2019

@author: FernandoZeng
"""

import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
os.chdir('E:\\oc')
GPLname='GPL8179-tbl-1.txt'#platform
dataname='GSE69002_series_matrix.txt'#GES_data
datarows=63#data_begin
zcf='GSM1689916'#the first GSM(non-oscc)
zcl='GSM1689919'#the end GSM(non-oscc)
occsf='GSM1689913'#the first GSM(oscc)
occsl='GSM1689915'#the end GSM(oscc)
gsmdl=4
gsmend=8
mygenedna='ILMN_Gene'


data=pd.read_csv(dataname,delimiter='\t',
                 skiprows=datarows)
#data=data.drop(data.index[0:34])
data.rename(columns={'ID_REF':'gene_id'},inplace=True)
data.isna().sum()
data=data.dropna(axis=0)
data.dtypes
data.dtypes.eq(object)
cols=data.columns[data.dtypes.eq(object)]
cols=cols[1:]
data[cols]=data[cols].apply(pd.to_numeric,errors='coerce',axis=0)
data.dtypes.eq(object)
#plt.show(data.plot(kind = 'density', title = dataname+' Density'))
cancer = data.loc[:, occsf :occsl].mean(axis = 1)
no = data.loc[:,zcf:zcl].mean(axis = 1)
fold = np.log2(cancer/no)
'''
plt.hist(fold)
plt.title("Histogram of log2(FC)")
plt.show()
'''

from scipy import stats

pvalue = []
for i in range(0, len(data)):
    ttest = stats.ttest_ind(data.iloc[i][1:gsmdl], data.iloc[i][gsmdl:gsmend])
    pvalue.append(ttest[1])
    
myarray = np.asarray(pvalue)

'''
# Histogram of the p-values
plt.hist(-np.log(pvalue))
plt.title("Histogram of p-value")
plt.show()

'''
result = pd.DataFrame({'pvalue':myarray,'FoldChange':fold})
result['log(pvalue)'] = -np.log10(result['pvalue'])
result['log2Fc']=np.log2(result['FoldChange'])
result['sig'] = 'normal'
result['size']  =np.abs(result['FoldChange'])/10
result.loc[(result.FoldChange> 0.125 )&(result.pvalue < 0.05),'sig'] = 'up'
result.loc[(result.FoldChange< -0.125 )&(result.pvalue < 0.05),'sig'] = 'down'
'''
ax = sns.scatterplot(x="FoldChange", y="log(pvalue)",
                      hue='sig',
                      hue_order = ('down','normal','up'),
                      palette=("#377EB8","grey","#E41A1C"),
                      data=result)
ax.set_ylabel('-log(pvalue)',fontweight='bold')
ax.set_xlabel('log2(FC)',fontweight='bold')
plt.title(dataname[0:8])
'''

MyGene=pd.read_csv(GPLname,delimiter='\t')
ge=pd.DataFrame()
ge['ID']=data.iloc[:]['gene_id']
#ge['ID']=ge['ID'].apply(pd.to_numeric,errors='coerce')
ge['gename']=data.iloc[:]['gene_id']

j=0
test=[]
for i in range(0,len(MyGene)):
    if(ge.iloc[j]['ID']==MyGene.iloc[i]['ID']):    
        #ge.iloc[j]['gename']=MyGene.iloc[i][mygenedna]
        test.append(MyGene.iloc[i][mygenedna])
        j=j+1
        print(i)
        

result['gene_id']=data['gene_id']
result['gename']=test
#result.to_csv(dataname+'_re.csv')
'''
#2
for i in result[~result['sig'].isin(['normal'])]:
    result.loc[i,'gename']=MyGene[MyGene['ID'].isin([result.iloc[i]['gename']])][mygenedna]

for i in result[~result['sig'].isin(['normal'])]['gename']:
    result[result['gene_id'].isin([i])]['gename']=MyGene[MyGene['ID'].isin([i])][mygenedna]
'''
'''
#fold-change
wt = data.loc[:, occsf :occsl].mean(axis = 1)
ko = data.loc[:,zcf:zcl].mean(axis = 1)
fold = ko - wt
plt.hist(fold)
plt.title("Histogram of fold-change")
plt.show()
'''

'''
from scipy import stats

pvalue = []
for i in range(0, len(data)):
    ttest = stats.ttest_ind(data.iloc[i][1:gsmdl], data.iloc[i][gsmdl:gsmend])
    pvalue.append(ttest[1])

# Histogram of the p-values
plt.hist(-np.log(pvalue))
plt.title("Histogram of p-value")
plt.show()
'''

'''
#火山图
myarray = np.asarray(pvalue)
result = pd.DataFrame({'pvalue':myarray,'FoldChange':fold})
result['log(pvalue)'] = -np.log10(result['pvalue'])

result['sig'] = 'normal'
result['size']  =np.abs(result['FoldChange'])/10
result.loc[(result.FoldChange> 1 )&(result.pvalue < 0.05),'sig'] = 'up'
result.loc[(result.FoldChange< -1 )&(result.pvalue < 0.05),'sig'] = 'down'

ax = sns.scatterplot(x="FoldChange", y="log(pvalue)",
                      hue='sig',
                      hue_order = ('down','normal','up'),
                      palette=("#377EB8","grey","#E41A1C"),
                      data=result)
ax.set_ylabel('-log(pvalue)',fontweight='bold')
ax.set_xlabel('log2(FC)',fontweight='bold')
plt.title(dataname[0:8])
'''
'''
#筛选差异基因
fold_cutoff = 0.125
pvalue_cutoff = 0.5

filtered_ids = []
lfold=list(fold)
for i in range(0, len(data)):
    if (abs(lfold[i]) >= fold_cutoff) and (pvalue[i] <= pvalue_cutoff):
        filtered_ids.append(i)
        
filtered = data.iloc[filtered_ids]
print("Number of DE genes: ")
print(len(filtered.index))
del filtered['gene_id']
filtered.head()
'''
'''

#热图
sns.clustermap(filtered, cmap='RdYlGn_r', standard_scale = 0)
'''

'''
#基因通过NPL平台命名
MyGene=pd.read_csv('GPL21572-124634.txt',delimiter='\t')
ge=pd.DataFrame()
ge['ID']=data.iloc[:]['gene_id']
#ge['ID']=ge['ID'].apply(pd.to_numeric,errors='coerce')
ge['gename']=data.iloc[:]['gene_id']
j=0
for i in range(0,len(MyGene)):
    if(ge.iloc[j]['ID']==MyGene.iloc[i]['ID']):    
        ge['gename'][j:j+1]=MyGene.iloc[i]['Transcript ID(Array Design)']
        j=j+1
        print(i)
result['gene_id']=data['gene_id']
result['gename']=ge['gename']


'''

'''
#交集
total=[]
for i in range(0,max(len(gea),len(geb))):
    for j in range(0,min(len(gea),len(geb))):
        if(gea.iloc[i]['gename']==geb.iloc[j]['gename'] and gea.iloc[i]['sig']!='normal'):
            total.append(gea.iloc[j]['gename'])
            break
        
'''

'''
#交集2
total=[]
for i in range(0,len(gea[~gea['sig'].isin(['normal'])])):
    for j in range(0,len(geb[~geb['sig'].isin(['normal'])])):
        if(gea.iloc[i]['gename']==geb.iloc[j]['gename']):
            total.append(gea.iloc[i]['gename'])
            break
'''

'''
#交集3
gea=pd.read_csv('GSE45238_series_matrix.txt_re.csv')
geb=pd.read_csv('GSE69002_series_matrix.txt_re.csv')
gec=pd.read_csv('GSE98463_series_matrix.txt_re.csv')
a=gea[~gea['sig'].isin(['normal'])].iloc[:]['gename'].tolist()
b=geb[~geb['sig'].isin(['normal'])].iloc[:]['gename'].tolist()
c=gec[~gec['sig'].isin(['normal'])].iloc[:]['gename'].tolist()
total=[]
for i in a:
    if i in b:
        if i in c:
            total.append(i)
for i in b:
    if i in a:
        if i in c:
            total.append(i)
for i in c:
    if i in a:
        if i in b:
            total.append(i)
total=list(set(total))

'''

'''
df_np = np.array((a+b+c)) # 将3个list合并在一起，形成数组
df_count = pd.Series(df_np).value_counts() # 获取每个元素出现的次数
r2 = list(df_count[df_count == 1].index)  # 获取在3个list中只出现过1次的元素
r2.sort(reverse=False) # 对r2进行排序，方便观察结果
print('r2 -->', r2)  # 输出：r2 --> [-1, 1, 3, 4, 5, 6, 8]

'''