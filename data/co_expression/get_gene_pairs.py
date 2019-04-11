import pandas as pd
import numpy as np
import os
fold_path='/Share/home/xixiaochen/project/training/adjp/'
file_li=os.listdir(fold_path)
adjp_fi_li=[]
for i in range(len(file_li)):
    adjp_pre_li=[]
    file_name=str(i+1)+'.txt'
    file_path=fold_path+file_name
    #print file_path
    data=pd.read_csv(file_path,sep='\t')
    data=data.set_index('order_index').sort_index()
    adjp=data['adjp']
    adjp_pre_li=list(adjp)
    adjp_fi_li.append(adjp_pre_li)
fi_adjp=pd.DataFrame(adjp_fi_li)
adjp_T=fi_adjp.T
comp_adjp=adjp_T[adjp_T<0.01]
#add gene names, each row represent a protein-coding gene, each column represent a lncRNA
col_file='/Share/home/xixiaochen/project/training/Pcc_pvalue_colnames.txt'
row_file='/Share/home/xixiaochen/project/training/Pcc_pvalue_rownames.txt'
col=pd.read_csv(col_file)
colnames=list(col['x'])[:5196]
row=pd.read_csv(row_file)
rownames=list(row['x'])[5196:]
def com_get_pair(df):
    df.columns=colnames
    df.index=rownames
    a=[]
    for i in range(df.index.size):
        for j in range(df.columns.size):
            b=[]
            if not np.isnan(df.iloc[i][j]):
                #add the row names
                b.append(df.iloc[i].name)
                #add the column names
                b.append(df.iloc[:,j].name)
                #add value
                b.append(df.iloc[i][j])
                a.append(b)
    df_df=pd.DataFrame(a)
    df_df.columns=['coding genes','lncRNA','Pcc']
    return df_df
#filter Pcc in top 5%
pcc_file='/Share/home/xixiaochen/project/training/lnc_coding_pcc.txt'
pcc=pd.read_csv(pcc_file,sep='\t',header=None)
#filter Pcc in bottom 5%
pcc=abs(pcc)
all_pcc=[]
for i in range(pcc.columns.size):
    all_pcc.extend(list(pcc.iloc[:,i]))
all_pcc_se=pd.Series(all_pcc)
cut_pcc=all_pcc_se.quantile(0.95)
pcc_cut=pcc[pcc>cut_pcc]
pcc_df=com_get_pair(pcc_cut)
adjp_df=com_get_pair(comp_adjp)
pcc_out_file='/Share/home/xixiaochen/project/training/pcc_cut.txt'
pcc_df.to_csv(pcc_out_file,sep='\t',index=False)
adjp_out_file='/Share/home/xixiaochen/project/training/adjp_cut.txt'
adjp_df.to_csv(adjp_out_file,sep='\t',index=False)

