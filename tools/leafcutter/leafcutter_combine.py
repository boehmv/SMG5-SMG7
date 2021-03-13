#!/usr/bin/python3
import pandas as pd
import numpy as np
import xlsxwriter

# importing sys module
import sys

# reading leafcutter cluster and effect file
with open(sys.argv[1], 'rb') as cluster, open(sys.argv[2], 'rb') as effect, open(sys.argv[3], 'wb') as outf:

    #generate tab delimited pandas csv reader
    A=pd.read_csv(cluster, delimiter='\t')
    B=pd.read_csv(effect, delimiter='\t')

    #split cluster information into seperate columns
    A['chr'], A['clusterID'] = zip(*A['cluster'].apply(lambda x: x.split(':', 1)))
    B['chr'], B['start'], B['end'], B['cluster'] = zip(*B['intron'].apply(lambda x: x.split(':', 3)))


    #drop unnecessary columns
    A = A.drop(columns=['cluster', 'chr'])
    B = B.drop(columns=['intron'])

    #merge cluster and effect data frames
    C = pd.merge(B, A, how='left',
        left_on='cluster', right_on='clusterID')
    C = C.drop(columns=['cluster'])  
    mygenes = C['genes'].str.split(" ", n = 1, expand = True)
    C['gene_name'] =  mygenes[0]
    C['gene_id'] = mygenes[1]
    C = C.drop(columns=['genes'])
    altgenes = C['gene_id'].str.split(",", expand = True)
    C['gene_id'] = altgenes[0]
    altgenes = altgenes.drop(altgenes.columns[0], 1)
    C['altgenes'] = altgenes.apply(lambda x: x.str.cat(sep=','), axis=1)
    

    #filter for p.adjust<0.05
    #D = C.query('p.adjust < 0.05')

    #write to Excel file
    writer = pd.ExcelWriter(outf, engine='xlsxwriter')
    C.to_excel(writer, sheet_name='total')
    #D.to_excel(writer, sheet_name='P.adjust < 0.05')
    writer.save()
