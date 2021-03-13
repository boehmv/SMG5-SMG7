#!/usr/bin/python3
import pandas as pd
import numpy as np
import xlsxwriter
import glob
from functools import reduce

def get_col_widths(dataframe):
    # First we find the maximum length of the index column   
    idx_max = max([len(str(s)) for s in dataframe.index.values] + [len(str(dataframe.index.name))])
    # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
    return [idx_max] + [max([len(str(s)) for s in dataframe[col].values] + [len(col)]) for col in dataframe.columns]

# importing sys module
import sys

# Import first argument as IRFinder/Combined folder
myfolder = sys.argv[1]

# Import third argument as condition string, gets converted to characters
cond=sys.argv[3]
cond=list(cond.split(","))

# reading second argument as final output file 
with open(sys.argv[2], 'wb') as outf:

    # Define final data frame as reference to align all data to
    df_final = pd.DataFrame(columns=['coordinates', 'strand', 'gene_ID', 'symbol', 'baseMean'])

    # Get all conditions    
    print('# of conditions:')
    print(len(cond))
    print('Conditions:')
    print(cond)

    # Loop over each condition in IRFinder/combined folder
    for f in cond[1:]:
        print(f)
        myfilename=('{folder}/{condition}_vs_control_IRFinder_results.csv'.format(folder=myfolder, condition=f))
        df = pd.read_csv(myfilename)

        # Report data before filtering
        print('# rows before filtering:')
        print(len(df.index))     
   
        # Drop unnecessary columns
        df = df.drop(columns=['lfcSE', 'stat', 'pvalue', 'status'])
        condition=f

        # Reorder columns, split coordinates and give p.adjust and log2FoldChange condition names
        df['chr'], df['coordinates'], df['strand'] = df['coordinates'].str.split(':', 2).str
        df['coordinates'] = df['chr'] + ':' + df['coordinates']
        df = df[['coordinates', 'strand', 'gene_ID', 'symbol', 'baseMean', 'log2FoldChange', 'padj']]
        df.columns = ['coordinates', 'strand', 'gene_ID', 'symbol', 'baseMean', "%s_log2FoldChange" % condition, "%s_p.adjust" % condition]

        writer = pd.ExcelWriter(outf, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='total', index=False)
        workbook  = writer.book
        worksheet = writer.sheets['total']
        writer.save()

        # Filter the dataframe for padjust
        df_filtered = df[(df["%s_p.adjust" % condition] < 0.5)]

        print(list(df_filtered.columns))

        # Report data after filtering
        print('# rows after filtering')
        print(len(df_filtered.index))

        # Do the real merge: Keep integral data (coordinates, etc.) on left side
        df_final = pd.merge(df_final,df_filtered,on=['coordinates', 'strand', 'gene_ID', 'symbol', 'baseMean'], how='outer')

        # Report total number of introns in the final file
        print('# rows of final df')
        print(len(df_final.index))
        print("Condition %s processed" % condition)
       
    # Get table range
    end_row = len(df_final.index)
    end_column = len(df_final.columns)-1
    cell_range = xlsxwriter.utility.xl_range(0, 0, end_row, end_column)

    #write to Excel file
    writer = pd.ExcelWriter(outf, engine='xlsxwriter')
    df_final.to_excel(writer, sheet_name='total', index=False)
    workbook  = writer.book
    worksheet = writer.sheets['total']

    # Hack for preserving column headers when inserting table
    header = [{'header': di} for di in df_final.columns.tolist()]
    worksheet.add_table(cell_range,{'header_row': True,'columns':header})

    # Formating the output excel file
    worksheet.set_zoom(100)
    for i, width in enumerate(get_col_widths(df_final)):
        worksheet.set_column(i, i, width)
    worksheet.set_column(0, 0, 25)
    worksheet.set_column(1, 1, 8)
    worksheet.set_column(2, 2, 20)
    worksheet.set_column(3, 3, 20)
    for f in range(4,end_column-1):
        worksheet.conditional_format(0, f, end_row-1, f, {'type':'3_color_scale', 'min_color': "red", 'mid_color': "white", 'max_color': "green"})
    writer.save()
