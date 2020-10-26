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
    df_final = pd.DataFrame(columns=['uniqueID', 'gene_name', 'gene_id', 'altgenes'])
    
    # Get all conditions    
    print('# of conditions:')
    print(len(cond))
    print('Conditions:')
    print(cond)

    # Loop over each condition in leafcutter/Combined folder
    for f in cond[1:]:
        myfilename=('{folder}/control_vs_{condition}_final.xlsx'.format(folder=myfolder, condition=f))
        df = pd.read_excel(myfilename)

        # Report data before filtering
        print('# rows before filtering:')
        print(len(df.index))     
   
        # Define uniqueID
        df["uniqueID"] = df['chr'].map(str) + ":" + df['start'].map(str) + "-" + df['end'].map(str) + ":" + df['gene_id'].map(str)

        #drop unnecessary columns
        df = df.drop(columns=['Unnamed: 0', 'chr', 'start', 'end', 'logef', 'status', 'loglr', 'df', 'p'])
        condition=f
        # Reorder columns, split coordinates and give p.adjust and log2FoldChange condition names
        df = df[['uniqueID', 'gene_name', 'gene_id', 'altgenes', 'deltapsi', 'p.adjust', 'control', condition]]
        df.columns = ['uniqueID', 'gene_name', 'gene_id', 'altgenes', "%s_dPSI" % condition, "%s_p.adjust" % condition, "PSI_control_%s" % condition, "PSI_%s" % condition,]
        df_filtered = df[(df["%s_p.adjust" % condition] < 0.05)]
        print('# rows after filtering')
        print(len(df_filtered.index))
        df_final = pd.merge(df_final,df_filtered,on=['uniqueID', 'gene_name', 'gene_id', 'altgenes'], how='outer')
        print('# rows of final df')
        print(len(df_final.index))
        print("Condition %s processed" % condition)

    # Get proper coordinates
    IDs = df_final["uniqueID"].str.split(":", n = 2, expand = True)
    Coords = IDs[1].str.split("-", expand = True)
    df_final["coordinates"] = IDs[0] +':'+ Coords[0] +'-'+ Coords[1]
    cols = df_final.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_final = df_final[cols]    

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
    worksheet.set_column(1, 1, 10)
    worksheet.set_column(2, 2, 15)
    worksheet.set_column(3, 3, 20)
    worksheet.set_column(4, 4, 10)
    for f in range(5,end_column-1):
        worksheet.conditional_format(0, f, end_row-1, f, {'type':'3_color_scale', 'min_color': "red", 'mid_color': "white", 'max_color': "green"})
    writer.save()
