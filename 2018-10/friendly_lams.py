import pandas as pd
from openpyxl import load_workbook


# First read in Excel file as DataFrame. The file should be in the format:
# Radial Poloidal z total EF180 EF182 EF183 EF184 EF185

filename = '/mnt/c/Users/Shawn/Google Drive/School/Tennessee/Research/' + \
           'Polodial_Scans/CD09_Map.xlsx'

df = pd.read_excel(filename, usecols=[0,1,2], sheet_name='Sheet1')

# Make the index the Radial location, grouped under poloidal locations.
df = df.pivot('Radial [mm]', 'Poloidal [mm]', 'total')

# Little bit of lines so we can write to the existing file and not overwrite it.
book = load_workbook(filename)
writer = pd.ExcelWriter(filename, engine='openpyxl')
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

df.to_excel(writer, sheet_name='Sheet2')
writer.save()
