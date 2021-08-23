import json
import pandas as pd
import re

with open('atomphys/data/tabula-nucdata.json') as f:
    nucdata = json.load(f)

header_row = nucdata[0]['data'][0]
col_headers = []
for col in header_row:
    col_headers.append(col['text'])


table_data = []
p = 0
for page in nucdata:
    #print("page = %d" % p)
    r = 0
    page['data'].pop(0)
    for row in page['data']:
        #print("r = %d" % r)
        this_row_data = {key:None for key in col_headers}
        c = 0
        for col in row:
            #print("c = %d" % c)
            try:
                this_row_data[col_headers[c]] = col['text']
            except Exception:
                pass
                #print("skipping page %d row %d col %d" % (page,r,c))
            c = c+1
        if not(all([len(item)==0 for item in this_row_data.values()])):
            table_data.append(this_row_data)
        r = r+1
    p = p+1

df = pd.DataFrame(table_data)

new_nuc_string = []
r = 0
keep_rows = []
for nuc_string in df.Nucleus:
    nuc_re = re.search('^(\d+) (\w+) (\d+)',nuc_string)
    if nuc_re is not None:
        new_nuc_string.append(nuc_re.group(3) + nuc_re.group(2))
        keep_rows.append(r)
    r = r+1

df = df.iloc[keep_rows]
df['Nucleus'] = new_nuc_string

df.to_csv(r'atomphys/data/nucdata.csv',index=False)