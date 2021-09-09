import json
import re
import os

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)
NuclearDataJSON_file = os.path.join(directory, "raw_NuclearDataJSON.json")

with open(NuclearDataJSON_file) as f:
    NuclearDataJSON = json.load(f)

header_row = NuclearDataJSON[0]['data'][0]
col_headers = []
for col in header_row:
    col_headers.append(col['text'])

table_data = []
p = 0
for page in NuclearDataJSON:
    r = 0
    page['data'].pop(0)
    for row in page['data']:
        this_row_data = {key: None for key in col_headers}
        c = 0
        for col in row:
            try:
                this_row_data[col_headers[c]] = col['text']
            except Exception:
                pass
            c = c + 1
        if not (all([len(item) == 0 for item in this_row_data.values()])):
            table_data.append(this_row_data)
        r = r + 1
    p = p + 1

table_data_cut = []
for row in table_data:
    nuc_re = re.search('^(\d+) (\w+) (\d+)', row['Nucleus'])
    if nuc_re is not None:
        row['Nucleus'] = nuc_re.group(3) + nuc_re.group(2)
        row['Z'] = nuc_re.group(1)
        table_data_cut.append(row)

output_data_file = os.path.join(directory, "NuclearDataJSON.json")
with open(output_data_file, 'w') as f:
    json.dump(table_data_cut, f)
