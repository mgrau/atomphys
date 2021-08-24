import urllib
from bs4 import BeautifulSoup
import os
import json
import re

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)
nuc_table_filename = os.path.join(directory, "NuclearDataJSON.json")
with open(nuc_table_filename) as f:
    nuc_table = json.load(f)

def fetch_isotopes(name, symbol):
    url = 'https://en.wikipedia.org/wiki/Isotopes_of_' + name.lower()
    with urllib.request.urlopen(url) as response:
        response = response.read()
    soup = BeautifulSoup(response, 'lxml')

    tables = soup.find_all('table', {'class': 'wikitable'})
    last_table = tables[-1]

    col_headers = []
    for col in last_table.find_all('th'):
        try:
            col_headers.append(col.a.get('title'))
        except AttributeError:
            pass
    table_proto = dict.fromkeys(col_headers)
    col_headers[0] = "Nuclide"

    data = []
    colspans = []
    for row in last_table.find_all('tr')[1::1]:
        row_data = table_proto
        col_data = []
        colspan = []
        for col in row.find_all('td'):
            try:
                col_data.append(col.text.strip('\n'))
                if col.get('colspan') is not None:
                    colspan.append(int(col.get('colspan')))
                else:
                    colspan.append(1)
            except IndexError:
                pass
        data.append(col_data)
        colspans.append(colspan)

    # assume most freq number of total cols is the right one
    cols_count = [sum(i) for i in colspans]
    most_freq_num_cols = max(set(cols_count), key=cols_count.count)

    # only keep rows with right number of total cols
    rows_keep = [i for i, e in enumerate(colspans) if sum(e) == most_freq_num_cols]
    data = [data[i] for i in rows_keep]
    colspans = [colspans[i] for i in rows_keep]

    # weird indexing to span multiple columns in output table
    data_out = []
    for i in range(len(data)):
        col_names_cut = []
        col_data_cut = []
        data_index = 0
        j = 0
        cols_remaining_with_this_header = colspans[i][0]
        while j < most_freq_num_cols - 1:
            if symbol in data[i][0]:
                col_names_cut.append(col_headers[j])
                col_data_cut.append(data[i][data_index])
            j = j + 1
            cols_remaining_with_this_header = cols_remaining_with_this_header - 1
            if cols_remaining_with_this_header == 0:
                data_index = data_index + 1
                cols_remaining_with_this_header = colspans[i][data_index]
        result = dict(zip(col_names_cut, col_data_cut))
        if len(result) > 0:
            data_out.append(result)
        i = i + 1

    # clean up nuclide column so we can compare to nuc_table
    [[row.update({'Nuclide': re.search('^\d+[m]*\d*[A-Za-z]+',row['Nuclide']).group()})] for row in data_out]

    nuc_table_nuclides = [row['Nucleus'] for row in nuc_table]
    
    for row in range(len(data_out)):
        if data_out[row]['Nuclide'] in nuc_table_nuclides:
            indx = nuc_table_nuclides.index(data_out[row]['Nuclide'])
            data_out[row].update( {"mag_moment_μN": nuc_table[indx]['μ(nm)']})    
        else:
            data_out[row].update( {"mag_moment_μN": '0.0'})

    return data_out
