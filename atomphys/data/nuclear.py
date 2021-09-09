import urllib
from bs4 import BeautifulSoup
import os
import json
import re

current_file = os.path.realpath(__file__)
directory = os.path.dirname(current_file)

# various data file locations
raw_nuc_table_file = os.path.join(directory, "raw_NuclearMomentDataJSON.json")
periodic_table = os.path.join(directory, "PeriodicTableJSON.json")
output_raw_nuc_data_file = os.path.join(directory, "NuclearMomentDataJSON.json")
nuc_ptable_file = os.path.join(directory, "NuclearPeriodicTableJSON.json")

# nuclear data file is huge so load it only once
nuc_table_file = os.path.join(directory, "NuclearMomentDataJSON.json")
with open(nuc_table_file) as f:
    nuc_table = json.load(f)

def build_nuclear_ptable():
    wiki_data = get_all_wiki_data()
    nuc_ptable_data = wiki_append_nuc_mag_moments(wiki_data)

    with open(nuc_ptable_file, 'w') as f:
        json.dump(nuc_ptable_data, f, indent=4, ensure_ascii=False)

def get_all_wiki_data():
    with open(periodic_table) as f:
        pt = json.load(f)
    wiki_nuc_data = []
    for element in pt['elements']:
        try:
            print('loading isotope data for %s...' % element['name'])
            new_wiki_data = get_wiki_isotope_data(element['name'],element['symbol'])
            wiki_nuc_data = wiki_nuc_data + new_wiki_data
            print('done')
        except:
            print('failed to get data for %s' % element['name'])
    return wiki_nuc_data

def get_wiki_isotope_data(name,symbol):
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
    spin_col = [i for i,s in enumerate(col_headers) if 'spin' in str(s).lower()]
    if len(spin_col) > 0:
        col_headers[spin_col[0]] = 'Spin'

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

    return data_out

def wiki_append_nuc_mag_moments(wiki_nuc_data):
    nuc_table_nuclides = [row['Nucleus'] for row in nuc_table]
    
    for row in range(len(wiki_nuc_data)):
        if wiki_nuc_data[row]['Nuclide'] in nuc_table_nuclides and 'Spin' in wiki_nuc_data[row].keys() and not('0' in wiki_nuc_data[row]['Spin']):
            indx = nuc_table_nuclides.index(wiki_nuc_data[row]['Nuclide'])
            wiki_nuc_data[row].update( {"mag_moment_μN": nuc_table[indx]['μ(nm)']})    
        else:
            wiki_nuc_data[row].update( {"mag_moment_μN": '0.0'})

    return wiki_nuc_data

def massage_nuclear_data():
    with open(raw_nuc_table_file) as f:
        NuclearMomentDataJSON = json.load(f)

    header_row = NuclearMomentDataJSON[0]['data'][0]
    col_headers = []
    for col in header_row:
        col_headers.append(col['text'])

    table_data = []
    p = 0
    for page in NuclearMomentDataJSON:
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

    table_data_after_cut = []
    for row in table_data:
        nuc_re = re.search('^(\d+) (\w+) (\d+)', row['Nucleus'])
        if nuc_re is not None:
            row['Nucleus'] = nuc_re.group(3) + nuc_re.group(2)
            row['Z'] = nuc_re.group(1)
            table_data_after_cut.append(row)

    with open(output_raw_nuc_data_file, 'w') as f:
        json.dump(table_data_after_cut, f, indent=4, ensure_ascii=False)
