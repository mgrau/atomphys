import urllib
import csv
import io


def fetch_states(atom):
    url = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
    values = {
        'spectrum': atom,
        'units': 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        'format': 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        'multiplet_ordered': 1,  # energy ordred
        'term_out': 'on',  # output the term symbol string
        'conf_out': 'on',  # output the configutation string
        'level_out': 'on',  # output the energy level
        'unc_out': 0,  # uncertainty on energy
        'j_out': 'on',  # output the J level
        'g_out': 'on',  # output the g-factor
        'lande_out': 'off'  # output experimentally measured g-factor
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + '?' + get_postfix) as response:
        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()),
                          dialect='excel-tab', restkey='None')

    return data


def fetch_transitions(atom):
    # the NIST url and GET options.
    url = 'http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
    values = {
        'spectra': atom,
        'format': 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        'en_unit': 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        'line_out': 2,  # only with {1: transition , 2: level classifications}
        'show_av': 5,
        'allowed_out': 1,
        'forbid_out': 1,
        'enrg_out': 'on'
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + '?' + get_postfix) as response:
        # when there are no transitions ASD returns a texl/html page with the
        # error message "No lines are available in ASD with the parameters selected"
        # rather than the expected text/plain when using format=3
        if response.headers.get_content_type() != 'text/plain':
            print(response.headers.get_content_type())
            return []

        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()), dialect='excel-tab')
    data_with_transition_probabilities = [
        transition for transition in data if transition['Aki(s^-1)']]

    return data_with_transition_probabilities
