import urllib.request
import csv
import io
import re
import math
from fractions import Fraction


def get_nist_states(atom):
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

    data = csv.DictReader(io.StringIO(response.decode()), dialect='excel-tab')

    levels = [row for row in data]
    for level in levels:
        del level[None]

    return levels


def parse_nist_states(atom):

    levels = get_nist_states(atom)

    for level in levels:

        level['configuration'] = level['Configuration']
        del level['Configuration']

        level['term'] = level['Term']
        del level['Term']

        try:
            level['g'] = float(Fraction(level['g']))
        except:
            del level['g']

        level['energy'] = float(level['Level (Ry)'].strip('[]a+?'))
        del level['Level (Ry)']

        try:
            level['J'] = float(Fraction(level['J'].strip('?')))
        except:
            del level['J']

    levels = [{**level, **parse_term(level['term'])} for level in levels]

    return levels


LS_term = re.compile(r'^(?P<S>\d+)(?P<L>[A-Z])\*?$')
JJ_term = re.compile(r'^\((?P<J1>\d+/?\d*),(?P<J2>\d+/?\d*)\)\*?$')
SK_term = re.compile(r'^(?P<S>\d+)\[(?P<K>\d+/?\d*)\]\*?$')

L = {
    'S': 0,
    'P': 1,
    'D': 2,
    'F': 3,
    'G': 4,
    'H': 5,
    'I': 6,
    'K': 7,
    'L': 8,
    'M': 9,
    'N': 10,
    'O': 11,
    'Q': 12,
    'R': 13,
    'T': 14,
    'U': 15,
    'V': 16
}


def parse_term(term):
    '''
    parse term symbol string
    '''
    if term == 'Limit':
        return {}

    parity = -1 if '*' in term else 1

    match = LS_term.match(term)
    if match is None:
        match = JJ_term.match(term)
    if match is None:
        match = SK_term.match(term)
    if match is None:
        return {'partiy': parity}

    def convert(key, value):
        if key == 'S':
            return int(value)
        if key == 'J1' or key == 'J2' or key == 'K':
            return float(Fraction(value))
        if key == 'L':
            return L[value]

    term = {key: convert(key, value)
            for key, value in match.groupdict().items()}

    return {**term, 'parity': parity}


def get_nist_transitions(atom):
    # the NIST url and GET options.
    url = 'http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
    values = {
        'spectra': atom,
        'format': 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        'en_unit': 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        'line_out': 1,  # only with transition probabilities
        'show_av': 5,
        'allowed_out': 1,
        'forbid_out': 1,
        'enrg_out': 'on'
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + '?' + get_postfix) as response:
        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()), dialect='excel-tab')

    transitions = [row for row in data]

    for transition in transitions:
        del transition['']

    return transitions


def parse_nist_transitions(atom):

    transitions = get_nist_transitions(atom)

    for transition in transitions:
        del transition['Acc']

        if transition['Type'] != '':
            transition['type'] = transition['Type']
        del transition['Type']

        transition['Gamma'] = float(transition['Aki(s^-1)'])
        del transition['Aki(s^-1)']

        transition['Ei'] = float(transition['Ei(Ry)'].strip('[]a +?'))
        del transition['Ei(Ry)']

        transition['Ek'] = float(transition['Ek(Ry)'].strip('[]a +?'))
        del transition['Ek(Ry)']

    return transitions
