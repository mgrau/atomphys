# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 12:32:11 2016

@author: Matt Grau
"""

import urllib.request
from bs4 import BeautifulSoup
import re
import fractions

import numpy as np
import pandas as pd
import scipy.constants as c_
import matplotlib.pyplot as plt
from sympy.physics import wigner


class Atom:
    """An atom class. 
    # Detailed description
    You can initialize atom with your favourit atomic species. You can then download the spectral data for that atom from the NIST atomic spectra database.
    After you have done that once, the data is saved locally in the folder parsed_data/ from which it is loaded everytime you reinitialize the class for that atom

    Example for neutral magnesium
    ```python
    atom = Atom('mg')
    atom.populate_states() # parse data from NIST database and store it locally. This has to be done only once for each atom.
    ```

    """

    def __init__(self, atom, data_folder='/home/tiqi/Projects/OpticalTrap/Calculations/Software library/Python/atom/parsed_data/'):
        """Constructor. 
        The atom name `atom` is what is used for retreiving the data from the NIST data base.
        For neutral you can use either mg i, mg 0, mg
        For singly ionized you can use mg ii, mg 1, mg+
        And so on.
        """

        self.name = atom
        self.states = None
        self.data_folder = data_folder

        try: 
            self.load_states()
        except IOError:
            pass

    def __repr__(self):
        rstr = "Atom: " + self.name + '\n'
        if self.states is not None:
            rstr += "%d states found" % len(self.states)
        return rstr

    def __call__(self, index=0):
        return self.states.iloc[index]

    def term(self, index=0, include_n=False):
        state = self.states.iloc[index]

        if include_n:
            return '$%d^{%d}%s_{%s}$' % (state.n,
                                         2*state.S+1,
                                         L_str(state.L),
                                         fractions.Fraction(state.J))
        else:
            return '$^{%d}%s_{%s}$' % (2*state.S+1,
                                       L_str(state.L),
                                       fractions.Fraction(state.J))

    def lifetime(self, state_index):
        transitions = self.transitions_down(state_index)
        if transitions.empty:
            return -1
        else:
            return 1/sum(transitions.gamma)

    def transitions_up(self, state_index):
        state = self.states.iloc[state_index]
        transitions = self.transitions[self.transitions.E_i == state.E]
        transitions = transitions.sort_values('E_k')
        transitions.reset_index(drop=True, inplace=True)
        return transitions

    def transitions_down(self, state_index):
        state = self.states.iloc[state_index]
        transitions = self.transitions[self.transitions.E_k == state.E]
        transitions = transitions.sort_values('E_i')
        transitions.reset_index(drop=True, inplace=True)
        return transitions

    def ac_stark(self, state_index, omega, I=1):
        state = self.states.iloc[state_index]
        alpha = self.scalar_polarizability(state, omega)
        dE = -alpha/(2*c_.c*c_.epsilon_0)*I
        return dE

    def r_sc(self, state_index, omega, I=1):
        state = self.states.iloc[state_index]
        alpha = self.scalar_polarizability(state, omega)
        R = omega**3*alpha**2/(6*c_.hbar*c_.pi*c_.epsilon_0**2*c_.c**4)*I
        return R

    def scalar_polarizability(self, state, omega):
        t_up = self.transitions[self.transitions.E_i == state.E]
        omega0 = (t_up.E_k - t_up.E_i)*c_.e/c_.hbar
        RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*(2*t_up.J_k+1)/(2*t_up.J_i+1)*t_up.gamma
        RME2 = np.tile(RME2, (np.size(omega), 1))
        omega0_m = np.tile(omega0, (np.size(omega), 1))
        omega_m = np.tile(omega, (len(omega0), 1)).T
        det = omega0_m/(omega0_m**2 - omega_m**2)
        alpha_up = np.sum((2/3)*det*RME2, axis=1)

        t_down = self.transitions[self.transitions.E_k == state.E]
        omega0 = (t_down.E_k - t_down.E_i)*c_.e/c_.hbar
        RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*t_down.gamma
        RME2 = np.tile(RME2, (np.size(omega), 1))
        omega0_m = np.tile(omega0, (np.size(omega), 1))
        omega_m = np.tile(omega, (len(omega0), 1)).T
        det = omega0_m/(omega0_m**2 - omega_m**2)
        alpha_down = np.sum((2/3)*det*RME2, axis=1)

        return alpha_up - alpha_down

    def vector_polarizability(self, state, omega):
        J = state.J

        w6j = np.vectorize(lambda Jp: wigner.wigner_6j(1, 1, 1, J, J, Jp))
        pre = np.sqrt(6*J*(2*J+1)/(J+1))*(2*J+1)

        t_up = self.transitions[self.transitions.E_i == state.E]
        if len(t_up) > 0:
            omega0 = (t_up.E_k - t_up.E_i)*c_.e/c_.hbar
            RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*(2*t_up.J_k+1)/(2*t_up.J_i+1)*t_up.gamma
            RME2 = np.tile(RME2, (np.size(omega), 1))
            omega0_m = np.tile(omega0, (np.size(omega), 1))
            omega_m = np.tile(omega, (len(omega0), 1)).T
            det = omega0_m/(omega0_m**2 - omega_m**2)
            decJ = pre*(-1)**(-J-t_up.J_k)*w6j(t_up.J_k)
            decJ = np.tile(decJ, (np.size(omega), 1))
            alpha_up = np.sum(det*RME2*decJ, axis=1)
        else:
            alpha_up = 0

        t_down = self.transitions[self.transitions.E_k == state.E]
        if len(t_down) > 0:
            omega0 = (t_down.E_k - t_down.E_i)*c_.e/c_.hbar
            RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*t_down.gamma
            RME2 = np.tile(RME2, (np.size(omega), 1))
            omega0_m = np.tile(omega0, (np.size(omega), 1))
            omega_m = np.tile(omega, (len(omega0), 1)).T
            det = omega0_m/(omega0_m**2 - omega_m**2)
            decJ = pre*(-1)**(-J-t_down.J_k)*w6j(t_down.J_k)
            decJ = np.tile(decJ, (np.size(omega), 1))
            alpha_down = np.sum(det*RME2*decJ, axis=1)
        else:
            alpha_down = 0

        return alpha_up - alpha_down

    def tensor_polarizability(self, state, omega):
        J = state.J

        w6j = np.vectorize(lambda Jp: wigner.wigner_6j(1, 1, 2, J, J, Jp))
        pre = np.sqrt(40*J*(2*J+1)*(2*J-1)/(3*(J+1)*(2*J+3)))

        t_up = self.transitions[self.transitions.E_i == state.E]
        if len(t_up) > 0:
            omega0 = (t_up.E_k - t_up.E_i)*c_.e/c_.hbar
            RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*(2*t_up.J_k+1)/(2*t_up.J_i+1)*t_up.gamma
            RME2 = np.tile(RME2, (np.size(omega), 1))
            omega0_m = np.tile(omega0, (np.size(omega), 1))
            omega_m = np.tile(omega, (len(omega0), 1)).T
            det = omega0_m/(omega0_m**2 - omega_m**2)
            decJ = pre*(-1)**(-J-t_up.J_k)*w6j(t_up.J_k)
            decJ = np.tile(decJ, (np.size(omega), 1))
            alpha_up = np.sum(det*RME2*decJ, axis=1)
        else:
            alpha_up = 0

        t_down = self.transitions[self.transitions.E_k == state.E]
        if len(t_down) > 0:
            omega0 = (t_down.E_k - t_down.E_i)*c_.e/c_.hbar
            RME2 = 3*c_.pi*c_.epsilon_0*c_.c**3/omega0**3*t_down.gamma
            RME2 = np.tile(RME2, (np.size(omega), 1))
            omega0_m = np.tile(omega0, (np.size(omega), 1))
            omega_m = np.tile(omega, (len(omega0), 1)).T
            det = omega0_m/(omega0_m**2 - omega_m**2)
            decJ = pre*(-1)**(-J-t_down.J_k)*w6j(t_down.J_k)
            decJ = np.tile(decJ, (np.size(omega), 1))
            alpha_down = np.sum(det*RME2*decJ, axis=1)
        else:
            alpha_down = 0

        return alpha_up - alpha_down

    def set_I(self, I=0):
        self.states['I'] = I
        new_states = pd.DataFrame([])
        for i in range(len(self.states)):
            state = a.states.iloc[i]
            F_min = abs(state.J - state.I)
            F_max = abs(state.J + state.I)
            for F in np.arange(F_min, F_max+1):
                new_states = new_states.append(state, ignore_index=True)
                new_states.loc[new_states.tail(1).index,'F'] = F
        self.states = new_states

    def plot_states(self, num_states=0, min_spacing=0.0, min_gamma=1e3,
                    show_transitions=True, show_lifetimes=False,
                    show_energy=False, include_n=False, label_size=12):
        states = self.states
        if num_states == 0:
            num_states = len(states)
        if type(num_states) is list:
            states = states.iloc[num_states]
        elif type(num_states) is range:
            states = states.iloc[list(num_states)]
        else:
            states = states.iloc[range(num_states)]

        line_width = 0.75

        energy = np.array(states.E)
        max_energy = max(energy)
        energy /= max_energy
        while any(np.diff(energy) < min_spacing):
            fix_index = np.nonzero(np.diff(energy) < min_spacing)[0][0]
            energy[fix_index+1:] += min_spacing*1.01 - np.diff(energy[fix_index:fix_index+2])
        energy *= max_energy

        for i in range(len(states)):
            state = states.iloc[i]
            plt.hlines(energy[i], state.L-line_width/2, state.L+line_width/2)
            label = self.term(states.index[i], include_n)
            if show_lifetimes:
                label += ',  $\\tau = %ss$' % to_SI(self.lifetime(i))
            t = plt.text(state.L+line_width/2, energy[i], label,
                         size=label_size, verticalalignment = 'center')

        if show_transitions:
            for i in range(len(states)):
                state = states.iloc[i]
                # get all transitions emanating from state i
                transitions = self.transitions_up(states.index[i])
                transitions = transitions[transitions.gamma > min_gamma]
                # select only those transitions that go to states we have plotted
                if transitions.empty == False:
                    transitions = transitions[[E_k in np.array(states.E) for E_k in transitions.E_k]]

                if transitions.empty == False:
                    angle = np.arctan2(transitions.E_k - transitions.E_i, transitions.L_k - transitions.L_i)
                    transitions = transitions.iloc[np.argsort(angle)[::-1]]
                    xpos = np.linspace(state.L-line_width/3, state.L+line_width/3, len(transitions)+2)[1:-1]
                    for k in range(len(transitions)):
                        transition = transitions.iloc[k]
                        upper_state = states[states.E == transition.E_k]
                        upper_state = upper_state.head(1)
                        x0 = float(xpos[k])
                        y0 = float(energy[i])
                        x1 = float(upper_state.L)
                        en = energy[np.where(states.E ==transition.E_k)]
                        y1 = float(en[0])
                        ax = plt.gca()
                        s0 = ax.transData.transform_point((x0, y0))
                        s1 = ax.transData.transform_point((x1, y1))
                        angle = (180/np.pi)*np.arctan2(s1[1]-s0[1], s1[0]-s0[0])
                        if angle>90:
                            angle -= 180
                        plt.annotate('',
                                     xy=(x0, y0), xytext=(x1, y1),
                                     arrowprops={'arrowstyle': '<|-|>', 'lw': 0.5})
                        label = '$\\lambda = %s m$\n' % to_SI(transition.wl)
                        label += '$\\Gamma = 2\\pi \\times %sHz$' % to_SI(transition.gamma/(2*np.pi))
                        t = plt.annotate(label,
                                         xy=((x0+x1)/2, (y0+y1)/2),
                                         xytext=((x0+x1)/2, (y0+y1)/2),
                                         size=8,
                                         rotation=angle,
                                         horizontalalignment='center',
                                         verticalalignment='center')
                        t.draggable()

        if show_energy:
            plt.ylabel('Energy (eV)')

        plt.ylim(plt.ylim()[0] - 0.05*np.diff(plt.ylim()),
                 plt.ylim()[1] + 0.05*np.diff(plt.ylim()))

        plt.xticks([])
        if ~show_energy:
            plt.yticks([])

        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)

        return [(states.L.iloc[i], energy[i]) for i in range(len(states))]

    def get_nist_transitions(self, atom):
        """Helper function used in `populate_states` that extracts the atomic state data for a certain atom from the NIST atomic spectra database.

        # Detailed description
        Function used by `populate_states` to get atomic data. Sends a http request and parses the response. The data is contained in a <table> html contaier. 
        We can select it by finding the one table htat has the color #FFFEEE in it's style.

        # Arguments
        * `atom` -- Atom name string passed to NIST data base. You can use mg as shorthand for mg i (neutral magnesium) and mg+ as shorthand for mg ii. 
        
        # Returns
        * `data` -- a pandas dataframe containing all transitions.
        """

        # if the atom string doesn't match the form 'atom+', the shorthand for
        # a singly charge ion, check to see if a charge is specified. If not,
        # we will assume the neutral (e.g. 'atom i').
        # mg -> mg i
        # mg+ -> mg ii
        # mg++ -> mg iii

        if re.search(' [i, I, 0-9]{1,999}$', atom) is None: # check whether there is at least one i, I or a number at the end of the atom name
            if re.search('\+{1,999}$', atom) is not None: # check whether atom+ was used as a shorthand for atom ii
                s = re.search('\+{1,999}$', atom) # convert atom+ into atom i
                atom = atom[:s.start()] + ' ' + ''.join(['i' for i in range(s.end()-s.start())])
            else: # mg -> mg i
                atom += ' i'     

        # the NIST url and GET options.
        url = 'http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
        values = {
            'spectra': atom,
            'low_wl': '',
            'upp_wn': '',
            'upp_wl': '',
            'low_wn': '',
            'unit': 1,
            'de': 0,
            'java_window': 3,
            'java_mult': '',
            'format': 0,
            'line_out': 1,
            'remove_js': 'on',
            # en_unit 0 will return wavenumbers, and 1 will return eV
            'en_unit': 1,
            'output': 0,
            'page_size': 15000,
            'order_out': 0,
            'max_low_enrg': '',
            'show_av': 2,
            'max_upp_enrg': '',
            'tsb_value': 0,
            'min_str': '',
            'A_out': 0,
            'max_str': '',
            'allowed_out': 1,
            'forbid_out': 1,
            'min_accur': '',
            'min_intens': '',
            'conf_out': 'on',
            'term_out': 'on',
            'enrg_out': 'on',
            'J_out': 'on',
            'submit': 'Retrieve Data',
        }

        get_postfix = urllib.parse.urlencode(values)
        with urllib.request.urlopen(url + '?' + get_postfix) as response:
            html_doc = response.read()

        # parsing with lxml is slightly, but noticably faster
        soup = BeautifulSoup(html_doc, 'lxml')
        # the NIST page has 6 HTML tables on it. The tables of transitions can
        # be most easily identified by the bgcolor attribute
        for s in soup.find_all('table'):
            if re.search('#FFFEEE', s.get('style')) is not None:
                table = s

        # grab the first table row to get the table headers
        header_cols = table.find('tr').find_all('th')
        data_header = []
        for col in header_cols:
            try:
                span = int(col['colspan'])
            except:
                span = 1
            for i in range(span):
                data_header.append(col.text.encode('ascii', errors='ignore').
                                   decode('ascii').strip())

        # This parses out some of the spaces and characeters we don't care for
        regex = re.compile('[^a-zA-Z\.\d\/\+\-[]]')
        data = []

        rows = table.findAll('tr', {'class': ['evn', 'odd']})
        for row in rows:
            cols = row.find_all('td')
            data_row = []
            for col in cols:
                data_row.append(regex.sub('', col.text))
            data.append(data_row)

        return pd.DataFrame(data, columns=data_header)

    def parse_nist_transitions(self, data):
        df = data.copy()
        df = df[['Ei(eV)', 'Lower LevelConf.,Term,J', 'Aki(s-1)',
                 'Ek(eV)', 'Upper LevelConf.,Term,J']]
        df.columns = ['E_i', 'conf_i', 'term_i', 'J_i', 'gamma',
                      'E_k', 'conf_k', 'term_k', 'J_k']
        df.dropna(inplace=True)
        df.drop(df[df.E_i == ''].index, inplace=True)
        df.drop(df[df.E_k == ''].index, inplace=True)
        df.drop(df[df.term_i == ''].index, inplace=True)
        df.drop(df[df.term_k == ''].index, inplace=True)
        df.drop(df[df.J_i == ''].index, inplace=True)
        df.drop(df[df.J_k == ''].index, inplace=True)
        df.E_i = df.E_i.map(lambda x: re.sub('[^\de\-\.]+', '', x))
        df.E_k = df.E_k.map(lambda x: re.sub('[^\de\-\.]+', '', x))
        df.E_i = pd.to_numeric(df.E_i)
        df.E_k = pd.to_numeric(df.E_k)
        df.gamma = df.gamma.map(lambda x: re.sub('[^\de\-\+\.]+', '', x))
        df.gamma = pd.to_numeric(df.gamma)

        df['S_i'] = df.term_i.map(lambda x: (float(re.sub('[^\d]+', '', x))-1)/2)
        df['S_k'] = df.term_k.map(lambda x: (float(re.sub('[^\d]+', '', x))-1)/2)
        df['L_i'] = df.term_i.map(lambda x: re.search('([A-Z]|\[\d+\/\d+\])', x).group())
        df['L_k'] = df.term_k.map(lambda x: re.search('([A-Z]|\[\d+\/\d+\])', x).group())
        df.L_i = df.L_i.map(lambda x: re.sub('[\[\]]+', '', x))
        df.L_k = df.L_k.map(lambda x: re.sub('[\[\]]+', '', x))
        df.L_i = df.L_i.map(lambda x: L(x))
        df.L_k = df.L_k.map(lambda x: L(x))
        df.J_i = df.J_i.map(lambda x: float(fractions.Fraction(x)))
        df.J_k = df.J_k.map(lambda x: float(fractions.Fraction(x)))
        df.conf_i = df.conf_i.map(lambda x: re.sub('[^a-z\d]+', '', x))
        df.conf_k = df.conf_k.map(lambda x: re.sub('[^a-z\d]+', '', x))
        df['n_i'] = pd.to_numeric(df.conf_i.map(
                        lambda x: re.search('(\d+)[a-z]\d?$', x).group(1)))
        df['n_k'] = pd.to_numeric(df.conf_k.map(
                        lambda x: re.search('(\d+)[a-z]\d?$', x).group(1)))
        df['wl'] = c_.h*c_.c/(df.E_k - df.E_i)/c_.e

        def trim_digit(x):
            d = 10**np.floor(np.log10(x))
            return int(x-d*np.floor(x/d))
        df.n_i = df.n_i.map(lambda x: x if x < 20 else int(trim_digit(x)))
        df.n_k = df.n_k.map(lambda x: x if x < 20 else int(trim_digit(x)))
        self.transitions = df

    def get_states(self):
        s_i = self.transitions[['E_i', 'n_i', 'S_i', 'L_i', 'J_i']].drop_duplicates()
        s_i.columns = ['E', 'n', 'S', 'L', 'J']
        s_k = self.transitions[['E_k', 'n_k', 'S_k', 'L_k', 'J_k']].drop_duplicates()
        s_k.columns = ['E', 'n', 'S', 'L', 'J']
        self.states = pd.concat([s_i, s_k]).drop_duplicates()
        self.states.sort_values('E', inplace=True)
        self.states.reset_index(drop=True, inplace=True)

    def populate_states(self, local=False):
        """Retreive atomic state data.

        # Keyword arguments
        * If `local` is false the data is collected from the NIST atomic specta database. Otherwise a previously saved csv file is loaded.

        # Raises
        IOError when csv file does not exist.
        """

        if local:
            self.load_states()
        else:
            self.parse_nist_transitions(self.get_nist_transitions(self.name))
            self.get_states()
            self.save_states()

    def save_states(self, filename=''):
        """Save state information as csv file. """

        if filename == '':
            filename = self.data_folder + self.name + '.csv'
        self.transitions.to_csv(filename, encoding='utf8')

    def load_states(self, filename=''):
        """Load state information from a csv file. """

        if filename == '':
            filename = self.data_folder + self.name + '.csv'
        try:
            self.transitions = pd.read_csv(filename, encoding='utf8')
            self.get_states()
        except:
            raise IOError('Could not load atomic transitions file')


    def _flip_operator(self, symbol):
        """Convert '>' to '<', '>=' to '<=', and vice versa. Doesn't do anything to = sign. """

        if symbol[0] == '<':
            symbol = '>' + symbol[1:]
        elif symbol[0] == '>':
            symbol = '<' + symbol[1:]
        return symbol

    def _find_conditions(self, search_string, match_pattern):
        """Find the condition `match_pattern` in a string `search_string`. 

        # Detailed description
        Find all occurances of `match_pattern` in `search_string` and convert it into a relation <number> <operator> <term>.
        Examples
        * 1<n<4 will be converted to [1, 4], ['<', '>'], n which has to be understood as (1 < n) and (4 > n). 
        * 1>=L will be converted to [1], ['>='], L
        * 1<S will be converted to [1], ['<'], S
        """

        matches = []
        symbols = []
        a = None
        b = None
        term = None

        for s in re.finditer(match_pattern, search_string):
            matches.append((s.start(), s.end()))
            symbols.append(search_string[s.start(): s.end()])
        
        if len(matches) == 2: 
            # a < term < b or a > term > b
            # rewrite as (a < term, b > term) or (a > term, b < term)
            term = str(search_string[matches[0][1]:matches[1][0]]) # the thing between the two conditions
            a = float(search_string[:matches[0][0]])
            b = float(search_string[matches[1][1]:]) 
            symbols[1] = self._flip_operator(symbols[1])
        elif len(matches) == 1: # c < term, term < c, c > term or term > x
            ex1 = search_string[:matches[0][0]]
            ex2 = search_string[matches[0][1]:]
            try: # c < term or c > term
                a = float(ex1)
                term = str(ex2)
            except ValueError: # term < c or term > c -> turn operator around
                try:
                    a = float(ex2)
                    term = str(ex1)
                    symbols[0] = self._flip_operator(symbols[0])
                except: # <= or >= symbol was mistaken for == symbol. Reset symbols
                    symbols = []

        return [x for x in [a, b] if x is not None], term, symbols

    def select_states(self, condition):
        """Select state data according to given conditions.

        # Detailed description
        `condition` is a string or list of strings specifying the conditions with which state data is selected.

        `conditions` can have the form 
        <number> <operator> <term> <operator> <number> where the left or right side of <term> can also be omitted.

        <term> must be in {E, n, L, S, J} and <operator> in {<, <=, >, >=, =, ==}. 

        Example conditions are
        * 1 < n < 5
        * 7 > n > 1
        * E >= 2
        * L = 1
        * L == 1
        and you can select data with several conditions with `atom.select_states(['n<=4', 'S=0', '1<=J<3', 'E>3'])`

        # Returns
        * `selected_data` -- The state data in self.states that satisfies `condition`
        """

        if not isinstance(condition, (list, tuple)):
            condition = [condition]

        selected_data = self.states

        for cond in condition:
            # remove white spaces
            cond = ''.join(cond.split(' '))
            term = None
            # search all less and greater (equal) conditions
            n_l, term_l, operators_l = self._find_conditions(cond, '<={0,1}') # matches < and <=
            n_g, term_g, operators_g = self._find_conditions(cond, '>={0,1}') # matches > and >=
            if term_l is not None:
                term = term_l
            elif term_g is not None:
                term = term_g
            n = n_l + n_g
            operators = operators_l + operators_g

            if term is None: # we couldn't find a <, <=, >, >= condition. Hence it must be an == condition.
                # We don't want to call this when there are also <= and >= signs since the matching pattern will also match those.
                n, term, operators = self._find_conditions(cond, '={1,2}') # matches = and == but also <= and >=.     
    
            for (num, op) in zip(n, operators):
                if op == '<':
                    selected_data = selected_data.loc[num < selected_data[term]]
                elif op == '<=':
                    selected_data = selected_data.loc[num <= selected_data[term]]
                elif op == '>':
                    selected_data = selected_data.loc[num > selected_data[term]]
                elif op == '>=':
                    selected_data = selected_data.loc[num >= selected_data[term]]
                elif op in ['=', '==']:
                    selected_data = selected_data.loc[num == selected_data[term]]
        
        return selected_data




def to_SI(d):
    inc_prefix = ['k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    dec_prefix = ['m', 'u', 'n', 'p', 'f', 'a', 'z', 'y']

    degree = int(np.floor(np.log10(np.fabs(d))/3))
    prefix = ''
    if abs(degree)>0:
        if degree>0:
            if degree < len(inc_prefix)+1:
                prefix = inc_prefix[degree-1]
        else:
            if degree < len(dec_prefix)+1:
                prefix = dec_prefix[-degree-1]
        return '%.3g' % (d*10**(-3*degree)) + prefix
    else:
        return '%.3g' % d


L_dict = {'S': 0, 'P': 1, 'D': 2,
          'F': 3, 'G': 4, 'H': 5,
          'I': 6, 'J': 7, 'K': 8}
L_inv = {v: k for k, v in L_dict.items()}


def L(l):
    if l in L_dict:
        return L_dict[l]
    else:
        return float(fractions.Fraction(l))


def L_str(l):
    if l in L_inv:
        return L_inv[l]
    else:
        return str(fractions.Fraction(l))

#a = Atom('mg+')
#a.load_states()
#a.states = a.states[0:8]