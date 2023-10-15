import csv
import io
import re
import urllib
import urllib.request
from typing import List

from atomphys.term import print_term
from atomphys.util import disk_cache

monovalent = re.compile(r"^[a-z0-9]*\.(?P<n>\d+)[a-z]$")


def remove_annotations(s: str) -> str:
    """remove annotations from energy strings in NIST ASD"""
    # re_energy = re.compile("-?\\d+\\.\\d*|$")
    # return re_energy.findall(s)[0]

    # this is about 3.5× faster than re.findall, but it's less flexible
    # overall this can make a several hundred ms difference when loading
    return s.strip("()[]aluxyz +?").replace("&dagger;", "")


@disk_cache
def fetch_states(atom, refresh_cache=False):
    url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
    values = {
        "spectrum": atom,
        "units": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "multiplet_ordered": 1,  # energy ordred
        "term_out": "on",  # output the term symbol string
        "conf_out": "on",  # output the configutation string
        "level_out": "on",  # output the energy level
        "unc_out": 1,  # uncertainty on energy
        "j_out": "on",  # output the J level
        "g_out": "on",  # output the g-factor
        "lande_out": "on",  # output experimentally measured g-factor
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + "?" + get_postfix) as response:
        response = response.read()

    data = list(
        csv.DictReader(
            io.StringIO(response.decode()), dialect="excel-tab", restkey="None"
        )
    )

    return data


def parse_states(data: List[dict]):
    return [
        {
            **{
                "energy": remove_annotations(state["Level (Ry)"]) + " Ry",
                "term": print_term(state["Term"], J=state["J"]),
                "configuration": state["Configuration"],
                "g": None if state["g"] == "" else float(state["g"]),
            },
            **(
                {"n": int(monovalent.match(state["Configuration"])["n"])}
                if monovalent.match(state["Configuration"])
                else {}
            ),
        }
        for state in data
        if print_term(state["Term"], J=state["J"])
    ]


@disk_cache
def fetch_transitions(atom, refresh_cache=False):
    # the NIST url and GET options.
    url = "http://physics.nist.gov/cgi-bin/ASD/lines1.pl"
    values = {
        "spectra": atom,
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "en_unit": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "line_out": 2,  # only with {1: transition , 2: level classifications}
        "show_av": 5,
        "allowed_out": 1,
        "forbid_out": 1,
        "enrg_out": "on",
        "term_out": "on",
        "J_out": "on",
        "no_spaces": "on",
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + "?" + get_postfix) as response:
        # when there are no transitions ASD returns a texl/html page with the
        # error message "No lines are available in ASD with the parameters selected"
        # rather than the expected text/plain when using format=3
        if response.headers.get_content_type() != "text/plain":
            return []

        response = response.read()

    data = list(csv.DictReader(io.StringIO(response.decode()), dialect="excel-tab"))

    return data


def parse_transitions(data: List[dict]):
    return [
        {
            "state_i": {
                "energy": remove_annotations(transition["Ei(Ry)"]) + " Ry",
                "term": print_term(term=transition["term_i"], J=transition["J_i"]),
            },
            "state_f": {
                "energy": remove_annotations(transition["Ek(Ry)"]) + " Ry",
                "term": print_term(term=transition["term_k"], J=transition["J_k"]),
            },
            "A": transition["Aki(s^-1)"] + "s^-1",
            "type": transition["Type"],
        }
        for transition in data
        if (
            transition["Aki(s^-1)"]
            and print_term(term=transition["term_i"], J=transition["J_i"])
            and print_term(term=transition["term_k"], J=transition["J_k"])
        )
    ]
