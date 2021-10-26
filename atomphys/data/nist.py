import csv
import io
import json
import os
import re
import urllib
import urllib.request
from typing import List

current_file = os.path.realpath(__file__)
cache = os.path.join(os.path.dirname(current_file), ".cache")
try:
    os.stat(cache)
except FileNotFoundError:
    os.mkdir(cache)

monovalent = re.compile(r"^[a-z0-9]*\.(?P<n>\d+)[a-z]$")


def remove_annotations(s: str) -> str:
    """remove annotations from energy strings in NIST ASD"""
    # re_energy = re.compile("-?\\d+\\.\\d*|$")
    # return re_energy.findall(s)[0]

    # this is about 3.5× faster than re.findall, but it's less flexible
    # overall this can make a several hundred ms difference when loading
    return s.strip("()[]aluxyz +?").replace("&dagger;", "")


def fetch_states(atom, refresh_cache=False):
    if not refresh_cache:
        try:
            return json.load(open(os.path.join(cache, atom + " states.cache")))
        except FileNotFoundError:
            pass

    url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
    values = {
        "spectrum": atom,
        "units": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "multiplet_ordered": 1,  # energy ordred
        "term_out": "on",  # output the term symbol string
        "conf_out": "on",  # output the configutation string
        "level_out": "on",  # output the energy level
        "unc_out": 0,  # uncertainty on energy
        "j_out": "on",  # output the J level
        "g_out": "on",  # output the g-factor
        "lande_out": "off",  # output experimentally measured g-factor
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + "?" + get_postfix) as response:
        response = response.read()

    data = list(
        csv.DictReader(
            io.StringIO(response.decode()), dialect="excel-tab", restkey="None"
        )
    )

    json.dump(data, open(os.path.join(cache, atom + " states.cache"), "w+"))

    return data


def parse_states(data: List[dict]):
    return [
        {
            **{
                "energy": remove_annotations(state["Level (Ry)"]) + " Ry",
                "term": state["Term"] + state["J"],
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
    ]


def fetch_transitions(atom, refresh_cache=False):
    if not refresh_cache:
        try:
            return json.load(open(os.path.join(cache, atom + " transitions.cache")))
        except FileNotFoundError:
            pass

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
            print(response.headers.get_content_type())
            return []

        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()), dialect="excel-tab")
    data = [transition for transition in data if transition["Aki(s^-1)"]]

    json.dump(data, open(os.path.join(cache, atom + " transitions.cache"), "w+"))

    return data


def parse_transitions(data: List[dict]):
    return [
        {
            "state_i": {
                "energy": transition["Ei(Ry)"] + " Ry",
                "term": transition["term_i"].replace("*", "") + transition["J_i"],
            },
            "state_f": {
                "energy": transition["Ek(Ry)"] + " Ry",
                "term": transition["term_k"].replace("*", "") + transition["J_k"],
            },
            "A": transition["Aki(s^-1)"] + "s^-1",
            "type": transition["Type"],
        }
        for transition in data
    ]
