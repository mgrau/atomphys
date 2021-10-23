import enum
import re
from fractions import Fraction


class Coupling(enum.Enum):
    LS = "LS"  # Russell-Saunders coupling
    jj = "jj"
    LK = "LK"  # pair coupling


LS_term = re.compile(r"^(?P<S>\d+)(?P<L>[A-Z])\*?(?P<J>\d+(/\d)?\d*)?")
J1J2_term = re.compile(r"^\((?P<J1>\d+/?\d*),(?P<J2>\d+/?\d*)\)\*?(?P<J>\d+(/\d)?\d*)?")
LK_term = re.compile(r"^(?P<S2>\d+)\[(?P<K>\d+/?\d*)\]\*?(?P<J>\d+(/\d)?\d*)?")

L = {
    "S": 0,
    "P": 1,
    "D": 2,
    "F": 3,
    "G": 4,
    "H": 5,
    "I": 6,
    "K": 7,
    "L": 8,
    "M": 9,
    "N": 10,
    "O": 11,
    "Q": 12,
    "R": 13,
    "T": 14,
    "U": 15,
    "V": 16,
    "W": 17,
    "X": 18,
    "Y": 19,
}
L_inv = {value: key for key, value in L.items()}


def parse_term(term: str) -> dict:
    """
    parse term symbol string in NIST ASD
    """
    if "Limit" in term:
        return {"ionization_limit": True}

    parity = -1 if "*" in term else 1

    match = LS_term.match(term)
    if match is None:
        match = J1J2_term.match(term)
    if match is None:
        match = LK_term.match(term)
    if match is None:
        return {"parity": parity}

    def convert(key, value):
        if key in ["S", "S2"]:
            return float(Fraction((int(value) - 1) / 2))
        if key in ["J1", "J2", "J", "K"]:
            return float(Fraction(value))
        if key == "L":
            return L[value]

    term = {
        key: convert(key, value) for key, value in match.groupdict().items() if value
    }

    return {**term, "parity": parity}


def print_term(include_parity=False, **quantum_numbers) -> str:
    if "ionization_limit" in quantum_numbers:
        return "Ionization Limit"

    P = ""
    if include_parity:
        if "parity" in quantum_numbers and quantum_numbers["parity"] == -1:
            P = "*"

    J = ""
    if "J" in quantum_numbers:
        J = f"{Fraction(quantum_numbers['J'])}"

    # Russell-Saunders coupling
    if all(key in quantum_numbers.keys() for key in ["S", "L"]):
        return f"{2*quantum_numbers['S']+ 1:g}{L_inv[quantum_numbers['L']]}{P}{J}"

    # J1J2 coupling
    if all(key in quantum_numbers.keys() for key in ["J1", "J2"]):
        return f"({Fraction(quantum_numbers['J1'])},{Fraction(quantum_numbers['J2'])}){P}{J}"

    # LK coupling
    if all(key in quantum_numbers.keys() for key in ["S2", "K"]):
        return (
            f"{2*quantum_numbers['S2'] + 1:g}[{Fraction(quantum_numbers['K'])}]{P}{J}"
        )

    return ""
