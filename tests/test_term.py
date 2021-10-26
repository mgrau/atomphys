from atomphys.term import L, L_inv, parse_term, print_term


def test_L():
    assert L["S"] == 0
    assert L["D"] == 2

    assert L_inv[0] == "S"
    assert L_inv[3] == "F"


def test_parse_term():
    assert parse_term("2S1/2") == {"S": 1 / 2, "L": 0, "J": 1 / 2, "parity": 1}
    assert parse_term("3P*1") == {"S": 1, "L": 1, "J": 1, "parity": -1}
    assert parse_term("4D") == {"S": 3 / 2, "L": 2, "parity": 1}
    assert parse_term("z 4D") == {"S": 3 / 2, "L": 2, "parity": 1}
    assert parse_term("(1,3/2)") == {"J1": 1, "J2": 3 / 2, "parity": 1}
    assert parse_term("2[5/2]5/2") == {"S2": 1 / 2, "K": 5 / 2, "J": 5 / 2, "parity": 1}
    assert parse_term("3[7]") == {"S2": 1, "K": 7, "parity": 1}
    assert parse_term("") == {"parity": 1}
    assert parse_term("  ") == {"parity": 1}
    assert parse_term("*") == {"parity": -1}
    assert parse_term("  * ") == {"parity": -1}


def test_print_term():

    assert print_term(ionization_limit=True) == "Ionization Limit"
    assert print_term(term="3P1", ionization_limit=True) == "Ionization Limit"
    assert print_term(S=1, L=0, J=1, ionization_limit=True) == "Ionization Limit"

    assert print_term(J="") is None
    assert print_term(J="1/2 or 3/2") is None
    assert print_term(term="2D", J="1/2?") is None

    assert print_term(S=0, L=1) == "1P"
    assert print_term(S=1, L=1, J=0) == "3P0"
    assert print_term(S=1, L=1, J=0, parity=-1, include_parity=True) == "3P*0"
    assert print_term(S=1, L=1, J=0, parity=+1, include_parity=True) == "3P0"
    assert print_term(term="3P", J=0) == "3P0"
    assert print_term(term="3P1", J=0) == "3P0"

    assert print_term(J1=2, J2=1 / 2) == "(2,1/2)"
    assert print_term(J1=1, J2=1, J=0) == "(1,1)0"
    assert print_term("(2,2)") == "(2,2)"
    assert print_term("(2,2)", J1=1) == "(1,2)"

    assert print_term(S2=0, K=3, J=0) == "1[3]0"
    assert print_term(S2=1, K=5 / 2) == "3[5/2]"
    assert print_term("3[1]0") == "3[1]0"
    assert print_term("3[1]0", K=2) == "3[2]0"
    assert print_term("3[1]0", S2=3) == "7[1]0"

    assert print_term(L=0) is None
    assert print_term(L=0, J=0) is None
    assert print_term("S") is None
