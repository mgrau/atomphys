import atomphys


def test_states():
    states = atomphys.data.nist.fetch_states("li i", refresh_cache=True)
    assert len(states) == 182
    assert isinstance(states, list)
    assert all([isinstance(state, dict) for state in states])

    parsed_states = atomphys.data.nist.parse_states(states[:10])
    assert len(parsed_states) == 10
    assert isinstance(parsed_states, list)
    assert all([isinstance(state, dict) for state in parsed_states])
    assert all(
        all([attr in state for attr in ["energy", "term", "configuration", "g"]])
        for state in parsed_states
    )


def test_transitions():
    transitions = atomphys.data.nist.fetch_transitions("li i", refresh_cache=True)
    assert len(transitions) == 328
    assert isinstance(transitions, list)
    assert all([isinstance(transition, dict) for transition in transitions])

    parsed_transitions = atomphys.data.nist.parse_transitions(transitions[:50])
    assert isinstance(parsed_transitions, list)
    assert all([isinstance(transition, dict) for transition in parsed_transitions])
    assert all(
        all([attr in transition for attr in ["state_i", "state_f", "A", "type"]])
        for transition in parsed_transitions
    )


def test_no_transitions():
    transitions = atomphys.data.nist.fetch_transitions("rf i", refresh_cache=True)
    assert transitions == []
