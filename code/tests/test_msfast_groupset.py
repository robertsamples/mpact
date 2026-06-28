import pandas as pd

from MSFaST import groupset


class FakeQuery:
    def __init__(self, name='Set1', excl=None, incl=None, colour='#000000'):
        self.name = name
        self.excl = excl or []
        self.incl = incl or []
        self.colour = colour


def make_iondict():
    # Mirrors what run_MSFaST builds: indexed by Compound, with a 'groups'
    # column of space-joined group-name tokens (' ' + group per membership).
    return pd.DataFrame(
        {'groups': [' CtrlA', ' CtrlB Treated', ' Treated', None]},
        index=pd.Index(['c1', 'c2', 'c3', 'c4'], name='Compound'),
    )


def test_no_constraints_keeps_every_ion_with_groups():
    iondict = make_iondict()
    gs = groupset('Set1', FakeQuery(), iondict)
    assert set(gs.ionlist) == {'c1', 'c2', 'c3'}  # c4 excluded: groups is NaN


def test_exclusion_filters_out_matching_ions():
    iondict = make_iondict()
    gs = groupset('Set1', FakeQuery(excl=['CtrlA']), iondict)
    assert 'c1' not in gs.ionlist
    assert set(gs.ionlist) == {'c2', 'c3'}


def test_inclusion_requires_all_listed_groups():
    iondict = make_iondict()
    gs = groupset('Set1', FakeQuery(incl=['CtrlB', 'Treated']), iondict)
    assert gs.ionlist == ['c2']  # only c2 has both CtrlB and Treated


def test_does_not_mutate_the_shared_iondict():
    """Multiple groupsets are constructed from the same iondict reference
    (run_MSFaST builds it once and reuses it per Plot Feature Set) -- one
    groupset's filtering must not affect another's view of the data."""
    iondict = make_iondict()
    groupset('Set1', FakeQuery(excl=['CtrlA']), iondict)
    gs2 = groupset('Set2', FakeQuery(), iondict)
    assert set(gs2.ionlist) == {'c1', 'c2', 'c3'}
