import pytest

from tests.flat_dict import FlatDict

BIG_D = {
    1: {
        2: {
            3: [
                {'+':'a1','-':'a2'},
                {'+':'b1'},
                {'-':'c2'},
            ],
            4: {
                9: [
                    {'-':'d2'},
                ],
            },
        },
        5: [
            {'+':'e1','-':'e2'},
        ],
    },
    6: [
        {'+':'f1'},
    ],
    7: [
    ],
    8: {
    }
}

def test_flat_dict():
    fd = FlatDict(BIG_D)
    assert fd._fd._d == {
        (1,): {'+': ['a1', 'b1', 'e1'], '-': ['a2', 'c2', 'd2', 'e2']},
        (1, 2): {'+': ['a1', 'b1'], '-': ['a2', 'c2', 'd2']},
        (1, 2, 3): {'+': ['a1', 'b1'], '-': ['a2', 'c2']},
        (1, 2, 4): {'-': ['d2']},
        (1, 2, 4, 9): {'-': ['d2']},
        (1, 5): {'+': ['e1'], '-': ['e2']},
        (6,): {'+': ['f1']},
        (7,): {},
    }  # FIXME
    assert fd._fd._keys_tree == {
        (1,): {2, 5},
        (1, 2): {3, 4},
        (1, 2, 3): set(),
        (1, 2, 4): {9},
        (1, 2, 4, 9): set(),
        (1, 5): set(),
        (6,): set(),
        (7,): set(),
        (): {1, 6, 7}
    }  # FIXME
    assert fd[1,2,3] == {'+': ['a1', 'b1'], '-': ['a2', 'c2']}
    assert set(fd.keys((1, 2))) == {3, 4}
    assert set(fd.keys()) == {1, 6, 7}
    assert dict(fd.reads((1, 2))) == {
        3: {'+': ['a1', 'b1'], '-': ['a2', 'c2']},
        4: {'-': ['d2']},
    }
    assert dict(fd.reads()) == {
        1: {'+': ['a1', 'b1', 'e1'], '-': ['a2', 'c2', 'd2', 'e2']},
        6: {'+': ['f1']},
        7: {},
    }
    with pytest.raises(KeyError):
        fd[1,2,5]


def test_flat_dict2():
    fd = FlatDict(BIG_D, deinterlaced_keys=['+','-'])
    assert fd._fd._d == {
        (1,): {'+': ['a1', 'b1', 'e1'], '-': ['a2', 'c2', 'd2', 'e2']},
        (1, 2): {'+': ['a1', 'b1'], '-': ['a2', 'c2', 'd2']},
        (1, 2, 3): {'+': ['a1', 'b1'], '-': ['a2', 'c2']},
        (1, 2, 4): {'+': [], '-': ['d2']},
        (1, 2, 4, 9): {'+': [], '-': ['d2']},
        (1, 5): {'+': ['e1'], '-': ['e2']},
        (6,): {'+': ['f1'], '-': []},
        (7,): {'+': [], '-': []},
    }  # FIXME
    assert fd._fd._keys_tree == {
        (1,): {2, 5},
        (1, 2): {3, 4},
        (1, 2, 3): set(),
        (1, 2, 4): {9},
        (1, 2, 4, 9): set(),
        (1, 5): set(),
        (6,): set(),
        (7,): set(),
        (): {1, 6, 7}
    }  # FIXME
    assert fd[1,2,3] == {'+': ['a1', 'b1'], '-': ['a2', 'c2']}
    assert set(fd.keys((1, 2))) == {3, 4}
    assert set(fd.keys()) == {1, 6, 7}
    assert dict(fd.reads((1, 2))) == {
        3: {'+': ['a1', 'b1'], '-': ['a2', 'c2']},
        4: {'+': [], '-': ['d2']},
    }
    assert dict(fd.reads()) == {
        1: {'+': ['a1', 'b1', 'e1'], '-': ['a2', 'c2', 'd2', 'e2']},
        6: {'+': ['f1'], '-': []},
        7: {'+': [], '-': []},
    }
    assert fd[1,2,5] == {'+': [], '-': []}


def test_sub_flat_dict():
    fd = FlatDict(BIG_D)
    d1 = {k: s for k, s in fd.items(1)}
    d2 = {k: s for k, s in fd.items((1,))}
    for d in (d1,d2):
        assert set(d.iterkeys()) == {2, 5}
        sfd2 = d[2]
        sfd5 = d[5]
        assert set(sfd2.keys()) == {3, 4}
        assert set(sfd2.keys(4)) == {9}
        assert set(sfd5.keys()) == set()
        assert sfd2[()] == {'+': ['a1', 'b1'], '-': ['a2', 'c2', 'd2']}
        assert sfd5[()] == {'+': ['e1'], '-': ['e2']}
        assert dict(sfd2.reads()) == {
            3: {'+': ['a1', 'b1'], '-': ['a2', 'c2']},
            4: {'-': ['d2']},
        }
