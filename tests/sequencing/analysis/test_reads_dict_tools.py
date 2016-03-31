
from .reads_dict_tools import ReadIdGen


def test_ReadIdGen():
    rig = ReadIdGen()
    l = [str(rig) for i in range(3)]
    assert len(set(l)) == len(l)
    assert len(set([len(x) for x in l])) == 1


def test_ReadIdGen2():
    rig = ReadIdGen()
    l = [str(rig) for i in range(2)]
    assert l[1] == "M00321:123:000000000-ABCDE:1:1111:22222:00002"
