import re


def abc2num(abc_string):
    """
    abc2num converts numbers represented in ABC (pseudo-base26) to their matching integer values

    Examples:
    abc2num('A') == 1
    abc2num('Z') == 26
    abc2num('AA') == 27
    abc2num('aA') == 27
    """
    assert re.match('[a-z]+$', abc_string, re.IGNORECASE)
    num = 0
    for char in abc_string:
        num = num * 26 + (ord(char.upper()) - ord('A')) + 1
    return num


def num2abc(num):
    """
    num2abc converts integer values to numbers represented in ABC characters (pseudo-base26)

    Examples:
    num2abc(1) == 'A'
    num2abc(26) == 'Z'
    num2abc(27) == 'AA'
    """
    assert num > 0
    num -= 1
    div = num/26
    if div == 0:
        return chr(65+num)
    else:
        return num2abc(div)+chr(65 + num % 26)


def str_2_xy_index(well_string):
    """
    str2well turns an alphanumeric string into x, y coordinates

    Examples:
    str_2_xy_index('A1') == (1, 1)
    str_2_xy_index('H1') == (8, 1)
    str_2_xy_index('A2') == (1, 2)
    str_2_xy_index('H12') == (8, 12)
    """
    # validate well_string proper format (e.g. a-zA-Z row followed by 0-9 column)
    assert re.match('[a-z]+[0-9]+$', well_string, re.IGNORECASE)

    # get row segment
    x_string = re.search('[a-z]+', well_string, re.IGNORECASE).group()

    # get column segment
    y_string = re.search('[0-9]+', well_string, re.IGNORECASE).group()

    #convert row,column strings to integers
    x = abc2num(x_string)
    y = int(y_string)

    return x, y


def xy_index_2_str(x, y):
    """
    xy_index_2_str turns x, y coordinates into alphanumeric string

    Examples:
    xy_index_2_str(1, 1) == 'A1'
    xy_index_2_str(8, 1) == 'H1'
    xy_index_2_str(1, 2) == 'A2'
    xy_index_2_str(8, 12) == 'H12'
    """

    assert x and y > 0  # validate 1-based indexes

    # get row segment
    x_string = num2abc(x)

    # get column segment
    y_string = str(y)

    return x_string+y_string


def xy_index_2_index(x, y, dimentions=(8, 12)):
    """
    xy_index_2_index turns x, y coordinates into a single integer plate index

    Examples:
    xy_index_2_index(1, 1) == 1
    xy_index_2_index(8, 1) == 8
    xy_index_2_index(1, 2) == 9
    xy_index_2_index(8, 12) == 96

    @param x: rows
    @param y: columns
    @param dimentions: plate dimentions, tuple of (rows, columns)
    @return: single integer plate index
    """
    rows, columns = dimentions
    assert x <= rows and y <= columns
    return x + (y-1)*rows  # global plate index is 1-based (e.g. 1-96 instead of 0-95)


def index_2_xy_index(index, dimentions=(8, 12)):
    """

    Examples:
    index_2_xy_index(1) == (1, 1)
    index_2_xy_index(8) == (8, 1)
    index_2_xy_index(9) == (1, 2)
    index_2_xy_index(96) == (8, 12)

    @param index: single integer plate index
    @param dimentions: plate dimentions, tuple of (rows, columns)
    @return:
    """
    rows, columns = dimentions
    assert index <= rows*columns
    index -= 1  # factor for 1 based indexing
    x = index % rows + 1  # factor for 0 based indexing
    y = index / rows + 1  # factor for 0 based indexing
    return x, y


def index2str(well_index_num, dimentions=(8, 12)):
    """

    Examples:
    index2str(1) == 'A1'
    index2str(8) == 'H1'
    index2str(9) == 'A2'
    index2str(96) == 'H12'

    @param well_index_num: single integer plate index
    @param dimentions: plate dimentions, tuple of (rows, columns)
    @return: well index string (e.g. 'A1', 'H12, ..)
    """
    return xy_index_2_str(*index_2_xy_index(well_index_num, dimentions=dimentions))


def str2index(well_index_string, dimentions=(8, 12)):
    """

    Examples:
    str2index('A1') == 1
    str2index('H1') == 8
    str2index('A2') == 9
    str2index('H12') == 96

    @param well_index_string: well index string (e.g. 'A1', 'H12, ..)
    @param dimentions: plate dimentions, tuple of (rows, columns)
    @return: single integer plate index
    """
    return xy_index_2_index(*str_2_xy_index(well_index_string), dimentions=dimentions)


def test():
    #test abc2num
    assert abc2num('A') == 1
    assert abc2num('Z') == 26
    assert abc2num('AA') == 27
    assert abc2num('aA') == 27

    #test num2abc
    assert num2abc(1) == 'A'
    assert num2abc(26) == 'Z'
    assert num2abc(27) == 'AA'

    #test str_2_xy_index
    assert str_2_xy_index('A1') == (1, 1)
    assert str_2_xy_index('H1') == (8, 1)
    assert str_2_xy_index('A2') == (1, 2)
    assert str_2_xy_index('H12') == (8, 12)

    #test xy_index_2_str
    assert xy_index_2_str(1, 1) == 'A1'
    assert xy_index_2_str(8, 1) == 'H1'
    assert xy_index_2_str(1, 2) == 'A2'
    assert xy_index_2_str(8, 12) == 'H12'

    assert xy_index_2_index(1, 1) == 1
    assert xy_index_2_index(8, 1) == 8
    assert xy_index_2_index(1, 2) == 9
    assert xy_index_2_index(8, 12) == 96

    assert index_2_xy_index(1) == (1, 1)
    assert index_2_xy_index(8) == (8, 1)
    assert index_2_xy_index(9) == (1, 2)
    assert index_2_xy_index(96) == (8, 12)

    assert index2str(1) == 'A1'
    assert index2str(8) == 'H1'
    assert index2str(9) == 'A2'
    assert index2str(96) == 'H12'

    assert str2index('A1') == 1
    assert str2index('H1') == 8
    assert str2index('A2') == 9
    assert str2index('H12') == 96
    print 'Test passed successfully'


if __name__ == '__main__':
    test()