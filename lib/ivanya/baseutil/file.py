import fileinput


def get_all_ints(fn):
    """Returns a list  of integers from 'fn'.

    Arguments:
        - fn - a string representing input file name. If None or equal to '-' - read
          from STDIN; fn is treated as a single input integer if it is a digit.
    """
    if fn is None:
        fn = '-'

    if fn.isdigit():
        return [int(fn)]   # just a single integer

    all_ints = []
    for line in fileinput.input(files=fn):
        line = line.rstrip()
        if line.isspace():
            continue
        if not line.isdigit():
            raise Exception("Wrong integer '%s'!" % line)
        all_ints.append(int(line))

    return all_ints

