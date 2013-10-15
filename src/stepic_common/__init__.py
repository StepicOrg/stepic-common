import random
import collections


def rand_N(N, low=0.8):
    return max(1, random.randint(int(low * N), N))


def equal_int(reply, clue):
    return int(reply) == int(clue)


def equal_ints(reply, clue):
    return [int(x) for x in reply.split()] == [int(x) for x in clue.split()]


def equal_string_multiset(reply, clue):
    answer = collections.Counter(reply.split())
    output = collections.Counter(clue.split())
    return answer == output


def equal_int_multiset(reply, clue):
    reply = collections.Counter(int(x) for x in reply.split())
    clue = collections.Counter(int(x) for x in clue.split())
    return reply == clue


def nice(*args):
    """Format args nicely.

    >>> nice("foo")
    'foo'
    >>> nice("foo", 42)
    'foo\\n42'
    >>> nice([1,2,3], [4,5,6])
    '1 2 3\\n4 5 6'
    """
    def flat_nice(x):
        if isinstance(x, list):
            return ' '.join(map(str, x))
        else:
            return str(x)

    return '\n'.join(map(flat_nice, args))
