import math
import os
import logging

def is_prime(n):
    '''
    Tells (somewhat inefficiently) of a number is prime.

    >>> from pysnptools.util.mapreduce1.examples import is_prime
    >>> is_prime(5)
    True
    >>> is_prime(9)
    False
    '''
    assert n == int(n) and n>1, "Expect integers greater than 1"
    for j in range(2,int(math.sqrt(n))+1):
        if n % j == 0:
            return False
    return True

def prime_search0(start,stop):
    '''
    Iterative algorithm for finding prime numbers in a range, but just testing them one-by-one.

    >>> from pysnptools.util.mapreduce1.examples import prime_search0
    >>> prime_search0(2,10)
    [2, 3, 5, 7]
    '''
    assert start < stop, "start must be less than stop"
    prime_list = []
    for i in range(start,stop):
        if is_prime(i):
            prime_list.append(i)
    return prime_list

def prime_search1(start,stop,runner):
    '''
    Distributed algorithm for finding prime numbers in a range, but just testing each number.

    >>> from pysnptools.util.mapreduce1.examples import prime_search1
    >>> from pysnptools.util.mapreduce1.runner import LocalMultiProc
    >>> prime_search1(2,10,LocalMultiProc(4))
    [2, 3, 5, 7]
    '''
    from pysnptools.util.mapreduce1 import map_reduce

    def mapper(i):
        if is_prime(i):
            return i
        else:
            return None

    def reducer(sequence):
        result = []
        for i in sequence:
            if i is not None:
                result.append(i)
        return result

    return map_reduce(range(start,stop),
                        mapper=mapper,
                        reducer=reducer, #lambda sequence: [i for i in sequence if i is not None], #Filter out the None's
                        name="prime_search1",
                        runner=runner)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
