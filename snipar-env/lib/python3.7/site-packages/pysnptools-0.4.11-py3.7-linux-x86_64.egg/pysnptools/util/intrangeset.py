from __future__ import absolute_import
import sys
import re
from bisect import bisect_left
import logging
import numpy as np
import unittest
import doctest
import six
from six.moves import range
import numbers
import io

class IntRangeSet(object):
    '''
    A class for efficiently manipulating ranges of integers (including negatives and longs) using set operations such as :meth:`union`, :meth:`intersection`, and difference.

    The class differs from the built-in *set* class (and from Boolean numpy arrays) because it does not need to store every element in the set, only for every contiguous range of elements.
    It differs from other Python interval libraries (that we know of) by being specialized and optimized for integer elements.

    :Example:

    Here we take the union (operator "|") of two IntRangeSets:

    .. image:: example1.png

    >>> from __future__ import print_function #Python 2 & 3 compatibility
    >>> a = IntRangeSet("100:500,501:1000") # a is the set of integers from 100 to 500 (exclusive) and 501 to 1000 (exclusive)
    >>> b = IntRangeSet("-20,400:600")      # b is the set of integers -20 and the range 400 to 600 (exclusive)
    >>> c = a | b                           # c is the union of a and b, namely -20 and 100 to 1000 (exclusive)
    >>> print(c)
    IntRangeSet('-20,100:1000')

    :Example:

    Suppose we want to find the intron regions of a gene but we are given only the transcription region and the exon regions.

    .. image:: example2.png

    >>> from pysnptools.util import IntRangeSet
    >>> line = "chr15   29370   37380   29370,32358,36715   30817,32561,37380"
    >>> chr,trans_start,trans_last,exon_starts,exon_lasts = line.split() # split the line on white space
    >>> trans_start = int(trans_start)
    >>> trans_stop = int(trans_last) + 1 # add one to convert the inclusive "last" value into a Pythonesque exclusive "stop" value
    >>> int_range_set = IntRangeSet((trans_start,trans_stop)) # creates a IntRangeSet from 29370 (inclusive) to 37381 (exclusive)
    >>> print(int_range_set) # print at any time to see the current value
    IntRangeSet('29370:37381')
    
    Parse the exon start and last lists from strings to lists of integers (converting 'last' to 'stop')

    >>> exon_starts = [int(start) for start in exon_starts.strip(",").split(',')]
    >>> exon_stops = [int(last)+1 for last in exon_lasts.strip(",").split(',')]
    >>> assert len(exon_starts) == len(exon_stops)

    Zip together the two lists to create an iterable of exon_start,exon_stop tuples. Then
    'set subtract' all these ranges from int_range_set

    >>> from __future__ import absolute_import #Python 2 & 3 compatibility
    >>> from six.moves import range #Python 2 & 3 compatibility
    >>> int_range_set -= zip(exon_starts,exon_stops)
    >>> print(int_range_set) # See what it looks like
    IntRangeSet('30818:32358,32562:36715')

    Create the desired output by iterating through each contiguous range of integers

    >>> for start, stop in int_range_set.ranges():
    ...    print("{0}\t{1}\t{2}".format(chr, start, stop-1))
    chr15   30818     32357
    chr15   32562     36714


    **Ranges Input**

    The input to the :class:`IntRangeSet` constructor and many of its methods is a *ranges input*.
   
    A *ranges input* is one of the following:

        * A comma-separated string of integers or integer ranges, e.g., ``'100:500,500:1000,2000'``
        * An integer or integer expression, e.g., ``3``
        * A tuple of two integers, *start* and *stop*, e.g., ``(2,8)``. *stop* is exclusive.
        * A slice with non-negative values, e.g. ``slice(2,8)``
        * A :class:`IntRangeSet` (or any class with a :meth:`ranges` method), e.g., ``IntRangeSet(3)``
        * A list or iterable (but not tuple) of *ranges inputs*, e.g., ``[1,6,7,(100,200)]``

        :Example:
        
        Strings:

        >>> a = IntRangeSet("100:500,500:1000,2000")
        >>> a = IntRangeSet('-10:-8,1,2,3:11')
        >>> a = IntRangeSet('')

        The ranges in an *ranges input* can overlap and can be in any order.

        >>> assert IntRangeSet("2000,100:1500,500:1000,2000") == IntRangeSet('100:1500,2000')

        Integers and Integer Expressions:

        >>> a = IntRangeSet(7)
        >>> a = IntRangeSet(151000000000) # longs are OK
        >>> a = IntRangeSet(2*3+1) # integer expressions are OK

        A tuple must have exactly two integers. They represent the *start* (inclusive) and *stop* (exclusive) integers in a range.

        >>> assert IntRangeSet((2,8)) == IntRangeSet('2:8') # check 'set equality'
        >>> assert 7 in IntRangeSet((2,8)) # The integer 7 is an element of the set
        >>> assert IntRangeSet((2,3)) == IntRangeSet(2)

        Lists and iterables of *ranges inputs*:

        >>> assert IntRangeSet([1,6,7,(100,200)]) == IntRangeSet('1,6:8,100:200')
        >>> assert IntRangeSet(range(0,100)) == IntRangeSet('0:100')
        >>> assert IntRangeSet([range(100,200),3,'1000:2000',[4,6],(20,30)]) == IntRangeSet('3:5,6,20:30,100:200,1000:2000')



        Some methods accept zero or more *ranges input* as their input. This is called a *\*ranges_inputs*. For example:

            >>> a = IntRangeSet() # zero ranges inputs
            >>> a = IntRangeSet(3) # one ranges input
            >>> a = IntRangeSet('3:11,5',14,100) # Three ranges inputs, a string and two integers.


    **Most Important Methods and Operators**

    ================================== ================================== ======================
    *Description*                       *Method or Operator*              *For Details*
    is a 'set equal' to b?             ``a == b``                         :meth:`__eq__`
    is a is superset of b?             ``b in a``,                        :meth:`__contains__`
                                       ``a.issuperset(b)``,
                                       ``a >= b``
    is a is subset of b?               ``a <= b``,                        :meth:`__le__`
                                       ``a.issubset(b)``
    remove *i* th smallest element(s)  ``del a[i]``                       :meth:`__delitem__`
    get *i* th smallest element(s)     ``a[i]``                           :meth:`__getitem__`
    is a not 'set equal' to b?         ``a != b``                         :meth:`__ne__`
    is a is proper superset of b?      ``a > b``                          :meth:`__gt__`
    is a is proper subset of b?        ``a < b``                          :meth:`__lt__`
    union b into a                     ``a |= b``,                        :meth:`add`
                                       ``a += b``,
                                       ``a.add(b,c)``,
                                       ``a.update(b,c)``
    union a and b                      ``a | b``,                         :meth:`__or__`
                                       ``a b``,
                                       ``a.union(b,c)``
    intersect b into a                 ``a &= b``,                        :meth:`__iand__`
                                       ``a.intersection_update(b,c)``,
    intersect a and b                  ``a & b``,                         :meth:`intersection`
                                       ``a.intersection(b,c)``
    subtract b from in a               ``a -= b``,                        :meth:`__isub__`
                                       ``a.difference_update(b)``,
                                       ``a.discard(b)``,
    remove b from a, error if missing  ``a.remove(b)``                    :meth:`remove`
    a 'set difference' b                ``a - b``,                        :meth:`__sub__`
                                        ``a.difference(b)``
    iterate integers in a, from low    ``for element in a:``              :meth:`__iter__`
    iterate integers in a, from high   ``for element in reverse(a):``     :meth:`__reversed__`
    symmetric difference into a        ``a ^= b``                         :meth:`__ixor__`
    a 'symmetric difference' b         ``a ^ b``                          :meth:`__xor__`
    count of integer elements in a     ``len(a)``                         :meth:`__len__`
    iterate ranges in a, from low      ``for start,stop in a.ranges():``  :meth:`ranges`
    count of ranges in a               ``a.ranges_len``                   :meth:`ranges_len`
    index of range containing b        ``a.ranges_index(b)``              :meth:`ranges_index`
    get *i* th smallest range          ``a.ranges_getitem(i)``            :meth:`ranges_getitem`
    a as a string                      ``str(a)``                         :meth:`__str__`
    remove all elements from a         ``a.clear()``                      :meth:`clear`
    copy a                             ``a.copy()``                       :meth:`copy`
    find index of b in a               ``a.index(b)``                     :meth:`index`
    are a and b disjoint?              ``a.isdisjoint(b)``                :meth:`isdisjoint`
    is a empty?                        ``a.isempty(b)``                   :meth:`isempty`
    smallest element of a              ``a.min()``                        :meth:`min`
    largest element of a               ``a.max()``                        :meth:`max`
    sum of elements in a               ``a.sum()``                        :meth:`sum`
    remove and return an element       ``a.pop()``                        :meth:`pop`
    ================================== ================================== ======================
    
    **Examples, Tips, and Warnings**

    The argument to a method and the right-hand side of an operator can be a *ranges input* rather than :class:`IntRangeSet`. For example,

        >>> big = IntRangeSet('-1000:2000')
        >>> print(big - '10:20,500') # We can subtract this string because it is a legal ranges input.
        IntRangeSet('-1000:10,20:500,501:2000')

        This even works for equality testing:

        >>> assert big - '10:21,500' == '-1000:10,21:500,501:2000'

        The Python *in* operator is backwards from other operators, so the left-hand side can be any *ranges input*.

        >>> assert 501 in big

    Like most other Python libraries you specify a range with an inclusive *start* integer and an exclusive *stop*.

        >>> print(7 in IntRangeSet('4:8')) # includes 7
        True
        >>> print(8 in IntRangeSet('4:8')) # excludes 8
        False
        >>> print(8 in IntRangeSet(range(4,7))) # also excludes 8
        False

    Be careful with the *ranges inputs* specified with tuples. Suppose we want a :class:`IntRangeSet` containing 4,5,6,7.
        >>> assert IntRangeSet((4,8)) == '4:8' #OK
        >>> assert IntRangeSet([(4,8)]) == '4:8' # List of a tuple is OK
        >>> assert IntRangeSet(4,8) == '4,8' #No. This is not a tuple. It is two integer inputs, so we get only 4 and 8
        >>> assert IntRangeSet([4,8]) == '4,8' #No. List of two integers, so get only 4 and 8
        >>> #Illegal: IntRangeSet((4,8,10)) # tuples must be pairs


    **Methods and Operators**
    '''
    _rangeExpression = re.compile(r"^(?P<start>-?\d+)(\:(?P<stop>-?\d+))?$")

    def __init__(self, *ranges_inputs):
        '''
        Create a :class:`IntRangeSet`.
        '''
        if len(ranges_inputs) > 0 and isinstance(ranges_inputs[0],IntRangeSet): #Because we know self is empty, optimize for the case in which the first item is a IntRangeSet
            self._start_items = list(ranges_inputs[0]._start_items)
            self._start_to_length = dict(ranges_inputs[0]._start_to_length)
            ranges_inputs = ranges_inputs[1:]
        else:
            self._start_items = []
            self._start_to_length = {}
        self.add(*ranges_inputs) 

    def add(self, *ranges_inputs):
        '''
        Union zero or more ranges inputs into the current IntRangeSet.

        These are the same:
        
        * ``a |= b``
        * ``a += b``
        * ``a.add(b)``
        * ``a.update(b)``

        :Example:

        >>> a = IntRangeSet('0:5,6:10')
        >>> a |= 5
        >>> print(a)
        IntRangeSet('0:10')


        The 'add' and 'update' methods also support unioning multiple ranges inputs,

        :Example:

        >>> a = IntRangeSet('0:5,6:10')
        >>> a.add('5','100:200')
        >>> print(a)
        IntRangeSet('0:10,100:200')
        '''

        #!!consider special casing the add of a single int. Anything else?
        for start,stop in IntRangeSet._static_ranges(*ranges_inputs):
            self._internal_add(start, stop-start)
    def __iadd__(self, *ranges_inputs):
        '''See :meth:`IntRangeSet.add`
        '''
        self.add(*ranges_inputs)
        return self
    __ior__ = __iadd__
    update = add


    def __imul__(self, n):
        '''
        Put produces n copies of a unioned together in a, where a is an IntRangeSet.
        Because a is a set, the result will either be an empty IntRangeSet (n is 0 or less) or the original a.

        * ``a *= n``

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> a *= 5
        >>> print(a) # should be unchanged
        IntRangeSet('0:5,6:11')
        >>> a *= 0
        >>> print(a) # should be empty
        IntRangeSet('')
        '''
        if n <= 0:
            self.clear()
        return self


    def copy(self):
        '''
        Create a deep copy of a IntRangeSet.
        '''
        return IntRangeSet(self)

    def ranges(self):
        '''
        Iterate, in order, the ranges of a IntRangeSet as (start,stop) tuples.

        :Example:

        >>> for start,stop in IntRangeSet('0:10,100:200').ranges():
        ...       print("start is {0}, stop is {1}".format(start,stop))
        start is 0, stop is 10
        start is 100, stop is 200

        '''
        for item in self._start_items:
            stop = item + self._start_to_length[item]
            yield item, stop

    def __iter__(self):
        '''
        Iterate, in order from smallest to largest, the integer elements of the IntRangeSet

        :Example:

        >>> for i in IntRangeSet('1:4,10'):
        ...    print(i)
        1
        2
        3
        10

        '''
        for (first, stop) in self.ranges():
            for i in range(first,stop):
                yield i

    def clear(self):
        '''
        Remove all ranges from this IntRangeSet.

        >>> a = IntRangeSet('0:10,12')
        >>> a.clear()
        >>> print(a)
        IntRangeSet('')

        '''
        del self._start_items[:]
        self._start_to_length.clear()

    def __len__(self):
        '''
        The number of integer elements in the IntRangeSet

        >>> print(len(IntRangeSet('0:10,12')))
        11

        Note: This is computed in time linear in the number of ranges, rather than integer elements.

        '''
        return sum(self._start_to_length.values())


    @property
    def ranges_len(self):
        '''
        The number of contiguous ranges in the IntRangeSet

        >>> print(IntRangeSet('0:10,12').ranges_len)
        2

        '''
        return len(self._start_items)

    def ranges_getitem(self, index):
        '''
        Returns the index-th range in the collection of ranges

        >>> print(IntRangeSet('-30:-20,0:10,12').ranges_getitem(1))
        (0, 10)

        '''
        start = self._start_items[index]
        stop = start+self._start_to_length[start]
        return (start,stop)

    def ranges_index(self, element):
        '''
        Returns the ranges index of range containing the element

        >>> int_range_set = IntRangeSet('-30:-20,0:10,12')
        >>> index = int_range_set.ranges_index(5)
        >>> print(index)
        1
        >>> int_range_set.ranges_getitem(index)
        (0, 10)
        '''
        index = bisect_left(self._start_items, element)
        if index == len(self._start_items) or self._start_items[index] != element:
            index -= 1
        start, stop = self.ranges_getitem(index)
        if element < start or stop <= element:
            raise ValueError("element Not Found")
        return index

    def sum(self):
        '''
        The sum of the integer elements in the IntRangeSet

        >>> print(IntRangeSet('0:10,12').sum())
        57

        Note: This is more efficient than ``sum(IntRangeSet('0:10,12'))`` because is computed
        in time linear in the number of ranges, rather than integer elements.
        '''
        result = 0
        for start in self._start_items:
            length = self._start_to_length[start]
            result += (start + start + length - 1)*length//2
        return result

    def __eq__(self, other):#!!  'others' to ranges_input
        '''
        True exactly when the IntRangeSet on the left is *set equivalent* to the ranges input on the right.

        * ``a == b``


        >>> print(IntRangeSet('0:10,12') == IntRangeSet('0:10,12'))
        True
        >>> print(IntRangeSet('0:10,12') == IntRangeSet('0:10'))
        False
        >>> print(IntRangeSet('0:10,12') == IntRangeSet('12,0:5,5:10'))
        True
        >>> print(IntRangeSet('0:10,12') == '0:10,12') # The right-hand can be any ranges input
        True
        '''
        self, other = IntRangeSet._make_args_range_set(self, other)
        if other is None or len(self._start_items)!=len(other._start_items):
            return False
        for i, start in enumerate(self._start_items):
            if start != other._start_items[i] or self._start_to_length[start] != other._start_to_length[start]:
                return False
        return True

    def __ne__(self, other):
        '''
        *False* exactly when the IntRangeSet on the left is *set equivalent* to the ranges input on the right.

        These are the same.

        * ``a != b``
        * ``not a == b``
        '''
        #Don't need to call _make_args_range_set because __eq__ will call it
        return not self==other

    #Same: a >= b, a.issuperset(b,...), b in a
    #Returns True iff item is within the ranges of this IntRangeSet.
    def __contains__(self, *ranges_inputs):
        '''
        True exactly when all the ranges input is a subset of the IntRangeSet.

        These are the same:

        * ``b in a``
        * ``a.issuperset(b)``
        * ``a >= b``


        :Example:

        >>> print(3 in IntRangeSet('0:5,6:11'))
        True
        >>> print(IntRangeSet('4:7') in IntRangeSet('0:5,6:11'))
        False
        >>> '6:9' in IntRangeSet('0:5,6:11') # The left-hand of 'in' can be any ranges input
        True
        >>> print(IntRangeSet('0:5,6:11') >= '6:9') # The right-hand of can be any ranges input
        True

        The 'issuperset' method also supports unioning multiple ranges inputs.

        :Example:

        >>> print(IntRangeSet('0:5,6:11').issuperset(4,7,8))
        True
        >>> print(IntRangeSet('0:5,6:11').issuperset(4,7,8,100))
        False

        Note: By definition, any set is a superset of itself.
        '''
        for start_in,stop_in in IntRangeSet._static_ranges(*ranges_inputs):
            start_self,length_self,index,contains = self._best_start_length_index_contains(start_in)
            if not contains or stop_in > start_self+length_self:
                return False
        return True
    def __ge__(self,other):
        '''See :meth:`IntRangeSet.__contains__`
        '''
        return other in self
    issuperset = __contains__

    
    @property
    def isempty(self):
        '''
        True exactly when the IntRangeSet is empty.

        >>> print(IntRangeSet().isempty)
        True
        >>> print(IntRangeSet(4).isempty)
        False
        '''
        return len(self._start_items) == 0

    def __str__(self):
        '''
        Use the standard str(a) function to create a string representation of a, an IntRangeSet.

        >>> print("Hello " + str(IntRangeSet(2,3,4,10)))
        Hello IntRangeSet('2:5,10')
        '''
        return repr(self)


    def __repr__(self):
        '''
        Use the standard repr(a) function to create a string representation of a, an IntRangeSet.

        >>> print("Hello " + repr(IntRangeSet(2,3,4,10)))
        Hello IntRangeSet('2:5,10')
        '''
        return "IntRangeSet('{0}')".format(self._repr_internal(":", ","))

    def _repr_internal(self, seperator1, separator2):
        if self.isempty:
            return ""

        fp = io.StringIO() if sys.version_info >= (3,0) else io.BytesIO()

        for index, (start, stop) in enumerate(self.ranges()):
            if index > 0:
                fp.write(separator2)

            if start == stop-1:
                fp.write(str(start))
            else:
                fp.write("{0}{1}{2}".format(start, seperator1, stop))
        return fp.getvalue()

    @staticmethod
    def _test():


        assert IntRangeSet('100:110,1000').index('109,100:104') == '0:4,9'
        assert IntRangeSet("-3") == "-3"
        int_range_set = IntRangeSet()
        int_range_set.add(0)
        assert "IntRangeSet('0')" == str(int_range_set)
        int_range_set.add(1)
        assert "IntRangeSet('0:2')" == str(int_range_set)
        int_range_set.add(4)
        assert "IntRangeSet('0:2,4')" == str(int_range_set)
        int_range_set.add(5)
        assert "IntRangeSet('0:2,4:6')" == str(int_range_set)
        int_range_set.add(7)
        assert "IntRangeSet('0:2,4:6,7')" == str(int_range_set)
        int_range_set.add(2)
        assert "IntRangeSet('0:3,4:6,7')" == str(int_range_set)
        int_range_set.add(3)
        assert "IntRangeSet('0:6,7')" == str(int_range_set)
        int_range_set.add(6)
        assert "IntRangeSet('0:8')" == str(int_range_set)
        int_range_set.add(-10)
        assert "IntRangeSet('-10,0:8')" == str(int_range_set)
        int_range_set.add(-5)
        assert "IntRangeSet('-10,-5,0:8')" == str(int_range_set)

        assert IntRangeSet("-10:-4") == "-10:-4"
        assert IntRangeSet("-10:-4,-3") == "-10:-4,-3"
        assert IntRangeSet("-10:-4,-3,-2:2") == "-10:-4,-3:2"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6") == "-10:-4,-3:6"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6,7:13") == "-10:-4,-3:6,7:13"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6,7:13,13:16") == "-10:-4,-3:6,7:16"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6,7:13,13:16,14:17") == "-10:-4,-3:6,7:17"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6,7:13,13:16,14:17,20:26") == "-10:-4,-3:6,7:17,20:26"
        assert IntRangeSet("-10:-4,-3,-2:2,1:6,7:13,13:16,14:17,20:26,22:24") == "-10:-4,-3:6,7:17,20:26"

        s = "-10:-4,-3,-2:2,1:6,7:13,13:16,14:17,20:26,22:24"
        int_range_set = IntRangeSet(s)
        assert int_range_set == "-10:-4,-3:6,7:17,20:26"

        s = "1:6,0,4:11,-10:-4,-12:-2,15:21,12:22,-13"
        int_range_set = IntRangeSet(s)
        assert int_range_set == "-13:-2,0:11,12:22"

        assert len(int_range_set) == 32

        int_range_set1 = IntRangeSet("-10:-4")
        int_range_set2 = int_range_set1.copy()
        assert int_range_set1 is not int_range_set2
        assert int_range_set1 == int_range_set2
        int_range_set2.add(7)
        assert int_range_set1 != int_range_set2

        assert str(IntRangeSet(7)) == "IntRangeSet('7')"
        assert str(IntRangeSet((7,8))) == "IntRangeSet('7')"
        assert str(IntRangeSet((7,11))) == "IntRangeSet('7:11')"
        assert str(IntRangeSet(list(range(7,11)))) == "IntRangeSet('7:11')"
        assert str(IntRangeSet(np.s_[7:11])) == "IntRangeSet('7:11')"
        assert str(IntRangeSet(np.s_[7:11:2])) == "IntRangeSet('7,9')"
        assert str(IntRangeSet(list(range(7,11,2)))) == "IntRangeSet('7,9')"
        assert str(IntRangeSet(None)) == "IntRangeSet('')"
        assert str(IntRangeSet()) == "IntRangeSet('')"
        assert [e for e in IntRangeSet("-10:-4,-3")] == [-10,-9,-8,-7,-6,-5,-3]
        int_range_set3 = IntRangeSet(7,10)
        int_range_set3.clear()
        assert str(int_range_set3) == "IntRangeSet('')" 
        assert len(IntRangeSet("-10:-4,-3")) == 7

        int_range_set4 = IntRangeSet("-10:-4,-3")
        int_range_set4.add(-10,-7)
        assert int_range_set4 == "-10:-4,-3"
        int_range_set4.add(-10,-4)
        assert int_range_set4 == "-10:-2"

        int_range_set5 = IntRangeSet("-10:-4,-3")
        assert -11 not in int_range_set5
        assert -10 in int_range_set5
        assert -5 in int_range_set5
        assert -4 not in int_range_set5
        assert -3 in int_range_set5
        assert -2 not in int_range_set5
        assert 19999 not in int_range_set5
        assert "-11" not in int_range_set5
        assert "-10" in int_range_set5
        assert "-10:-4" in int_range_set5
        assert "-10:-3" not in int_range_set5
        assert "-10:-4,-3" in int_range_set5
        assert "-10:-4,-3,100" not in int_range_set5
        assert IntRangeSet("-11") not in int_range_set5
        assert IntRangeSet("-10") in int_range_set5
        assert IntRangeSet("-10:-4") in int_range_set5
        assert IntRangeSet("-10:-3") not in int_range_set5
        assert IntRangeSet("-10:-4,-3") in int_range_set5
        assert IntRangeSet("-10:-4,-3,100") not in int_range_set5
        assert [-11] not in int_range_set5
        assert [-10] in int_range_set5
        assert list(range(-10,-6)) in int_range_set5
        assert list(range(-10,-3)) not in int_range_set5
        assert [-10,-9,-8,-7,-3] in int_range_set5
        assert [-10,-9,-8,-7,-3,-100] not in int_range_set5

        assert IntRangeSet("-11,-10,-9") == IntRangeSet("-11:-8")
        a = IntRangeSet("-11,-10")
        b = IntRangeSet("-11:-8")
        assert a!=b
        assert IntRangeSet("-11,1") != IntRangeSet("-11:1")

        assert IntRangeSet("1:4") | IntRangeSet("2,4,6") == IntRangeSet("1:5,6")
        assert IntRangeSet("1:4").union("2,4,6","1,6,-100") == IntRangeSet("-100,1:5,6")


        assert IntRangeSet() & IntRangeSet() == IntRangeSet()
        assert IntRangeSet().intersection("","") == IntRangeSet()
        assert IntRangeSet() & IntRangeSet() & IntRangeSet() == IntRangeSet()

        assert IntRangeSet("1") & IntRangeSet() == IntRangeSet()
        assert IntRangeSet() & IntRangeSet("1") == IntRangeSet()
        assert IntRangeSet("1:6") & IntRangeSet() == IntRangeSet()
        assert IntRangeSet() & IntRangeSet("1:6") == IntRangeSet()
        assert IntRangeSet("1:6,7") & IntRangeSet() == IntRangeSet()
        assert IntRangeSet() & IntRangeSet("1:6,7") == IntRangeSet()

        assert IntRangeSet("1") & IntRangeSet("1") == IntRangeSet("1")
        assert IntRangeSet("1:6") & IntRangeSet("1") == IntRangeSet("1")
        assert IntRangeSet("1") & IntRangeSet("1:6") == IntRangeSet("1")
        assert IntRangeSet("1:6,7") & IntRangeSet("1") == IntRangeSet("1")

        assert IntRangeSet("2") & IntRangeSet("1:6,7") == IntRangeSet("2")
        assert IntRangeSet("1:6") & IntRangeSet("2") == IntRangeSet("2")
        assert IntRangeSet("2") & IntRangeSet("1:6") == IntRangeSet("2")
        assert IntRangeSet("1:6,7") & IntRangeSet("2") == IntRangeSet("2")


        assert IntRangeSet("-2") & IntRangeSet("1:6,7") == IntRangeSet()
        assert IntRangeSet("1:6") & IntRangeSet("-2") == IntRangeSet()
        assert IntRangeSet("-2") & IntRangeSet("1:6") == IntRangeSet()
        assert IntRangeSet("1:6,7") & IntRangeSet("-2") == IntRangeSet()

        assert IntRangeSet("22") & IntRangeSet("1:6,7") == IntRangeSet()
        assert IntRangeSet("1:6") & IntRangeSet("22") == IntRangeSet()
        assert IntRangeSet("22") & IntRangeSet("1:6") == IntRangeSet()
        assert IntRangeSet("1:6,7") & IntRangeSet("22") == IntRangeSet()


        assert IntRangeSet("-2,1:4,20,25:100,101") & IntRangeSet("1:101") == IntRangeSet("1:4,20,25:100")
        assert IntRangeSet("2:4,90:110") & IntRangeSet("1:100") == IntRangeSet("2:4,90:100")
        assert IntRangeSet("1:100") & IntRangeSet("-2,1:4,20,25:100,100") == IntRangeSet("1:4,20,25:100")
        assert IntRangeSet("1:100") & IntRangeSet("2:4,90:110") == IntRangeSet("2:4,90:100")

        assert IntRangeSet("0:67") - "30:101" == IntRangeSet("0:30")
        assert IntRangeSet("0:30,51:67").difference("40:101") == IntRangeSet("0:30")
        assert IntRangeSet("0:67").difference("30:50","40:100") == IntRangeSet("0:30")
        assert IntRangeSet("0:67,99,200") - "30:100,300" == "0:30,200"
        assert IntRangeSet("30:100") - "0:67" == "67:100"
        assert IntRangeSet("30:100,300")-IntRangeSet("0:67,100,200") == IntRangeSet("67:100,300")

        assert IntRangeSet("30:100,300")^IntRangeSet("0:67,99,200") == IntRangeSet("0:30,67:99,200,300")
        assert IntRangeSet([1,2]).symmetric_difference([2,3]) == IntRangeSet([1,3])

        assert IntRangeSet("10:15")[0:1] == IntRangeSet("10")
        assert IntRangeSet("10:15,55:60")[0:5] == IntRangeSet("10:15")
        assert IntRangeSet("10:15,55:60")[1:5] == IntRangeSet("11:15")
        assert IntRangeSet("10:15,55:60")[8:9] == IntRangeSet("58")
        assert IntRangeSet("10:15,55:60")[-2:] == IntRangeSet("58:60")
        assert IntRangeSet("10:15,55:60")[0:5:2] == IntRangeSet("10,12,14")

        try:
            IntRangeSet("10:15,55:60")[100]
        except KeyError:
            pass
        try:
            IntRangeSet("10:15,55:60")[-100]
        except KeyError:
            pass

        assert IntRangeSet("10:15,55:60").max() == 59
        assert IntRangeSet("10:15,55:60").min() == 10

        assert IntRangeSet("10:15,55:60") + IntRangeSet("58:61,100") == IntRangeSet("10:15,55:61,100")
        assert IntRangeSet("10:15,55:60") + "58:61,100" == IntRangeSet("10:15,55:61,100")

        mult_test0 = IntRangeSet("10:15,55:60")
        mult_test1 = mult_test0 * 3
        assert mult_test0 == mult_test1
        assert mult_test0 is not mult_test1
        assert (mult_test0 * 0).isempty
        assert (mult_test0 * -1000).isempty

        assert IntRangeSet("10:15,55:60").index(10) == 0
        assert IntRangeSet("10:15,55:60").index(11) == 1
        assert IntRangeSet("10:15,55:60").index(55) == 5
        assert IntRangeSet("10:15,55:60").index(56) == 6
        try:
            IntRangeSet("10:15,55:60").index(9)
        except IndexError:
            pass
        try:
            IntRangeSet("10:15,55:60").index(15)
        except IndexError:
            pass

        assert IntRangeSet("10:15,55:60").index("57,56:57") == '6:8' #returns the index of the start of the contiguous place where 57,"56-56" occurs
        assert IntRangeSet("10:15,55:60").index([10,55]) == '0,5'
        try:
            IntRangeSet("10:15,55:60").index("100:111")
        except IndexError:
            pass
        try:
            IntRangeSet("10:15,55:60").index("14:17")
        except IndexError:
            pass



        try:
            IntRangeSet(3.34)
        except Exception:
            pass
        try:
            3.34 in IntRangeSet()
        except Exception:
            pass

        assert IntRangeSet("10:15,55:60").count(10) == 1
        assert IntRangeSet("10:15,55:60").count(100) == 0

        assert IntRangeSet("10:15,55:60").isdisjoint(100)
        assert not IntRangeSet("10:15,55:60").isdisjoint("57:100")

        assert IntRangeSet("10:15,55:60") <= "10:15,55:60"
        assert IntRangeSet("10:15,55:60").issubset("10:15,55:60")
        assert not IntRangeSet("10:15,55:60") < "10:15,55:60"
        assert IntRangeSet("10:15,55:60") < "9:15,55:60"

        assert IntRangeSet("10:15,55:60") >= "10:15,55:60"
        assert IntRangeSet("10:15,55:60").issuperset("10:15,55:60")
        assert not IntRangeSet("10:15,55:60") > "10:15,55:60"
        assert IntRangeSet("9:15,55:60") > "10:15,55:60"

        update0 = IntRangeSet("9:15,55:60")
        update0.update("10,30,100","30")
        assert update0 == IntRangeSet("9:15,30,55:60,100")

        assert IntRangeSet() != None

        update0 = IntRangeSet()
        update0 |= []
        assert update0 == IntRangeSet()

        update0 = IntRangeSet("9:15,55:60")
        update0.intersection_update("10,100","0:101")
        assert update0 == IntRangeSet("10")

        update0 = IntRangeSet("9:15,55:60")
        update0 |= IntRangeSet("10,30,100") | "30"
        assert update0 == IntRangeSet("9:15,30,55:60,100")

        update0 = IntRangeSet("9:15,55:60")
        update0 += IntRangeSet("10,30,100") + "30"
        assert update0 == IntRangeSet("9:15,30,55:60,100")


        update0 = IntRangeSet("9:15,55:60")
        update0 &= IntRangeSet("10,100") & "0:101"
        assert update0 == IntRangeSet("10")

        update0 = IntRangeSet("9:15,55:60")
        update0 -= IntRangeSet("10,100") - "30"
        assert update0 == IntRangeSet("9,11:15,55:60")

        update0 = IntRangeSet("9:15,55:60")
        update0.difference_update("10,100","30")
        assert update0 == IntRangeSet("9,11:15,55:60")


        update0 = IntRangeSet("9:15,55:60")
        update0.symmetric_difference_update("10,100,30")
        assert update0 == IntRangeSet("9,11:15,30,55:60,100")


        update0 = IntRangeSet("9:15,55:60")
        update0 ^= IntRangeSet("10,100,30")
        assert update0 == IntRangeSet("9,11:15,30,55:60,100")

        remove0 = IntRangeSet("9:15,55:60")
        remove0.remove(9)
        assert remove0 == IntRangeSet("10:15,55:60")

        remove0 = IntRangeSet("9:15,55:60")
        remove0.remove(10,13)
        assert remove0 == IntRangeSet("9,11:13,14,55:60")

        remove0 = IntRangeSet("9:15,55:60")
        try:
            remove0.remove("100")
        except KeyError:
            pass

        remove0 = IntRangeSet("9:15,55:60")
        try:
            remove0.remove(IntRangeSet("100"),101)
        except Exception:
            pass

        discard0 = IntRangeSet("9:15,55:60")
        discard0.discard(9)
        assert discard0 == IntRangeSet("10:15,55:60")

        discard0 = IntRangeSet("9:15,55:60")
        discard0.discard(10,13)
        assert discard0 == IntRangeSet("9,11:13,14,55:60")

        discard0 = IntRangeSet("9:15,55:60")
        discard0.discard("100")
        assert discard0 == "9:15,55:60"

        discard0 = IntRangeSet("9:15,55:60")
        discard0.discard(IntRangeSet("100"),101)
        assert discard0 == "9:15,55:60"

        pop0 = IntRangeSet("9:15,55:60")
        assert pop0.pop() == 59
        assert pop0 == "9:15,55:59"

        pop0 = IntRangeSet([1,3])
        assert pop0.pop() == 3
        assert pop0.pop() == 1
        try:
            pop0.pop()
        except KeyError:
            pass

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[0:5]
        delitem0 == "55:60"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[1:5]
        delitem0 == "10,55:60"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[9:10]
        delitem0 == "10:15,55:59"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[-2:]
        delitem0 == "10:15,55:58"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[0:5:2]
        delitem0 == "11,33,55:60"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[-2]
        delitem0 == "10:15,55:58,59"

        delitem0 = IntRangeSet("10:15,55:60")
        del delitem0[3]
        delitem0 == "10:13,14,55:60"

        delitem0 = IntRangeSet("10:15,55:60")
        try:
            del delitem0[100]
        except KeyError:
            pass
        try:
            del delitem0[-100]
        except KeyError:
            pass

        assert list(reversed(IntRangeSet("10:15,55:60"))) == list(reversed(list(IntRangeSet("10:15,55:60"))))

        IntRangeSet("10:15,55:60").sum() == sum(IntRangeSet("10:15,55:60"))


        add0 = IntRangeSet("1,12:15,55:61,71,102")
        add0.add("12:101")
        assert add0 == "1,12:101,102"

        assert IntRangeSet("1,12:15,55:61,71,102") - "5:72" == "1,102"
        assert IntRangeSet("1,12:15,55:61,71,102") - "12:66" == "1,71,102"
        assert IntRangeSet("1,12:15,55:61,71,102") - "13:57" == "1,12,57:61,71,102"

        a = IntRangeSet('100:200,1000')
        del a['2:11']
        assert a == '100:102,111:200,1000'

        assert IntRangeSet('0:6,6:10') - '3:100' == '0:3'


        a = IntRangeSet('100:200')
        a *= 3
        assert a == IntRangeSet('100:200')
        a *= 0
        assert a.isempty
        a = IntRangeSet('100:200')
        a *= -1
        assert a.isempty

        assert IntRangeSet('').ranges_len == 0
        assert IntRangeSet('3').ranges_len == 1
        assert IntRangeSet('3,6').ranges_len == 2
        assert IntRangeSet('3,4,5,6').ranges_len == 1


        assert IntRangeSet('100:200')['10:20,30'] == IntRangeSet('110:120,130')

        try:
            IntRangeSet('100..200') # not well formed
        except:
            pass

        assert IntRangeSet("1,12:15,55:61,71,102").ranges_getitem(1) == (12,15)
        assert IntRangeSet("1,12:15,55:61,71,102").ranges_getitem(-1) == (102,103)
        try:
            IntRangeSet("1,12:15,55:61,71,102").ranges_getitem(5)
        except:
            pass
        assert IntRangeSet("1,12:15,55:61,71,102").ranges_index(13) == 1
        try:
            assert IntRangeSet("1,12:15,55:61,71,102").ranges_index(13) == 1
        except:
            pass
        try:
            IntRangeSet("1,12:15,55:61,71,102").ranges_index(2)
        except:
            pass
        try:
            IntRangeSet("1,12:15,55:61,71,102").ranges_index(103)
        except:
            pass
        assert IntRangeSet("1,12:15,55:61,71,102").ranges_index(1) == 0
        assert IntRangeSet("1,12:15,55:61,71,102").ranges_index(102) == 4


    #s[i] ith item of s, origin 0 (3) 
    #s[i:j] slice of s from i to j (3)(4) 
    #s[i:j:k] slice of s from i to j with step k (3)(5) 
    def __getitem__(self, key):
        '''
        ``a[i]`` returns the ith integer in sorted order (origin 0) from a, an IntRangeSet

        >>> print(IntRangeSet('100:200,1000')[0])
        100
        >>> print(IntRangeSet('100:200,1000')[10])
        110

        If i is negative, the indexing goes from the end

        >>> print(IntRangeSet('100:200,1000')[-1])
        1000

        Python's standard slice notation may be used and returns IntRangeSets.
        (Remember that the Stop number in slice notation is exclusive.)

        >>> print(IntRangeSet('100:200,1000')[0:10]) # Integers 0 (inclusive) to 10 (exclusive)
        IntRangeSet('100:110')

        >>> print(IntRangeSet('100:200,1000')[0:10:2]) # Integers 0 (inclusive) to 10 (exclusive) with step 2
        IntRangeSet('100,102,104,106,108')

        >>> print(IntRangeSet('100:200,1000')[-3:]) # The last three integers in the IntRangeSet.
        IntRangeSet('198:200,1000')

        An IntRangeSet can also be accessed with any ranges input.

        >>> IntRangeSet('100:200,1000')['0:10,20']
        IntRangeSet('100:110,120')
        '''
        if isinstance(key, numbers.Integral):
            if key >= 0:
                for start in self._start_items:
                    length = self._start_to_length[start]
                    if key < length:
                        return start+key
                    key -= length
                raise KeyError()
            else:
                assert key < 0
                key = -key-1
                for start_index in range(len(self._start_items)):
                    start = self._start_items[-1-start_index]
                    length = self._start_to_length[start]
                    if key < length:
                        return start+length-1-key
                    key -= length
                raise KeyError()
        elif isinstance(key, slice):
            lenx = len(self)
            start_index,stop_index,step_index = key.start,key.stop,key.step
            start_index = start_index or 0
            stop_index = stop_index or lenx
            step_index = step_index or 1

            if step_index == 1:
                return self & (self[start_index],self[stop_index-1]+1)
            else:
                return IntRangeSet(self[index] for index in range(*key.indices(lenx)))
        else:
            start_and_stop_generator = (self._two_index(start_index,stop_index) for start_index,stop_index in IntRangeSet._static_ranges(key))
            return self.intersection(start_and_stop_generator)
            

    #max(s) largest item of s   
    def max(self):
        '''
        The largest integer element in the IntRangeSet

        >>> print(IntRangeSet('0:10,12').max())
        12

        Note: This is more efficient than max(IntRangeSet('0:10,12')) because is computed
        in constant time rather than in time linear to the number of integer elements.
        '''
        start = self._start_items[-1]
        return start + self._start_to_length[start] - 1

    #min(s) smallest item of s   
    def min(self):
        '''
        The smallest integer element in the IntRangeSet

        :Example:

        >>> print(IntRangeSet('0:10,12').min())
        0

        Note: This is more efficient than ``min(IntRangeSet('0:10,12'))`` because is computed
        in constant time rather than in time linear to the number of integer elements.
        '''
        return self._start_items[0]


    def _make_args_range_set(*args):
        for arg in args:
            if arg is None:
                yield None
            elif isinstance(arg,IntRangeSet):
                yield arg
            else:
                yield IntRangeSet(arg)

    def union(*ranges_inputs):
        '''
        Return the union of a IntRangeSet with zero or more ranges inputs. The original IntRangeSet is not changed.

        These are the same:
        
        * ``a | b``
        * ``a + b``
        * ``a.union(b)``

        :Example:

        >>> print(IntRangeSet('0:5,6:10') | 5)
        IntRangeSet('0:10')

        The 'union' method also support unioning multiple ranges inputs,

        :Example:

        >>> print(IntRangeSet('0:5,6:10').union(5,'100:200'))
        IntRangeSet('0:10,100:200')
        '''
        result = IntRangeSet()
        result.add(*ranges_inputs)
        return result
    __or__ = union
    __add__ = union


    #s * n, n shallow copies of s concatenated (2) 
    def __mul__(self, n):
        '''
        ``a * n``, produces n shallow copies of a unioned, where a is an IntRangeSet.
        Because a is a set, the result will either be an empty IntRangeSet (n is 0 or less) or a copy of
        the original IntRangeSet.

        * ``a * n``

        '''
        if n<=0:
            return IntRangeSet()
        else:
            return IntRangeSet(self)

    def index(self, other):
        '''
        If x is an integer, returns the index of x in a, where a is an IntRangeSet.
        If x is an ranges input, returns an IntRangeSet of index of every integer in x.
        Raises an IndexError is x not in a.

        ``* a.index(x)``


        >>> print(IntRangeSet('100:110,1000').index(109))
        9
        >>> print(IntRangeSet('100:110,1000').index('109,100:104'))
        IntRangeSet('0:4,9')
        '''
        if isinstance(other, numbers.Integral):
            return self._index_element(other)
        else:
            #If start and stop are adjacent, only call _index_element once
            start_and_stop_index_generator = ((self._index_element(start),self._index_element(stop-1)+1) for start,stop in IntRangeSet._static_ranges(other))
            return IntRangeSet(start_and_stop_index_generator)

    def _index_element(self, element):
        index = bisect_left(self._start_items, element)

        # Find the position_in_its_range of this element
        if index != len(self._start_items) and self._start_items[index] == element: #element is start value
            position_in_its_range = 0
        elif index == 0:   # item is before any of the ranges
            raise IndexError()
        else:
            index -= 1
            start = self._start_items[index]
            stop = start+self._start_to_length[start]
            if element >= stop: # we already know it's greater than start
                raise IndexError()
            else:
                position_in_its_range = element - start

        # Sum up the length of all preceding ranges
        preceeding_starts = self._start_items[0:index]
        preceeding_lengths = (self._start_to_length[start] for start in preceeding_starts)
        result = position_in_its_range + sum(preceeding_lengths)
        return result

    #s.count(x) total number of occurrences of x in s   
    def count(self, ranges):
        '''
        The number of times that the elements of ranges appears in the IntRangeSet. Because IntRangeSet is 
        a set, the number will be either 0 or 1.

        >>> print(IntRangeSet('100:110,1000').count('105:107,1000'))
        1
        '''
        if ranges in self:
            return 1
        else:
            return 0

    #Return True if the set has no elements in common with other. Sets are disjoint if and only if their intersection is the empty set.
    def isdisjoint(self, ranges):
        '''
        True exactly when the two sets have no integer elements in common.

        :Example:

        >>> print(IntRangeSet('100:110,1000').isdisjoint('900:2000'))
        False
        >>> print(IntRangeSet('100:110,1000').isdisjoint('1900:2000'))
        True
        '''
        isempty_generator = (IntRangeSet(tuple)._binary_intersection(self).isempty for tuple in IntRangeSet._static_ranges(ranges))
        return all(isempty_generator)

    def __le__(self, ranges):
        '''
        True exactly when the IntRangeSet is a subset of the ranges.

        These are the same:
        
        * ``a <= b``
        * ``a.issubset(b)``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') <= '-1:101') # The right-hand can be any ranges input
        True

        Note: By definition, any set is a subset of itself.
        '''
        self, ranges = IntRangeSet._make_args_range_set(self, ranges)
        return self in ranges
    issubset = __le__

    #set < other
    #Test whether the set is a proper subset of other, that is, set <= other and set != other.
    def __lt__(self, ranges):
        '''
        True exactly when the IntRangeSet is a proper subset of the ranges.

        * ``a < b``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') < '-1:101') # The right-hand can be any ranges input
        True

        Note: By definition, no set is a proper subset of itself.
        '''
        self, ranges = IntRangeSet._make_args_range_set(self, ranges)
        return self != ranges and self in ranges

    #set > other
    #Test whether the set is a proper superset of other, that is, set >= other and set != other.
    def __gt__(self, other):
        '''
        True exactly when the IntRangeSet is a proper superset of the ranges.

        * ``a > b``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') > '7:10') # The right-hand can be any ranges input
        True

        Note: By definition, no set is a proper superset of itself.
        '''
        self, other = IntRangeSet._make_args_range_set(self, other)
        return self != other and other in self

    def intersection(*ranges_inputs):
        '''
        Return the intersection of a IntRangeSet and zero or more ranges inputs. The original IntRangeSet is not changed.

        These are the same:
        
        * ``a & b``
        * ``a.intersection(b)``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') & '3:8')
        IntRangeSet('3:5,6:8')

        The 'intersection' method also support intersecting multiple ranges inputs,

        :Example:

        >>> print(IntRangeSet('0:5,6:11').intersection('3:8','4:7'))
        IntRangeSet('4,6')
        '''
        ranges_inputs = IntRangeSet._make_args_range_set(*ranges_inputs) #generator to made every ranges a IntRangeSet
        ranges_inputs = sorted(ranges_inputs,key=lambda int_range_set:len(int_range_set._start_items)) #sort so that IntRangeSet with smaller range_count is first
        result = ranges_inputs[0] #LATER what if no args, empty? The universe?
        for ranges in ranges_inputs[1:]:
            result = result._binary_intersection(ranges)
        return result
    __and__ = intersection

    def _binary_intersection(self,other):
        result = IntRangeSet()

        if self.isempty:
            return result

        index = 0
        start0 = self._start_items[index]
        length0 = self._start_to_length[start0]
        stop0 = start0+length0
        while True:
            start1,length1,index1,contains = other._best_start_length_index_contains(start0)
            stop1=start1+length1
            if contains:
                if stop0 <= stop1: #All of range0 fits inside some range1, so add it the intersection and next look at the next range0
                    result._internal_add(start0,length0)
                    index+=1
                    if index >= len(self._start_items):
                        break # leave the while loop
                    start0 = self._start_items[index]
                    length0 = self._start_to_length[start0]
                    stop0 = start0+length0
                else: #Only part of range0 fits inside some range0, so add that part and then next look at the rest of range0
                    result._internal_add(start0,stop1-start0)
                    start0 = stop1
                    length0=stop0-start0
            else: #start0 is not contained in any range1, so swap self and other and then next look at the next range0
                temp = other
                other = self
                self = temp
                index = index1+1
                if index >= len(self._start_items):
                    break # leave the while loop
                start0 = self._start_items[index]
                length0 = self._start_to_length[start0]
                stop0 = start0+length0
        return result


    #Same a-b, a.difference(b,...)
    #difference(other, ...)set - other - ...
    #Return a new set with elements in the set that are not in the others.
    #Changed in version 2.6: Accepts multiple input iterables.
    def __sub__(self, *ranges_inputs): #!!could be made faster by being more direct instead of using complements
        '''
        Return the set difference of a IntRangeSet with zero or more ranges inputs. The original IntRangeSet is not changed.

        These are the same:

        * ``a - b``
        * ``a.difference(b)``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') - 1)
        IntRangeSet('0,2:5,6:11')
        >>> print(IntRangeSet('0:5,6:11') - '3:100')
        IntRangeSet('0:3')

        The 'difference' method also supports subtracting multiple input ranges

        :Example:

        >>> print(IntRangeSet('0:5,6:11').difference('3:100',1))
        IntRangeSet('0,2')
        '''
        result = self.copy()
        result.difference_update(*ranges_inputs)
        return result
    difference = __sub__

    #same a^b, a.symmetric_difference(b)
    #symmetric_difference(other)set ^ other
    #Return a new set with elements in either the set or other but not both.
    def __xor__(self, ranges):
        '''
        Returns a new IntRangeSet set with elements in either the input IntRangeSet or the input range but not both.

        These are the same:

        * ``a ^ b``
        * ``a.symmetric_difference(b)``

        :Example:

        >>> print(IntRangeSet('0:5,6:11') ^ '3:9')
        IntRangeSet('0:3,5,9:11')
        '''
        result = self - ranges
        diff_generator = (IntRangeSet(tuple)-self for tuple in IntRangeSet._static_ranges(ranges))
        result += diff_generator
        return result
    symmetric_difference = __xor__

    def _clone_state(self, result):
        self._start_items = result._start_items
        self._start_to_length = result._start_to_length
        return self



    def __iand__(*ranges_inputs):
        '''
        Set the IntRangeSet to itself intersected with a input range

        These are the same:

        * ``a &= b``
        * ``a.intersection_update(b)``

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> a &= '3:8'
        >>> print(a)
        IntRangeSet('3:5,6:8')
        '''
        return ranges_inputs[0]._clone_state(IntRangeSet.intersection(*ranges_inputs))
    def intersection_update(*ranges_inputs):
        '''See :meth:`IntRangeSet.__iand__`
        '''
        IntRangeSet.__iand__(*ranges_inputs)

    #same a-=b, a.difference_update(b,...), a.discard(b,...), a.remove(b,...). Note that "remove" is the only one that raises an error if the b,... aren't in a.
    #difference_update(other, ...)set -= other | ...
    #Update the set, removing elements found in others.
    #Changed in version 2.6: Accepts multiple input iterables.
    def __isub__(self, *ranges_inputs):
        '''
        Remove the elements of the range inputs from the IntRangeSet

        These are the same:

        * ``a -= b``
        * ``a.difference_update(b)``
        * ``a.discard(b)``

        ``remove`` is almost the same except that it raises a KeyError if any element of b is not in a.

        * ``a.remove(b)``

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> a -= '3:7'
        >>> print(a)
        IntRangeSet('0:3,7:11')


        The 'difference_update', 'discard' and 'remove' methods also support subtracting multiple ranges inputs.

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> a.difference_update('3:7','8:100')
        >>> print(a)
        IntRangeSet('0:3,7')
        '''
        for start,stop in IntRangeSet._static_ranges(*ranges_inputs):
            self._internal_isub(start, stop-start)
        return self
    def difference_update(self, *ranges_inputs):
        '''See :meth:`IntRangeSet.__isub__`
        '''
        self.__isub__(*ranges_inputs)
    discard = difference_update
    def remove(self, *ranges_inputs):
        '''See :meth:`IntRangeSet__isub__`
        '''
        for start_stop_tuple in IntRangeSet._static_ranges(*ranges_inputs):
            if not start_stop_tuple in self:
                raise KeyError()
            self -= start_stop_tuple

    def __ixor__(self, ranges):
        '''
        Set the IntRangeSet to contains exactly those elements that appear in either itself or the input ranges but not both

        These are the same:

        * ``a ^= b``
        * ``a.symmetric_difference_update(b)``

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> a ^= '3:8'
        >>> print(a)
        IntRangeSet('0:3,5,8:11')
        '''
        return self._clone_state(self ^ ranges)
    def symmetric_difference_update(self, ranges):
        '''See :meth:`IntRangeSet.__ixor__`
        '''
        self.__ixor__(ranges)

    #pop()
    def pop(self):
        '''
        Remove and return the largest integer element from the IntRangeSet. Raises KeyError if the IntRangeSet is empty.

        :Example:

        >>> a = IntRangeSet('0:5,6:11')
        >>> print(a.pop())
        10
        >>> print(a)
        IntRangeSet('0:5,6:10')
        '''
        if self.isempty:
            raise KeyError()
        #Get the last range
        start = self._start_items[-1]
        length = self._start_to_length[start]
        if length == 1:
            del self._start_to_length[start]
            self._start_items.pop()
        else:
            self._start_to_length[start] = length - 1
        return start+length-1


    def __delitem__(self,key):
        '''
        Remove elements from the IntRangeSet by position index. Position index can be specified by an integer with
        negative integers counting from the end. Position indexes can also be specified with slices and a ranges input.

        :Example:

        Removing with an integer position index:

        >>> a = IntRangeSet('100:200,1000')
        >>> del a[2]
        >>> print(a)
        IntRangeSet('100:102,103:200,1000')
        >>> del a[-1]
        >>> print(a)
        IntRangeSet('100:102,103:200')

        :Example:
       
        Removing with a slice:

        >>> a = IntRangeSet('100:200,1000')
        >>> del a[2:11]
        >>> print(a)
        IntRangeSet('100:102,111:200,1000')

        :Example:

        Removing with a ranges input:

        >>> a = IntRangeSet('100:200,1000')
        >>> del a['2:11']
        >>> print(a)
        IntRangeSet('100:102,111:200,1000')
        '''
        if isinstance(key, numbers.Integral):
            if key >= 0:
                for start in self._start_items:
                    length = self._start_to_length[start]
                    if key < length:
                        self -= start+key 
                        return 
                    key -= length
                raise KeyError()
            else:
                assert key < 0
                key = -key-1
                for start_index in range(len(self._start_items)):
                    start = self._start_items[-1-start_index]
                    length = self._start_to_length[start]
                    if key < length:
                        self -= start+length-1-key
                        return 
                    key -= length
                raise KeyError()
        elif isinstance(key, slice):
            lenx = len(self)
            start,stop,step = key.start,key.stop,key.step
            start = start or 0
            stop = stop or lenx
            step = step or 1

            if step == 1:
                self -= (self[start],self[stop-1]+1)
            else:
                self -= (self[index] for index in range(*key.indices(lenx)))
        else:
            start_and_stop_generator = (self._two_index(start_index,stop_index) for start_index,stop_index in IntRangeSet._static_ranges(key))
            self -= (start_and_stop_generator)

    def _two_index(self,start_index,stop_index):
        start = self[start_index]
        if stop_index == start_index+1:
            stop = start+1
        else:
            stop = self[stop_index-1]+1
        return start,stop

    def __reversed__(self):
        '''
        reversed(a) is a generator that produces the integer elements of a in order from largest to smallest.

        :Example:
        
        >>> for i in reversed(IntRangeSet('1:4,10')):
        ...     print(i)
        10
        3
        2
        1
        '''
        for start in reversed(self._start_items):
            length = self._start_to_length[start]
            for item in range(start+length-1, start-1, -1):
                yield item

    @staticmethod
    #This will gather adjacent ranges together into one range, e.g.  1:4,4,5-6 -> 1-6
    def _static_ranges(*iterables):
        iter = IntRangeSet._inner_static_ranges(*iterables)
        try:
            start0,stop0 = next(iter)
            assert start0 < stop0, "Invalid range. Start " + str(start0) + " must be less than stop " + str(stop0) + "."
        except StopIteration:
            return
        while True:
            try:
                start1,stop1 = next(iter)
                assert start1 < stop1, "Invalid range. Start " + str(start1) + " must be less than stop " + str(stop1) + "."
            except StopIteration:
                yield start0,stop0
                return
            if stop0==start1: #We don't try to merge all cases, just the most common
                stop0=stop1
            elif stop1==start0:
                start0=start1
            else:
                yield start0,stop0
                start0,stop0=start1,stop1

    @staticmethod
    def _inner_static_ranges(*iterables):
        for iterable in iterables:
            if isinstance(iterable, numbers.Integral):
                yield iterable,iterable+1
            elif isinstance(iterable,tuple):
                assert len(iterable)==2 and isinstance(iterable[0], numbers.Integral) and isinstance(iterable[1],six.integer_types), "Tuples must contain exactly two int elements that represent the start (inclusive) and stop (exclusive) elements of a range."
                yield iterable[0],iterable[1]
            elif isinstance(iterable,slice):
                start = iterable.start
                stop = iterable.stop
                step = iterable.step or 1
                assert start is not None and start >=0 and stop is not None and start < stop and step > 0, "With slice, start and stop must be nonnegative numbers, stop must be more than start, and step, if given, must be at least 1"
                if step == 1:
                    yield start,stop
                else:
                    for start in range(start,stop,step):
                        yield start, start+1
            elif iterable is None:
                pass
            elif isinstance(iterable,str):
            # Parses strings of the form -10:-4,-2-10,12-12. Spaces are allowed, no other characters are.
            #  will return an empty range
                if iterable == "":
                    pass
                else:
                    for range_string in iterable.split(","):
                        match = IntRangeSet._rangeExpression.match(range_string) #LATER is there a good error message if it is not well formed?
                        if match is None:
                            raise Exception("The string is not well-formed. '{0}'".format(range_string))
                        start = int(match.group("start"))
                        stop = int(match.group("stop") or start + 1)
                        yield start,stop
            elif hasattr(iterable, 'ranges'):
                for start, stop in iterable.ranges():
                    yield start,stop
            elif hasattr(iterable, '__iter__'):
                for start, stop in IntRangeSet._static_ranges(*iterable):
                    yield start,stop
            else:
                raise Exception("Don't know how to construct a IntRangeSet from '{0}'".format(iterable))

    def _internal_add(self, start, length=1): #!! should it be "length" or "stop"
        assert length > 0, "Length must be greater than zero"
        assert len(self._start_items) == len(self._start_to_length)
        index = bisect_left(self._start_items, start)
        if index != len(self._start_items) and self._start_items[index] == start:
            if length <= self._start_to_length[start]:
                return
            else:
                self._start_to_length[start] = length
                index += 1      # index should point to the following range for the remainder of this method
                previous = start
                stop = start + length
        elif index == 0:
            self._start_items.insert(index, start)
            self._start_to_length[start] = length
            previous = start
            stop = start + length
            index += 1  # index_of_miss should point to the following range for the remainder of this method
        else:
            previous = self._start_items[index - 1]
            stop = previous + self._start_to_length[previous]

            if start <= stop:
                new_length = start - previous + length
                assert new_length > 0 # real assert
                if new_length < self._start_to_length[previous]:
                    return
                else:
                    self._start_to_length[previous] = new_length
                    stop = previous + new_length
            else: # after previous range, not contiguous with previous range
                self._start_items.insert(index, start)
                self._start_to_length[start] = length
                previous = start
                stop = start + length
                index += 1

        if index == len(self._start_items):
            return

        # collapse next range into this one
        next = self._start_items[index]
        while stop >= next:
            new_stop = max(stop, next + self._start_to_length[next])
            self._start_to_length[previous] = new_stop - previous #ItemToLength[previous] + ItemToLength[next]
            del self._start_to_length[next]
            del self._start_items[index]
            stop = new_stop
            if index >= len(self._start_items):
                break
            next = self._start_items[index]
        return

    # return the range that has the largest start and for which start<=element
    def _best_start_length_index_contains(self, element):
        index = bisect_left(self._start_items, element)
        if index != len(self._start_items) and self._start_items[index] == element: #element is the start element of some range
            return element, self._start_to_length[element], index, True
        elif index == 0: # element is before any of the ranges
            return element, 0, -1, False
        else:
            index -= 1
            start = self._start_items[index]
            length = self._start_to_length[start]
            return start, length, index, element < start+length # we already know element is greater than start

    def _delete_ranges(self,start_range_index,stop_range_index):
        for range_index in range(start_range_index,stop_range_index):
            del self._start_to_length[self._start_items[range_index]]
        del self._start_items[start_range_index:stop_range_index]

    def _shorten_last_range(self,stop_in,start1,length1,index1):
        delta_start = stop_in-start1
        self._start_items[index1] += delta_start
        del self._start_to_length[start1]
        self._start_to_length[start1+delta_start]=length1-delta_start
        assert len(self._start_items) == len(self._start_to_length)
           
    def _shorten_first_range(self,start_in,start0,length0):
        assert len(self._start_items) == len(self._start_to_length)
        self._start_to_length[start0] = start_in-start0
        assert len(self._start_items) == len(self._start_to_length)

    def _internal_isub(self, start_in, length_in=1): #!! should it be "length" or "stop"?
        assert length_in > 0, "Length must be greater than zero"
        assert len(self._start_items) == len(self._start_to_length)
        stop_in = start_in+length_in

        # return the range that has the largest start and for which start<=element
        start0,length0,index0,contains0 = self._best_start_length_index_contains(start_in)
        if length_in > 1:
            start1,length1,index1,contains1 = self._best_start_length_index_contains(stop_in-1)
        else:
            start1,length1,index1,contains1 = start0,length0,index0,contains0

        #Is the front of first range unchanged, changed, or deleted?
        if not contains0:#unchanged
            #Is the end of last range unchanged, changed, or deleted?
            if not contains1:#unchanged
                self._delete_ranges(index0+1,index1+1) #delete any middle range
            elif start1+length1 == stop_in: # deleted
                self._delete_ranges(index0+1,index1+1) #delete any middle and the last ranges
            else: #changed
                assert index0 < index1, "real assert"
                self._shorten_last_range(stop_in,start1,length1,index1)                #shorten last range
                self._delete_ranges(index0+1,index1)                    #delete any middle ranges
        elif start0 == start_in: # deleted
            #Is the end of last range unchanged, changed, or deleted?
            if not contains1:#unchanged
                self._delete_ranges(index0,index1+1) #delete start range and any middle ranges
            elif start1+length1 == stop_in: # deleted
                self._delete_ranges(index0,index1+1) #delete start range and any middle ranges and last range
            else: #changed
                assert index0 <= index1, "real assert"
                self._shorten_last_range(stop_in,start1,length1,index1)              #shorten last range
                self._delete_ranges(index0,index1)                    #delete start range and any middle ranges
        else: #changed
            #Is the end of last range unchanged, changed, or deleted?
            if not contains1:#unchanged
                self._shorten_first_range(start_in,start0,length0)              #shorten first range
                self._delete_ranges(index0+1,index1+1)                    #delete any middle ranges
            elif start1+length1 == stop_in: # deleted
                self._shorten_first_range(start_in,start0,length0)              #shorten first range
                self._delete_ranges(index0+1,index1+1)                    #delete any middle ranges and last range
            else: #changed
                if index0 == index1: #need to split into two ranges
                    self._start_items.insert(index0+1,stop_in)
                    self._start_to_length[stop_in] = start1+length1-stop_in
                    self._shorten_first_range(start_in,start0,length0)              #shorten first range
                else:
                    self._shorten_last_range(stop_in,start1,length1,index1)              #shorten last range
                    self._shorten_first_range(start_in,start0,length0)              #shorten first range
                    self._delete_ranges(index0+1,index1)                    #delete any middle ranges
        assert len(self._start_items) == len(self._start_to_length)


class TestIntRangeSet(unittest.TestCase):     

    def test_int_range_set(self):
        IntRangeSet._test()

    def test_doc(self):
        import pysnptools.util.intrangeset
        doctest.testmod(pysnptools.util.intrangeset)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    IntRangeSet._test()
    doctest.testmod()

