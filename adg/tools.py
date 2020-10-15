"""Miscellaneous diagram-unrelated tools for ADG."""


from builtins import object, zip
from itertools import count

greek_alphabet = ['\\alpha',
                  '\\beta',
                  '\\gamma',
                  '\\delta',
                  '\\epsilon',
                  '\\zeta',
                  '\\eta',
                  '\\theta',
                  '\\iota',
                  '\\kappa',
                  '\\lambda'
                  '\\mu',
                  '\\nu',
                  '\\omicron',
                  '\\pi',
                  '\\rho',
                  '\\sigma',
                  '\\tau',
                  '\\upsilon',
                  '\\phi',
                  '\\chi',
                  '\\psi',
                  '\\omega']


def greek_letter(number):
    """Returns the lower-case Greek letter in LaTeX format.

    Attributes:
        number (int): The index of the letter with start at 0.

    Returns:
        (str): The Greek letter in LaTeX math format.

    """
    return greek_alphabet[number]


def reversed_enumerate(data):
    """Return the index and item of the data in reversed order.

    Args:
        data (iterable data structure): The data to be used..

    Returns:
        (tuple): Index and item.

    >>> list(reversed_enumerate(['A', 'B', 'C']))
    [(2, 'C'), (1, 'B'), (0, 'A')]

    """
    for index, item in zip(count(len(data) - 1, -1), reversed(data)):
        yield index, item


class UniqueID(object):
    """Counter making sure of generating a unique ID number for diagrams.

    Attributes:
        current (int): The unique identifier to be attributedto a diagram.

    """

    __slots__ = 'current'

    def __init__(self):
        """Initialize counter with value 0."""
        self.current = 0

    def get(self):
        """Iterate on counter value and return current value.

        Returns:
            (int): A unique identifier for the diagram.

        """
        self.current += 1
        return self.current - 1
