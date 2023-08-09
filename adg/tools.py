"""Miscellaneous diagram-unrelated tools for ADG."""


from builtins import object, zip
from itertools import count


def reversed_enumerate(data):
    """Return the index and item of the data in reversed order.

    Args:
        data (iterable data structure): The data to be used..

    Returns:
        tuple: Index and item.

    >>> list(reversed_enumerate(['A', 'B', 'C']))
    [(2, 'C'), (1, 'B'), (0, 'A')]

    """
    for index, item in zip(count(len(data) - 1, -1), reversed(data)):
        yield index, item


class UniqueID(object):
    """Counter making sure of generating a unique ID number for diagrams."""

    __slots__ = 'current'

    def __init__(self):
        """Initialize counter with value 0."""
        self.current = 0
        """int: The unique identifier to be attributedto a diagram."""

    def get(self):
        """Iterate on counter value and return current value.

        Returns:
            int: A unique identifier for the diagram.

        """
        self.current += 1
        return self.current - 1
