"""Miscellaneous diagram-unrelated tools for ADG."""


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
