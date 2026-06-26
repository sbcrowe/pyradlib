# Style guide
The following are notes describing my general style, for my own reference to try and keep things consistent.

PEP-8 guidelines and the NumPy style of documentation are in use. The Ruff formatter is used to assist with this, with the 88 character line maximum used by Black.

Similar classes may be contained with the same file. Ordering within a class is as follows. These groupings are intended to be delimitered with region markers, for ease of collapse in an IDE. Alphabetical order is maintained within each grouping.
 - Class variables
 - Constructors, including both \_\_init\_\_ and relevant \@classmethods.
 - Magic methods, e.g., \_\_str\_\_, \_\_add\_\_.
 - Properties, including those with the \@cached_property decorator.
 - Public instance methods.
 - Public class methods.
 - Public static methods.
 - Private instance methods.
 - Private class methods.
 - Private static methods.
 - Nested or inner classes.