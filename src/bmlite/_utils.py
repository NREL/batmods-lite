import sys
import atexit
import warnings
from typing import Callable


class ExitHandler:
    """
    Exit handler.

    Use this class to register functions that you want to run just before a
    file exits. This is primarily used to register plt.show() so plots appear
    in both interactive and non-interactive environments, even if the user
    forgets to explicitly call it.

    """
    _registered = []

    @classmethod
    def register_atexit(cls, func: Callable) -> None:
        if func not in cls._registered:
            atexit.register(func)


class ProgressBar:
    """
    Progress bar utility.

    Prints a progress bar to the console based on a known maximum number of
    iterations. For example,

    .. code-block:: python

        import time

        progbar = ProgressBar(10)

        for i in range(10):
            time.sleep(0.5)
            progbar.update(i+1, prefix=f"iteration # {i+1}")

    """

    __slots__ = ('_max_iterations', '_decimals', '_width', '_fill',)

    def __init__(self, max_iterations: int, decimals: int = 1,
                 width: int = 20) -> None:
        """
        Initialize an instance of the class and use the ``update`` method
        inside a loop to display summary updates and overall progress.

        Parameters
        ----------
        max_iterations : int
            Maximum number of iterations that will occur in the loop.
        decimals : int, optional
            Number of decimals to print in percent complete. The default is 1.
        width : int, optional
            Character width of the printed progress bar. The default is 20.

        """

        if max_iterations <= 0:
            raise ValueError("'max_iterations' must be > 0.")
        if decimals < 0:
            raise ValueError("'decimals' must be >= 0.")
        if width <= 0:
            raise ValueError("'width' must be > 0.")

        self._max_iterations = max_iterations
        self._decimals = decimals
        self._width = width
        self._fill = '█'

    def __repr__(self) -> str:  # pragma: no cover

        data = {
            'max_iterations': self._max_iterations,
            'decimals': self._decimals,
            'width': self._width,
        }

        summary = "\n\t".join([f"{k}={v!r}," for k, v in data.items()])

        return f"ProgressBar(\n{summary}\n)"

    def update(self, iteration: int, prefix: str = '',
               suffix: str = '') -> None:
        """
        Print the updated progress bar based on the current iteration. Include
        a prefix and/or suffix to display other important summary info, e.g.,
        errors in a fitting routine.

        Parameters
        ----------
        iteration : int
            Current iteration, in reference to max_iterations.
        prefix : str, optional
            Prefix string. The default is ''.
        suffix : str, optional
            Suffix string. The default is ''.

        Returns
        -------
        None.

        """

        percent = iteration/self._max_iterations

        num_filled = int(percent*self._width)
        num_unfilled = self._width - num_filled

        bar = '|' + self._fill*num_filled + '-'*num_unfilled + '|'
        per_100 = f"{100*percent:.{self._decimals}f}"

        end = '\n' if percent >= 1 else '\r'

        sys.stdout.write("\r" + f"{prefix} {bar} {per_100}% {suffix}{end}")
        sys.stdout.flush()


def formatwarning(message, category, filename, lineno, line=None):
    """Shortened warning format - used for parameter/pre warnings."""
    return f"\n[bmlite {category.__name__}] {message}\n\n"


def short_warn(message, category=None, stacklevel=1, source=None):
    """Print a warning with the short format from ``formatwarning``."""
    original_format = warnings.formatwarning

    warnings.formatwarning = formatwarning
    warnings.warn(message, category, stacklevel, source)

    warnings.formatwarning = original_format
