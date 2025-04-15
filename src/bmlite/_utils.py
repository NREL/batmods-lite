import atexit
import warnings

from typing import Callable, Iterable

from tqdm import tqdm


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


class ProgressBar(tqdm):
    """Progress bar."""

    def __init__(self, iterable: Iterable = None, manual: bool = False,
                 desc: str = None, ncols: int = 80, **kwargs) -> None:
        """
        Wraps the progress bar from ``tqdm``, with different defaults. Also
        enables a custom "manual" mode in which the user manually sets the
        progress as a fraction in [0, 1] using ``set_progress``.

        Parameters
        ----------
        iterable : Iterable, optional
            The iterable to use to construct the "automatic" progress bar, by
            default None. 'manual' must be False if 'iterable' is not None.
        manual : bool, optional
            True enables a "manual" mode progress bar, allowing manual updates
            via 'set_progress'. If False (default), 'iterable' cannot be None.
        desc : str, optional
            Prefix description, by default None.
        ncols : int, optional
            Terminal column width, by default 80. The special case of zero will
            display limited stats and time, with no progress bar.

        Raises
        ------
        ValueError
            'iterable' and 'manual' values are conflicting.
        ValueError
            'iterable' cannot be None if 'manual' is False.

        """

        if manual and iterable is not None:
            raise ValueError("'iterable' and 'manual' values are conflicting.")
        elif iterable is None:
            raise ValueError("'iterable' cannot be None if 'manual' is False.")

        kwargs.setdefault('desc', desc)
        kwargs.setdefault('ncols', ncols)
        kwargs.setdefault('ascii', ' 2468█')
        kwargs.setdefault('iterable', iterable)

        self._iter = 0
        self._manual = manual
        if manual:
            kwargs['total'] = 1
            kwargs['bar_format'] = (
                "{l_bar}{bar}|{iter}[{elapsed}<{remaining}, {rate_fmt}]"
            )

        super().__init__(**kwargs)

    def set_progress(self, progress: float) -> None:
        """
        Updates the progress bar percentage and increments the tracked total
        number of iterations for the "manual" mode. Should be called once per
        "iteration", based on the user's definition of an iteration.

        Parameters
        ----------
        progress : float
            Progress fraction in [0, 1].

        Returns
        -------
        None.

        """
        self._iter += 1
        self.n = progress
        self.refresh()

    def format_meter(self, n: int | float, total: int | float, elapsed: float,
                     **kwargs) -> str:
        """
        Wraps the parent ``format_meter`` method to customize stats for the
        "manual" mode. Users should not need to call this method directly.

        Parameters
        ----------
        n : int or float
            Number of finished iterations.
        total : int or float
            The expected total number of iterations. If meaningless (None),
            only basic progress statistics are displayed (no ETA).
        elapsed : float
            Number of seconds passed since start.
        **kwargs : dict, optional
            Extra keyword arguments to pass through to the parent method.

        Returns
        -------
        out : str
            Formatted meter and stats, ready to display.

        """

        if self._manual:
            kwargs['rate'] = self._iter / elapsed if elapsed > 0 else 0

            try:
                perc = n / total
                t = elapsed*(1 - perc) / perc

                m, s = divmod(int(t), 60)
                h, m = divmod(int(m), 60)

                w_hours = f"{h:d}:{m:02d}:{s:02d}"
                wo_hours = f"{m:02d}:{s:02d}"

                kwargs['remaining'] = w_hours if h else wo_hours

            except ZeroDivisionError:
                kwargs['remaining'] = '?'

        kwargs['iter'] = f" {self._iter}it "
        return super().format_meter(n, total, elapsed, **kwargs)

    def reset(self) -> None:
        """
        Resets to 0 iterations for repeated use.

        Returns
        -------
        None.

        """
        self._iter = 0
        super().reset()


def formatwarning(message, category, filename, lineno, line=None):
    """Shortened warning format - used for parameter/pre warnings."""
    return f"\n[bmlite {category.__name__}] {message}\n\n"


def short_warn(message, category=None, stacklevel=1, source=None):
    """Print a warning with the short format from ``formatwarning``."""
    original_format = warnings.formatwarning

    warnings.formatwarning = formatwarning
    warnings.warn(message, category, stacklevel, source)

    warnings.formatwarning = original_format
