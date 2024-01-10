from ._base_solution import BaseSolution


class CCSolution(BaseSolution):
    """
    Constant current solution for P2D simulations.

    Base: :class:`~bmlite.P2D.solutions.BaseSolution`
    """

    __slots__ = ["postvars"]

    def __init__(self, sim: object, exp: dict) -> None:
        super().__init__(sim, exp)

        self.postvars = {}

    @property
    def classname(self) -> str:
        """
        Class name. Overwrites ``classname()`` from ``BaseSolution``.

        Returns
        -------
        classname : str
            Name of current class.
        """
        return "CCSolution"

    def post(self) -> None:
        from ..postutils import post

        self.postvars = post(self)

    def plot(self, *args) -> None:
        if "debug" in args or "all" in args:
            from ..postutils import debug

            debug(self)

        if "verify" in args or "all" in args:
            from ..postutils import verify

            verify(self)

        if "general" in args or "all" in args:
            from ..postutils import general

            general(self)

        if "contours" in args or "all" in args:
            from ..postutils import contours

            contours(self)
