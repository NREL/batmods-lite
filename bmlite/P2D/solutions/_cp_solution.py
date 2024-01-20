from ._base_solution import BaseSolution


class CPSolution(BaseSolution):
    """
    Constant power solution for P2D simuations.

    Base: :class:`~bmlite.P2D.solutions.BaseSolution`
    """

    __slots__ = ['postvars']

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
        return 'CPSolution'

    def post(self) -> None:
        from ..postutils import post

        self._sim._flags['BC'] = 'power'
        self.postvars = post(self)
        self._sim._flags['BC'] = None
