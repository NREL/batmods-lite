from ._base_solution import BaseSolution


class CVSolution(BaseSolution):
    """
    Constant voltage solution for SPM simuations.

    Base: :class:`~bmlite.SPM.solutions.BaseSolution`
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
        return 'CVSolution'

    def post(self) -> None:
        from ..postutils import post

        self._sim._flags['BC'] = 'voltage'
        self.postvars = post(self)
        self._sim._flags['BC'] = None

    def verify(self, plotflag: bool = False) -> bool:
        import numpy as np
        import matplotlib.pyplot as plt

        from ... import Constants
        from ...plotutils import format_ticks, show
        from ..postutils import voltage

        c = Constants()

        if len(self.postvars) == 0:
            self.post()

        sim, exp = self._sim, self._exp

        an, ca = sim.an, sim.ca

        V_ext = exp['V_ext']
        V_mod = self.y[:, ca.ptr['phi_ed']]
        i_mod = self.postvars['i_ext']

        i_an = -self.postvars['sdot_an'] * an.A_s * an.thick * c.F
        i_ca = self.postvars['sdot_ca'] * ca.A_s * ca.thick * c.F

        checks = []
        checks.append(np.min(V_mod / V_ext) >= 0.995)
        checks.append(np.max(V_mod / V_ext) <= 1.005)
        checks.append(np.min(i_an / i_mod) >= 0.995)
        checks.append(np.max(i_an / i_mod) <= 1.005)
        checks.append(np.min(i_ca / i_mod) >= 0.995)
        checks.append(np.max(i_ca / i_mod) <= 1.005)

        if plotflag:
            fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[15, 3.5])

            voltage(self, ax[0])

            ax[0].set_ylim([0.995 * V_ext, 1.005 * V_ext])

            ax[1].set_ylabel(r'$-i_{\rm an} / i_{\rm ext}$ [$-$]')
            ax[2].set_ylabel(r'$i_{\rm ca} / i_{\rm ext}$ [$-$]')

            ax[1].plot(self.t, i_an / i_mod, '-C3')
            ax[2].plot(self.t, i_ca / i_mod, '-C2')

            for i in range(1, 3):
                ax[i].set_ylim([0.995, 1.005])
                ax[i].set_xlabel(r'$t$ [s]')
                format_ticks(ax[i])

            fig.subplots_adjust(wspace=0.3)
            show(fig)

        return all(checks)
