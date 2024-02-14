from ._base_solution import BaseSolution


class CVSolution(BaseSolution):
    """
    Constant voltage solution for SPM simulations.

    Base: :class:`~bmlite.SPM.solutions.BaseSolution`
    """

    __slots__ = []

    def __init__(self, sim: object, exp: dict) -> None:
        super().__init__(sim, exp)

    @property
    def classname(self) -> str:
        """
        Class name. Overwrites ``BaseSolution``.

        Returns
        -------
        classname : str
            Name of current class.
        """
        return 'CVSolution'

    def verify(self, plotflag: bool = False, rtol: float = 5e-3,
               atol: float = 1e-3) -> bool:
        """
        Verifies the solution is mathematically consistent.

        Specifically, for a constant voltage experiment, this method checks
        that the calculated voltage was within tolerance of the boundary
        condition. In addition, this method checks the anodic and cathodic
        Faradaic currents against the current at each time step and the
        solid-phase lithium conservation.

        If the verification returns ``False``, you can see which checks failed
        using ``plotflag``. Any subplots shaded grey failed their test.
        Failures generally suggest that the solver's relative and/or absolute
        tolerance should be adjusted, and/or that the discretization should be
        modified to adjust the mesh.

        Parameters
        ----------
        plotflag : bool, optional
            A flag to see plots showing the verification calculations. The
            default is ``False``.

        rtol : float, optional
            Relative tolerance for comparisons. The default is 5e-3.

        atol : float, optional
            Absolute tolerance for comparisons. The default is 1e-3.

        Returns
        -------
        checks : bool
            ``True`` is all checks are satisfied, ``False`` otherwise.
        """

        import numpy as np
        import matplotlib.pyplot as plt

        from ... import Constants
        from ...plotutils import format_ticks, show
        from ..postutils import _solid_phase_Li, voltage

        c = Constants()

        if len(self.postvars) == 0:
            self.post()

        sim, exp = self._sim, self._exp

        an, ca = sim.an, sim.ca

        V_ext = exp['V_ext']
        V_mod = self.y[:, ca.ptr['phi_ed']]
        i_mod = self.postvars['i_ext']

        Li_ed_0, Li_ed_t = _solid_phase_Li(self)

        j_an_tot = self.postvars['sdot_an'] * an.A_s * an.thick * c.F
        j_ca_tot = self.postvars['sdot_ca'] * ca.A_s * ca.thick * c.F

        checks = []
        checks.append(np.allclose(V_ext, V_mod, rtol=rtol, atol=atol))
        checks.append(np.allclose(1., Li_ed_t / Li_ed_0, rtol=rtol, atol=atol))
        checks.append(np.allclose(i_mod, -j_an_tot, rtol=rtol, atol=atol))
        checks.append(np.allclose(i_mod, j_ca_tot, rtol=rtol, atol=atol))

        if plotflag:
            fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 6],
                                   layout='constrained')

            voltage(self, ax[0, 0])

            if V_mod.mean() != 0.:
                ylims = np.array([0.995 * V_mod.mean(), 1.005 * V_mod.mean()])
                ax[0, 0].set_ylim([min(ylims), max(ylims)])

            # Lithium conservation
            ax[0, 1].set_ylabel(r'$C_{\rm Li,s} \ / \ C_{\rm Li,s}^0$ [$-$]')
            ax[0, 1].plot(self.t, Li_ed_t / Li_ed_0, '-k')

            # Faradaic currents
            ax[1, 0].set_ylabel(r'$i_{\rm ext} + j_{\rm an, tot}$ [A/m$^2$]')
            ax[1, 1].set_ylabel(r'$i_{\rm ext} - j_{\rm ca, tot}$ [A/m$^2$]')

            ax[1, 0].plot(self.t, i_mod + j_an_tot, '-C3')
            ax[1, 1].plot(self.t, i_mod - j_ca_tot, '-C2')

            ymin = min([ax[1, j].get_ylim()[0] for j in range(2)])
            ymax = max([ax[1, j].get_ylim()[1] for j in range(2)])

            for j in range(2):
                ax[1, j].set_ylim([ymin, ymax])

            for i in range(2):
                for j in range(2):
                    ax[i, j].set_xlabel(r'$t$ [s]')
                    format_ticks(ax[i, j])

            for i in range(2):
                for j in range(2):
                    if not checks[2 * i + j]:
                        ax[i, j].patch.set_facecolor('grey')
                        ax[i, j].patch.set_alpha(0.5)

            fig.get_layout_engine().set(wspace=0.1, hspace=0.1)
            show(fig)

        return self.success and all(checks)
