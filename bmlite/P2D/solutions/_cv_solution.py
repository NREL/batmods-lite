from ._base_solution import BaseSolution


class CVSolution(BaseSolution):
    """
    Constant voltage solution for P2D simulations.

    Base: :class:`~bmlite.P2D.solutions.BaseSolution`
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
        Faradaic currents against the current at each time step, the
        through-current in each control volume, and the Li-ion and solid-phase
        lithium conservation.

        If the verification returns ``False``, you can see which checks failed
        using ``plotflag``. Any subplots shaded grey failed their test.
        Failures generally suggest that the solver's relative and/or absolute
        tolerance should be adjusted, and/or that the discretization should be
        modified to adjust the mesh. Note that the third row of the figure
        shows conservation of charge in each control volume in each domain.
        This row is not included in the checks and will never shade grey, but
        is useful in debugging.

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
        from ..postutils import _liquid_phase_Li, _solid_phase_Li, voltage

        c = Constants()

        if len(self.postvars) == 0:
            self.post()

        sim, exp = self._sim, self._exp

        an, sep, ca = sim.an, sim.sep, sim.ca

        V_ext = exp['V_ext']
        V_mod = self.y[:, ca.x_ptr['phi_ed'][-1]]
        i_mod = self.postvars['i_ext']
        i_sum = np.tile(i_mod, (an.Nx + sep.Nx + ca.Nx, 1)).T

        i_an = np.sum(  self.postvars['sdot_an'] * an.A_s
                      * (an.xp - an.xm) * c.F, axis=1)
        i_ca = np.sum(  self.postvars['sdot_ca'] * ca.A_s
                      * (ca.xp - ca.xm) * c.F, axis=1)

        sum_ip = self.postvars['sum_ip']

        Li_el_0, Li_el_t = _liquid_phase_Li(self)
        Li_ed_0, Li_ed_t = _solid_phase_Li(self)

        checks = []
        checks.append(np.allclose(V_ext, V_mod, rtol=rtol, atol=atol))
        checks.append(np.allclose(i_mod, -i_an, rtol=rtol, atol=atol))
        checks.append(np.allclose(i_mod, i_ca, rtol=rtol, atol=atol))
        checks.append(np.allclose(i_sum, -sum_ip, rtol=rtol, atol=atol))
        checks.append(np.allclose(1., Li_el_t / Li_el_0, rtol=rtol, atol=atol))
        checks.append(np.allclose(1., Li_ed_t / Li_ed_0, rtol=rtol, atol=atol))

        if plotflag:
            fig, ax = plt.subplots(nrows=3, ncols=3, figsize=[12, 9],
                                   layout='constrained')

            voltage(self, ax[0, 0])

            if V_mod.mean() != 0.:
                ylims = np.array([0.995 * V_mod.mean(), 1.005 * V_mod.mean()])
                ax[0, 0].set_ylim([min(ylims), max(ylims)])

            # Faradaic currents
            ax[0, 1].set_ylabel(r'$i_{\rm ext} + i_{\rm an}$ [A/m$^2$]')
            ax[0, 2].set_ylabel(r'$i_{\rm ext} - i_{\rm ca}$ [A/m$^2$]')

            ax[0, 1].plot(self.t, i_mod + i_an, '-C3')
            ax[0, 2].plot(self.t, i_mod - i_ca, '-C2')

            ymin = min([ax[0, j].get_ylim()[0] for j in range(1, 3)])
            ymax = max([ax[0, j].get_ylim()[1] for j in range(1, 3)])
            ylim = max(np.abs([ymin, ymax, 0.01]))

            for j in range(1, 3):
                ax[0, j].set_ylim([-ylim, ylim])

            # Throughput current and Li conservation
            ax[1, 0].plot(self.t, i_sum + self.postvars['sum_ip'])
            ax[1, 0].set_ylabel(r'$i_{\rm ext} + \sum \ \ i_{x}^{+}$'
                                ' [A/m$^2$]')

            ax[1, 1].plot(self.t, Li_el_t / Li_el_0, '-k')
            ax[1, 1].set_ylabel(r'$C_{\rm Li^+} \ / \ C_{\rm Li^+}^0$ [$-$]')

            ax[1, 2].plot(self.t, Li_ed_t / Li_ed_0, '-k')
            ax[1, 2].set_ylabel(r'$C_{\rm Li,s} \ / \ C_{\rm Li,s}^0$ [$-$]')

            ymin = min([ax[1, j].get_ylim()[0] for j in range(1, 3)])
            ymax = max([ax[1, j].get_ylim()[1] for j in range(1, 3)])
            ylim = max([1 - ymin, ymax - 1])

            for j in range(1, 3):
                ax[1, j].set_ylim([1 - ylim, 1 + ylim])

            # Divergence (conservation of charge)
            ax[2, 0].text(0.7, 0.85, 'Anode', transform=ax[2, 0].transAxes)
            ax[2, 0].plot(self.t, self.postvars['div_i_an'])

            ax[2, 1].text(0.7, 0.85, 'Separator', transform=ax[2, 1].transAxes)
            ax[2, 1].plot(self.t, self.postvars['div_i_sep'])

            ax[2, 2].text(0.7, 0.85, 'Cathode', transform=ax[2, 2].transAxes)
            ax[2, 2].plot(self.t, self.postvars['div_i_ca'])

            ymin = min([ax[2, j].get_ylim()[0] for j in range(3)])
            ymax = max([ax[2, j].get_ylim()[1] for j in range(3)])
            ylim = max(np.abs([ymin, ymax]))

            for j in range(3):
                ax[2, j].set_ylim([-ylim, ylim])
                ax[2, j].set_ylabel(r'$\nabla \cdot (i_{\rm el} + i_{\rm ed})$'
                                    ' [A/m$^3$]')

            for i in range(3):
                for j in range(3):
                    ax[i, j].set_xlabel(r'$t$ [s]')
                    format_ticks(ax[i, j])

            for i in range(2):
                for j in range(3):
                    if not checks[3 * i + j]:
                        ax[i, j].patch.set_facecolor('grey')
                        ax[i, j].patch.set_alpha(0.5)

            fig.get_layout_engine().set(wspace=0.1, hspace=0.1)
            show(fig)

        return self.success and all(checks)
