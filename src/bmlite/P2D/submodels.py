from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from numpy import ndarray

    from .domains import Electrode


class Hysteresis:
    """Hysteresis submodel."""

    def __init__(self, domain: Electrode, **options) -> None:
        """
        Sets up a hysteresis submodel to run within the given 'Electrode'. The
        hysteresis is solved in time according to:

        dh/dt = abs(sdot*g_hyst / 3600. / capacity)*(sign(sdot) - h)

        Here, h is unitless and is limited to the range [-1, +1]. When applied
        to modify the OCV, h is scaled to H = M_hyst*h. The equivalent OCV is
        shifted by +H after this scaling. Note that M_hyst must be defined in
        the material using 'get_Mhyst(xs)`.

        Parameters
        ----------
        domain : Electrode
            An instance of an anode or cathode class for which the hysteresis
            submodel will be added.
        **options : dict
            Required keyword options to set the parameters which are needed
            only when hysteresis is enabled. A full list is given below. Note
            there are no defaults, so all values must be provided.
        g_hyst : float
            The hysteresis transition rate. Larger values make h evolve between
            -1 and +1 more quickly.
        hyst0 : float
            The initial, unitless hysteresis. Should be in [-1, +1].

        Raises
        ------
        TypeError
            'Hysteresis' submodels can only be added to 'Electrode' class
            instances.
        ValueError
            The electrode's 'material' class requires a 'get_Mhyst' method to
            enable 'Hysteresis'.

        """

        from .domains import Electrode

        if not isinstance(domain, Electrode):
            raise TypeError("'Hysteresis' submodels can only be added to"
                            " 'Electrode' class instances.")

        elif not hasattr(domain._material, 'get_Mhyst'):
            raise ValueError("The electrode's 'material' class requires a"
                             " 'get_Mhyst' method to enable 'Hysteresis'.")

        self.domain = domain
        domain.g_hyst = options.pop('g_hyst')
        domain.hyst0 = options.pop('hyst0')

    def make_mesh(self, last_xvar: str, pshift: int = 0) -> str:
        """
        Set up mesh and pointer properties specific to the hysteresis submodel.

        Parameters
        ----------
        last_xvar : str
            The last x-only variable key that was added to the domain's pointer
            dictionary.
        pshift : int, optional
            A shifted index to use for the pointer to account for the fact that
            domains are stacked into a single 1D array, by default 0.

        Returns
        -------
        name : str
            The new name (key) of the variable added to the domain pointer.

        """

        domain = self.domain
        domain.ptr['hyst'] = domain.ptr[last_xvar] + 1
        return 'hyst'

    def sv0(self, sv0: ndarray) -> None:
        """
        Adds the initial hysteresis voltage state to sv0.

        Parameters
        ----------
        sv0 : 1D array
            An array instance with pre-allocated memory for the initial state
            vector.

        Returns
        -------
        None.

        """

        domain = self.domain
        start = domain.ptr['start']
        sv0[domain.x_ptr['hyst'] - start] = domain.hyst0

    def algidx(self, algdix: ndarray) -> None:
        """
        Adds algebraic indices to the algidx array, if necessary. Since the
        hysteresis submodel is purely algebraic, this function does nothing.

        Parameters
        ----------
        algdix : 1D array
            Indices of the system's algebraic governing equations.

        """
        pass

    def to_dict(self, soln: object) -> None:
        """
        Returns a dictionary of the state/output variables that are specific
        to the hysteresis submodel.

        Parameters
        ----------
        soln : StepSolution | CycleSolution
            A solution instance from an IVPExperiment.

        Returns
        -------
        soln : dict
            A sliced and organized dictionary with all of the hysteresis state
            and output variables.

        """

        domain = self.domain

        xs_ptr = domain.xr_ptr['xs']
        xs_R = soln.y[:, xs_ptr[:, -1]]

        hyst = soln.y[:, domain.x_ptr['hyst']]
        Hyst = domain.get_Mhyst(xs_R)*hyst

        return {'hyst': hyst, 'Hyst': Hyst}
