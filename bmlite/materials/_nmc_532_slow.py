from numpy import ndarray as _ndarray


class NMC532Slow(object):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally slow NMC532 kinetic and transport properties.

        Differs from ``NMC532Fast`` because the equilibrium potential is
        piecewise here, making it more accurate, but slower to evaluate.

        Parameters
        ----------
        alpha_a : float
            Anodic symmetry factor in Butler-Volmer expression [-].

        alpha_c : float
            Cathodic symmetry factor in Butler-Volmer expression [-].

        Li_max : float
            Maximum lithium concentration in solid phase [kmol/m^3].
        """

        self.alpha_a = alpha_a
        self.alpha_c = alpha_c
        self.Li_max = Li_max

    def get_Ds(self, x: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the lithium diffusivity in the solid phase given the local
        intercalation fraction ``x`` and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction [-].

        T : float
            Battery temperature [K].

        Returns
        -------
        Ds : float | 1D array
            Lithium diffusivity in the solid phase [m^2/s].
        """

        import numpy as np

        from .. import Constants

        c = Constants()

        A = np.array([
        -2.509010843479270e+2,  2.391026725259970e+3, -4.868420267611360e+3,
        -8.331104102921070e+1,  1.057636028329000e+4, -1.268324548348120e+4,
         5.016272167775530e+3,  9.824896659649480e+2, -1.502439339070900e+3,
         4.723709304247700e+2, -6.526092046397090e+1])

        Ds = np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * 2.25 * 10.0**(  A[0] * x**10 + A[1] * x**9 + A[2] * x**8
                           + A[3] * x**7 + A[4] * x**6 + A[5] * x**5
                           + A[6] * x**4 + A[7] * x**3 + A[8] * x**2
                           + A[9] * x + A[10]  )

        return Ds

    def get_i0(self, x: float | _ndarray, C_Li: float | _ndarray,
               T: float) -> float | _ndarray:
        """
        Calculate the exchange current density given the intercalation
        fraction ``x`` at the particle surface, the local lithium ion
        concentration ``C_Li``, and temperature ``T``. The input types for
        ``x`` and ``C_Li`` should both be the same (i.e., both float or both
        1D arrays).

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        C_Li : float | 1D array
            Lithium ion concentration in the local electrolyte [kmol/m^3].

        T : float
            Battery temperature [K].

        Returns
        -------
        i0 : float | 1D array
            Exchange current density [A/m^2].
        """

        import numpy as np

        from .. import Constants

        c = Constants()

        A = np.array([
         1.650452829641290e+1, -7.523567141488800e+1,  1.240524690073040e+2,
        -9.416571081287610e+1,  3.249768821737960e+1, -3.585290065824760e+0])

        i0 = 9 * (C_Li / 1.2)**self.alpha_a \
           * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * (  A[0] * x**5 + A[1] * x**4 + A[2] * x**3 + A[3] * x**2
              + A[4] * x**1 + A[5]  )

        return i0

    def get_Eeq(self, x: float | _ndarray, T: float) -> float | _ndarray:
        """
        Calculate the equilibrium potential given the intercalation
        fraction ``x`` at the particle surface and temperature ``T``.

        Parameters
        ----------
        x : float | 1D array
            Lithium intercalation fraction at ``r = R_s`` [-].

        T : float
            Battery temperature [K].

        Returns
        -------
        Eeq : float | 1D array
            Equilibrium potential [V].
        """

        import numpy as np

        conditions = []
        eval_funcs = []

        conditions.append((x >= 0.34) & (x < 0.342985075))
        eval_funcs.append(lambda x: np.polyval(
                          [ -347.97386056124896, -5.5252695033577215,
                            -2.229260693831332, 4.290455 ], x - 0.34))

        conditions.append((x >= 0.342985075) & (x < 0.345970149))
        eval_funcs.append(lambda x: np.polyval(
                          [ -347.97386056135764, -8.641453718802165,
                            -2.27154942515372, 4.283742 ], x - 0.342985075))

        conditions.append((x >= 0.345970149) & (x < 0.348955224))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7529.727270939448, -11.757636890326081,
                            -2.3324422201546735, 4.276875 ], x - 0.345970149))

        conditions.append((x >= 0.348955224) & (x < 0.351940299))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6800.255111880196, 55.672765009572274,
                            -2.201352269084113, 4.270008 ], x - 0.348955224))

        conditions.append((x >= 0.351940299) & (x < 0.354925373))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3299.553922406162, -5.225049574715127,
                            -2.050762054932407, 4.263752 ], x - 0.351940299))

        conditions.append((x >= 0.354925373) & (x < 0.357910448))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8532.037781534193, -34.773287450832655,
                            -2.1701600508306065, 4.257496 ], x - 0.354925373))

        conditions.append((x >= 0.357910448) & (x < 0.360895522))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7933.178960177786, 41.6330305913077,
                            -2.149683203075555, 4.250935 ], x - 0.357910448))

        conditions.append((x >= 0.360895522) & (x < 0.363880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ 6094.933167713825, -29.410348162813264,
                            -2.1131975915480017, 4.244678 ], x - 0.360895522))

        conditions.append((x >= 0.363880597) & (x < 0.366865672))
        eval_funcs.append(lambda x: np.polyval(
                          [ 659.2546815609609, 25.171149714025507,
                            -2.1258519168575156, 4.23827 ], x - 0.363880597))

        conditions.append((x >= 0.366865672) & (x < 0.369850746))
        eval_funcs.append(lambda x: np.polyval(
                          [ -14408.932506343654, 31.074923719706895,
                            -1.9579531692023153, 4.232166 ], x - 0.366865672))

        conditions.append((x >= 0.369850746) & (x < 0.372835821))
        eval_funcs.append(lambda x: np.polyval(
                          [ 22539.250522049388, -97.96026565761849,
                            -2.1576108644022867, 4.226215 ], x - 0.369850746))

        conditions.append((x >= 0.372835821) & (x < 0.375820896))
        eval_funcs.append(lambda x: np.polyval(
                          [ -18377.622455818673, 103.88379409870134,
                            -2.139928687741021, 4.219501 ], x - 0.372835821))

        conditions.append((x >= 0.375820896) & (x < 0.37880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ 16533.81610528034, -60.69194995820771,
                            -2.010997793593336, 4.21355 ], x - 0.375820896))

        conditions.append((x >= 0.37880597) & (x < 0.381791045))
        eval_funcs.append(lambda x: np.polyval(
                          [ -24749.099646326074, 87.37204377175249,
                            -1.9313557392329652, 4.207446 ], x - 0.37880597))

        conditions.append((x >= 0.381791045) & (x < 0.384776119))
        eval_funcs.append(lambda x: np.polyval(
                          [ 25016.712552582278, -134.26171210851814,
                            -2.071324915943336, 4.201801 ], x - 0.381791045))

        conditions.append((x >= 0.384776119) & (x < 0.387761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6442.838370811419, 89.76850251004466,
                            -2.2041404390922894, 4.195087 ], x - 0.384776119))

        conditions.append((x >= 0.387761194) & (x < 0.390746269))
        eval_funcs.append(lambda x: np.polyval(
                          [ -16464.11223456121, 32.07143526079611,
                            -1.8404390868510032, 4.189136 ], x - 0.387761194))

        conditions.append((x >= 0.390746269) & (x < 0.393731343))
        eval_funcs.append(lambda x: np.polyval(
                          [ 26395.273740935594, -115.36839422495214,
                            -2.0890867566309326, 4.18349 ], x - 0.390746269))

        conditions.append((x >= 0.393731343) & (x < 0.396716418))
        eval_funcs.append(lambda x: np.polyval(
                          [ -20204.476694735167, 121.00714187589907,
                            -2.07225467762553, 4.176928 ], x - 0.393731343))

        conditions.append((x >= 0.396716418) & (x < 0.399701493))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8443.428115266686, -59.92849293271051,
                            -1.8899303296314423, 4.171283 ], x - 0.396716418))

        conditions.append((x >= 0.399701493) & (x < 0.402686567))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2102.725302183366, 15.68430561082965,
                            -2.022002547101309, 4.165332 ], x - 0.399701493))

        conditions.append((x >= 0.402686567) & (x < 0.405671642))
        eval_funcs.append(lambda x: np.polyval(
                          [ 42.81021760242602, -3.1460662752390065,
                            -1.9845749748548613, 4.15938 ], x - 0.402686567))

        conditions.append((x >= 0.405671642) & (x < 0.408656716))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1893.7390621114525, -2.762691144310034,
                            -2.002213058909023, 4.153429 ], x - 0.405671642))

        conditions.append((x >= 0.408656716) & (x < 0.411641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1865.5231769699733, 14.196162566969068,
                            -1.9680833006354985, 4.147478 ], x - 0.408656716))

        conditions.append((x >= 0.411641791) & (x < 0.414626866))
        eval_funcs.append(lambda x: np.polyval(
                          [ -183.81184609088265, -2.510017225512028,
                            -1.9331992803303482, 4.14168 ], x - 0.411641791))

        conditions.append((x >= 0.414626866) & (x < 0.41761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2563.10341611601, -4.1560936649207605,
                            -1.9530981212966085, 4.135882 ], x - 0.414626866))

        conditions.append((x >= 0.41761194) & (x < 0.420597015))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4278.765875923847, 18.797066435356733,
                            -1.9093937341448717, 4.130083 ], x - 0.41761194))

        conditions.append((x >= 0.420597015) & (x < 0.42358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ 3085.3015850896486, -19.520244705864037,
                            -1.9115524755207045, 4.124437 ], x - 0.420597015))

        conditions.append((x >= 0.42358209) & (x < 0.426567164))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2385.614691715088, 8.10932518147088,
                            -1.945614926119984, 4.118639 ], x - 0.42358209))

        conditions.append((x >= 0.426567164) & (x < 0.429552239))
        eval_funcs.append(lambda x: np.polyval(
                          [ 6532.489978031815, -13.254383989298502,
                            -1.960973307395703, 4.11284 ], x - 0.426567164))

        conditions.append((x >= 0.429552239) & (x < 0.432537313))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12315.501007332612, 45.245533574221724,
                            -1.865477326548489, 4.107042 ], x - 0.429552239))

        conditions.append((x >= 0.432537313) & (x < 0.435522388))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14044.371417358027, -65.04251198766475,
                            -1.924572772089018, 4.101549 ], x - 0.432537313))

        conditions.append((x >= 0.435522388) & (x < 0.438507463))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9424.704209332518, 60.72799403834514,
                            -1.9374519317565824, 4.095598 ], x - 0.435522388))

        conditions.append((x >= 0.438507463) & (x < 0.441492537))
        eval_funcs.append(lambda x: np.polyval(
                          [ 6398.107202203057, -23.67235271467512,
                            -1.826838063232327, 4.090105 ], x - 0.438507463))

        conditions.append((x >= 0.441492537) & (x < 0.444477612))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10340.294555006047, 33.62411766085175,
                            -1.7971313084373806, 4.084611 ], x - 0.441492537))

        conditions.append((x >= 0.444477612) & (x < 0.447462687))
        eval_funcs.append(lambda x: np.polyval(
                          [ 11917.04610434036, -58.97554664550283,
                            -1.872807225313737, 4.079271 ], x - 0.444477612))

        conditions.append((x >= 0.447462687) & (x < 0.450447761))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2815.412337059408, 47.74428255423877,
                            -1.9063333909709672, 4.073472 ], x - 0.447462687))

        conditions.append((x >= 0.450447761) & (x < 0.453432836))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12234.63812102149, 22.531640054333707,
                            -1.6965545615661084, 4.068132 ], x - 0.450447761))

        conditions.append((x >= 0.453432836) & (x < 0.45641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ 17429.2556094758, -87.03229711299106,
                            -1.8890938604354803, 4.062943 ], x - 0.453432836))

        conditions.append((x >= 0.45641791) & (x < 0.459402985))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5863.77470573096, 69.05055616460885,
                            -1.942770687815231, 4.056992 ], x - 0.45641791))

        conditions.append((x >= 0.459402985) & (x < 0.46238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5516.003604493343, 16.53913432547932,
                            -1.687279042475531, 4.051652 ], x - 0.459402985))

        conditions.append((x >= 0.46238806) & (x < 0.465373134))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5069.742150285387, -32.857919053569546,
                            -1.7359918387977356, 4.046616 ], x - 0.46238806))

        conditions.append((x >= 0.465373134) & (x < 0.468358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3371.427642842522, 12.542747384993799,
                            -1.7966341295511379, 4.041276 ], x - 0.465373134))

        conditions.append((x >= 0.468358209) & (x < 0.471343284))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8453.498516684724, -17.6491457278804,
                            -1.811877111584531, 4.035935 ], x - 0.468358209))

        conditions.append((x >= 0.471343284) & (x < 0.474328358))
        eval_funcs.append(lambda x: np.polyval(
                          [ -13223.954974470604, 58.05383552619641,
                            -1.6912660821848253, 4.030594 ], x - 0.471343284))

        conditions.append((x >= 0.474328358) & (x < 0.477313433))
        eval_funcs.append(lambda x: np.polyval(
                          [ 10042.667003298588, -60.36961698819398,
                            -1.6981788612167148, 4.025711 ], x - 0.474328358))

        conditions.append((x >= 0.477313433) & (x < 0.480298507))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4051.2665697111183, 29.564725626421005,
                            -1.7901337722984605, 4.020371 ], x - 0.477313433))

        conditions.append((x >= 0.480298507) & (x < 0.483283582))
        eval_funcs.append(lambda x: np.polyval(
                          [ 448.02831254865794, -6.715265886519972,
                            -1.7219264441148363, 4.015183 ], x - 0.480298507))

        conditions.append((x >= 0.483283582) & (x < 0.486268657))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2259.0860188216384, -2.7030715412764863,
                            -1.7500408877121154, 4.009995 ], x - 0.483283582))

        conditions.append((x >= 0.486268657) & (x < 0.489253731))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3769.94151521248, 17.527552051625893,
                            -1.7057887015526851, 4.004807 ], x - 0.486268657))

        conditions.append((x >= 0.489253731) & (x < 0.492238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7143.901082451734, -16.233111144118595,
                            -1.7019246996551487, 3.999771 ], x - 0.489253731))

        conditions.append((x >= 0.492238806) & (x < 0.495223881))
        eval_funcs.append(lambda x: np.polyval(
                          [ -13414.332591094864, 47.74213042697942,
                            -1.6078679139193655, 3.994736 ], x - 0.492238806))

        conditions.append((x >= 0.495223881) & (x < 0.498208955))
        eval_funcs.append(lambda x: np.polyval(
                          [ 12188.797567243288, -72.38623615110794,
                            -1.6814324178138194, 3.990005 ], x - 0.495223881))

        conditions.append((x >= 0.498208955) & (x < 0.50119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1016.1382496677732, 36.76715197661725,
                            -1.7877580198869052, 3.984665 ], x - 0.498208955))

        conditions.append((x >= 0.50119403) & (x < 0.504179104))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8049.17892051984, 27.667405319736236,
                            -1.5954160337654957, 3.979629 ], x - 0.50119403))

        conditions.append((x >= 0.504179104) & (x < 0.507164179))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4490.103863916036, -44.41477883123986,
                            -1.645408183002975, 3.974899 ], x - 0.504179104))

        conditions.append((x >= 0.507164179) & (x < 0.510149254))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7307.378777596478, -4.204888456502246,
                            -1.790541536331932, 3.969711 ], x - 0.507164179))

        conditions.append((x >= 0.510149254) & (x < 0.513134328))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10786.506836528495, 61.23433265709917,
                            -1.6203043681848348, 3.964523 ], x - 0.510149254))

        conditions.append((x >= 0.513134328) & (x < 0.516119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1438.9905907832442, -35.36123066853167,
                            -1.5430712441394143, 3.959945 ], x - 0.513134328))

        conditions.append((x >= 0.516119403) & (x < 0.519104478))
        eval_funcs.append(lambda x: np.polyval(
                          [ 10782.573775435805, -22.474746155185,
                            -1.7157159726564697, 3.955062 ], x - 0.516119403))

        conditions.append((x >= 0.519104478) & (x < 0.522089552))
        eval_funcs.append(lambda x: np.polyval(
                          [ -15921.67876910849, 74.08562808294191,
                            -1.5616536192859696, 3.950027 ], x - 0.519104478))

        conditions.append((x >= 0.522089552) & (x < 0.525074627))
        eval_funcs.append(lambda x: np.polyval(
                          [ 12752.38065000946, -68.4965399071072,
                            -1.5449697774885793, 3.945602 ], x - 0.522089552))

        conditions.append((x >= 0.525074627) & (x < 0.528059701))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12154.794990734652, 45.70389809937455,
                            -1.6130075227327987, 3.940719 ], x - 0.525074627))

        conditions.append((x >= 0.528059701) & (x < 0.531044776))
        eval_funcs.append(lambda x: np.polyval(
                          [ 12971.323333167576, -63.1449894071433,
                            -1.6650704709272461, 3.935988 ], x - 0.528059701))

        conditions.append((x >= 0.531044776) & (x < 0.534029851))
        eval_funcs.append(lambda x: np.polyval(
                          [ 421.28645950633285, 53.016129589122,
                            -1.6953058771485248, 3.9308 ], x - 0.531044776))

        conditions.append((x >= 0.534029851) & (x < 0.537014925))
        eval_funcs.append(lambda x: np.polyval(
                          [ -20446.213964994102, 56.78884462345478,
                            -1.3675297937509172, 3.926223 ], x - 0.534029851))

        conditions.append((x >= 0.537014925) & (x < 0.54))
        eval_funcs.append(lambda x: np.polyval(
                          [ 29745.20841321683, -126.31154049256895,
                            -1.5750601855997215, 3.922103 ], x - 0.537014925))

        conditions.append((x >= 0.54) & (x < 0.542985075))
        eval_funcs.append(lambda x: np.polyval(
                          [ -29659.969928656123, 140.0634935196807,
                            -1.5340095744173148, 3.917067 ], x - 0.54))

        conditions.append((x >= 0.542985075) & (x < 0.545970149))
        eval_funcs.append(lambda x: np.polyval(
                          [ 25772.002026277245, -125.54821068466963,
                            -1.4906803665085928, 3.912947 ], x - 0.542985075))

        conditions.append((x >= 0.545970149) & (x < 0.548955224))
        eval_funcs.append(lambda x: np.polyval(
                          [ -16057.351274902258, 105.24578884508671,
                            -1.5512845980789622, 3.908064 ], x - 0.545970149))

        conditions.append((x >= 0.548955224) & (x < 0.551940299))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4057.591805558886, -38.55140472570019,
                            -1.352196859403784, 3.903944 ], x - 0.548955224))

        conditions.append((x >= 0.551940299) & (x < 0.554925373))
        eval_funcs.append(lambda x: np.polyval(
                          [ -210.6674283383761, -2.2147571487639386,
                            -1.4738869100612004, 3.899672 ], x - 0.551940299))

        conditions.append((x >= 0.554925373) & (x < 0.557910448))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2537.2792337451847, -4.101330737703725,
                            -1.4927408997928087, 3.895247 ], x - 0.554925373))

        conditions.append((x >= 0.557910448) & (x < 0.560895522))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4224.064617395509, 18.620575688312645,
                            -1.4493998646718715, 3.890822 ], x - 0.557910448))

        conditions.append((x >= 0.560895522) & (x < 0.563880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8682.180296529905, -19.206860702810182,
                            -1.4511499688252363, 3.886549 ], x - 0.560895522))

        conditions.append((x >= 0.563880597) & (x < 0.566865672))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7646.7316312048415, 58.54401734318209,
                            -1.3337256059669784, 3.882277 ], x - 0.563880597))

        conditions.append((x >= 0.566865672) & (x < 0.569850746))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6667.781674492003, -9.934184928874242,
                            -1.1886216104728387, 3.878614 ], x - 0.566865672))

        conditions.append((x >= 0.569850746) & (x < 0.572835821))
        eval_funcs.append(lambda x: np.polyval(
                          [ 17024.07721627204, -69.64565007148039,
                            -1.4261733068566806, 3.8748 ], x - 0.569850746))

        conditions.append((x >= 0.572835821) & (x < 0.575820896))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9810.141359169293, 82.8087918176095,
                            -1.3868803415088542, 3.870375 ], x - 0.572835821))

        conditions.append((x >= 0.575820896) & (x < 0.57880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6393.632396404129, -5.043231335557522,
                            -1.1547443110528925, 3.866712 ], x - 0.575820896))

        conditions.append((x >= 0.57880597) & (x < 0.581791045))
        eval_funcs.append(lambda x: np.polyval(
                          [ 18128.481823320406, -62.299628831749175,
                            -1.3557677320239576, 3.86305 ], x - 0.57880597))

        conditions.append((x >= 0.581791045) & (x < 0.584776119))
        eval_funcs.append(lambda x: np.polyval(
                          [ -20216.447245925327, 100.04500480449572,
                            -1.2430949538421114, 3.85893 ], x - 0.581791045))

        conditions.append((x >= 0.584776119) & (x < 0.587761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ 16871.01755682755, -80.99776833405603,
                            -1.1862375434823496, 3.855573 ], x - 0.584776119))

        conditions.append((x >= 0.587761194) & (x < 0.590746269))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12905.509851774299, 70.0859898662851,
                            -1.2188100205920307, 3.851759 ], x - 0.587761194))

        conditions.append((x >= 0.590746269) & (x < 0.593731343))
        eval_funcs.append(lambda x: np.polyval(
                          [ 11817.810306968111, -45.48575459607036,
                            -1.1453764732927947, 3.848402 ], x - 0.590746269))

        conditions.append((x >= 0.593731343) & (x < 0.596716418))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11394.85765728179, 60.34536025671468,
                            -1.1010194507849542, 3.844892 ], x - 0.593731343))

        conditions.append((x >= 0.596716418) & (x < 0.599701493))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5076.316893320174, -41.69815390721684,
                            -1.0453561412912267, 3.84184 ], x - 0.596716418))

        conditions.append((x >= 0.599701493) & (x < 0.602686567))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2556.135014168678, 3.7614060437665477,
                            -1.1586001789197162, 3.838483 ], x - 0.599701493))

        conditions.append((x >= 0.602686567) & (x < 0.605671642))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3834.1830342811654, 26.65216255762038,
                            -1.0678134260404986, 3.835126 ], x - 0.602686567))

        conditions.append((x >= 0.605671642) & (x < 0.608656716))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1351.5254683541073, -7.683809205550682,
                            -1.0111914686580679, 3.832074 ], x - 0.605671642))

        conditions.append((x >= 0.608656716) & (x < 0.611641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1647.0304457881111, 4.419401402214233,
                            -1.0209359675172043, 3.829023 ], x - 0.608656716))

        conditions.append((x >= 0.611641791) & (x < 0.614626866))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5274.153659581794, -10.330126821668475,
                            -1.0385799261986821, 3.825971 ], x - 0.611641791))

        conditions.append((x >= 0.614626866) & (x < 0.61761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2230.938923062282, 36.901105884459625,
                            -0.9592635608728196, 3.822919 ], x - 0.614626866))

        conditions.append((x >= 0.61761194) & (x < 0.620597015))
        eval_funcs.append(lambda x: np.polyval(
                          [ -13569.01165711545, 16.92255255999637,
                            -0.7985959574653975, 3.820325 ], x - 0.61761194))

        conditions.append((x >= 0.620597015) & (x < 0.62358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ 27821.68911959815, -104.59099985709562,
                            -1.060292847780786, 3.817731 ], x - 0.620597015))

        conditions.append((x >= 0.62358209) & (x < 0.626567164))
        eval_funcs.append(lambda x: np.polyval(
                          [ -23128.574967014978, 144.5584860889582,
                            -0.9409869038172091, 3.814374 ], x - 0.62358209))

        conditions.append((x >= 0.626567164) & (x < 0.629552239))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7284.509714642538, -62.56303728430586,
                            -0.6962244214721074, 3.812238 ], x - 0.626567164))

        conditions.append((x >= 0.629552239) & (x < 0.632537313))
        eval_funcs.append(lambda x: np.polyval(
                          [ -219.84620694244964, 2.6713862250042864,
                            -0.8750054917579532, 3.809796 ], x - 0.629552239))

        conditions.append((x >= 0.632537313) & (x < 0.635522388))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5061.525159895333, 0.7026146359767296,
                            -0.864933849511861, 3.807202 ], x - 0.632537313))

        conditions.append((x >= 0.635522388) & (x < 0.638507463))
        eval_funcs.append(lambda x: np.polyval(
                          [ -14311.792471707124, 46.02971128600058,
                            -0.7254343517103148, 3.804761 ], x - 0.635522388))

        conditions.append((x >= 0.638507463) & (x < 0.641492537))
        eval_funcs.append(lambda x: np.polyval(
                          [ 17748.26221766559, -82.13561045144291,
                            -0.8332131686615979, 3.802625 ], x - 0.638507463))

        conditions.append((x >= 0.641492537) & (x < 0.644477612))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10739.626454812706, 76.80401782196078,
                            -0.849128367198457, 3.799878 ], x - 0.641492537))

        conditions.append((x >= 0.644477612) & (x < 0.647462687))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2201.8514047952, -19.37175349683933,
                            -0.677688750768145, 3.797742 ], x - 0.644477612))

        conditions.append((x >= 0.647462687) & (x < 0.650447761))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2007.3860394852336, 0.34632124966787714,
                            -0.7344810929333708, 3.795605 ], x - 0.647462687))

        conditions.append((x >= 0.650447761) & (x < 0.653432836))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4554.438900694718, 18.322908872959086,
                            -0.6787520594942996, 3.793469 ], x - 0.650447761))

        conditions.append((x >= 0.653432836) & (x < 0.65641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4781.326016879791, -22.46311623151482,
                            -0.6911108889751404, 3.791485 ], x - 0.653432836))

        conditions.append((x >= 0.65641791) & (x < 0.659402985))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3104.2182795756626, 20.354719704019686,
                            -0.6974046086310559, 3.789349 ], x - 0.65641791))

        conditions.append((x >= 0.659402985) & (x < 0.66238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1845.8355399829754, -7.444253438693328,
                            -0.6588658985440867, 3.787366 ], x - 0.659402985))

        conditions.append((x >= 0.66238806) & (x < 0.665373134))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1510.5406698024829, 9.085619134850102,
                            -0.6539662988386318, 3.785382 ], x - 0.66238806))

        conditions.append((x >= 0.665373134) & (x < 0.668358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7925.5449291116, 22.6128461729599,
                            -0.5593440342083844, 3.783551 ], x - 0.665373134))

        conditions.append((x >= 0.668358209) & (x < 0.671343284))
        eval_funcs.append(lambda x: np.polyval(
                          [ 13048.129889377411, -48.36219191484335,
                            -0.6362077624488379, 3.781872 ], x - 0.668358209))

        conditions.append((x >= 0.671343284) & (x < 0.674328358))
        eval_funcs.append(lambda x: np.polyval(
                          [ -15656.927005509533, 68.48674707375665,
                            -0.5761344559578448, 3.779889 ], x - 0.671343284))

        conditions.append((x >= 0.674328358) & (x < 0.677313433))
        eval_funcs.append(lambda x: np.polyval(
                          [ 15179.835120446241, -71.724510098378,
                            -0.5857994181808034, 3.778363 ], x - 0.674328358))

        conditions.append((x >= 0.677313433) & (x < 0.680298507))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10625.103681740558, 64.21432886812039,
                            -0.6082178724167147, 3.776379 ], x - 0.677313433))

        conditions.append((x >= 0.680298507) & (x < 0.683283582))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4387.430097681434, -30.935833374884808,
                            -0.5088791007607385, 3.774853 ], x - 0.680298507))

        conditions.append((x >= 0.683283582) & (x < 0.686268657))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1247.7400790306128, 8.35459032162284,
                            -0.5762858048679517, 3.773175 ], x - 0.683283582))

        conditions.append((x >= 0.686268657) & (x < 0.689253731))
        eval_funcs.append(lambda x: np.polyval(
                          [ 678.6999451293758, -2.8192028276139767,
                            -0.5597622580442735, 3.771496 ], x - 0.686268657))

        conditions.append((x >= 0.689253731) & (x < 0.692238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4209.883198474608, 3.2587058524073544,
                            -0.5584503089920411, 3.769818 ], x - 0.689253731))

        conditions.append((x >= 0.692238806) & (x < 0.695223881))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11766.163565556068, 40.95915711846716,
                            -0.4264566716842576, 3.768292 ], x - 0.692238806))

        conditions.append((x >= 0.695223881) & (x < 0.698208955))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14207.087213553768, -64.40948499788996,
                            -0.4964576591789253, 3.767071 ], x - 0.695223881))

        conditions.append((x >= 0.698208955) & (x < 0.70119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10699.959360641938, 62.81813497284678,
                            -0.5012079567635811, 3.765393 ], x - 0.698208955))

        conditions.append((x >= 0.70119403) & (x < 0.704179104))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5697.139596567808, -33.00240859255807,
                            -0.4122057773389406, 3.764172 ], x - 0.70119403))

        conditions.append((x >= 0.704179104) & (x < 0.707164179))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6411.655402391743, 18.016741259695813,
                            -0.4569391032669161, 3.762799 ], x - 0.704179104))

        conditions.append((x >= 0.707164179) & (x < 0.710149254))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8558.072162724691, -39.40107549118777,
                            -0.5207729447729873, 3.761425 ], x - 0.707164179))

        conditions.append((x >= 0.710149254) & (x < 0.713134328))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4962.670270910762, 37.23838629224865,
                            -0.5272287342335108, 3.759747 ], x - 0.710149254))

        conditions.append((x >= 0.713134328) & (x < 0.716119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ -98.75913525832563, -7.203427696557833,
                            -0.43757216023843676, 3.758373 ], x - 0.713134328))

        conditions.append((x >= 0.716119403) & (x < 0.719104478))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5320.091629698413, -8.087837973601665,
                            -0.48321773510878796, 3.757 ], x - 0.716119403))

        conditions.append((x >= 0.719104478) & (x < 0.722089552))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9715.037829415383, 39.55477959096454,
                            -0.3892865543603388, 3.755627 ], x - 0.719104478))

        conditions.append((x >= 0.722089552) & (x < 0.725074627))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4854.809948991101, -47.44554090984925,
                            -0.4128410608135476, 3.754559 ], x - 0.722089552))

        conditions.append((x >= 0.725074627) & (x < 0.728059701))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7514.44406115188, -3.9696254843954755,
                            -0.5663191886378476, 3.753033 ], x - 0.725074627))

        conditions.append((x >= 0.728059701) & (x < 0.731044776))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11979.365018619374, 63.32388928979935,
                            -0.3891423189631999, 3.751507 ], x - 0.728059701))

        conditions.append((x >= 0.731044776) & (x < 0.734029851))
        eval_funcs.append(lambda x: np.polyval(
                          [ 288.74464021490786, -43.954019809066445,
                            -0.33132180582300114, 3.750591 ], x - 0.731044776))

        conditions.append((x >= 0.734029851) & (x < 0.737014925))
        eval_funcs.append(lambda x: np.polyval(
                          [ 16501.27122115043, -41.36824658839789,
                            -0.5860151701894124, 3.749218 ], x - 0.734029851))

        conditions.append((x >= 0.737014925) & (x < 0.74))
        eval_funcs.append(lambda x: np.polyval(
                          [ -20314.63321024551, 106.40430047921693,
                            -0.3918777366573277, 3.747539 ], x - 0.737014925))

        conditions.append((x >= 0.74) & (x < 0.742985075))
        eval_funcs.append(lambda x: np.polyval(
                          [ 13026.022987639257, -75.51781071100417,
                            -0.2996792482124799, 3.746777 ], x - 0.74))

        conditions.append((x >= 0.742985075) & (x < 0.745970149))
        eval_funcs.append(lambda x: np.polyval(
                          [ -14533.20858535682, 41.133155998477726,
                            -0.4023200213784749, 3.745556 ], x - 0.742985075))

        conditions.append((x >= 0.745970149) & (x < 0.748955224))
        eval_funcs.append(lambda x: np.polyval(
                          [ 22173.64761626816, -89.01495325569991,
                            -0.5452507294442818, 3.744335 ], x - 0.745970149))

        conditions.append((x >= 0.748955224) & (x < 0.751940299))
        eval_funcs.append(lambda x: np.polyval(
                          [ -16828.471629916636, 109.55505021869538,
                            -0.483936999502468, 3.742504 ], x - 0.748955224))

        conditions.append((x >= 0.751940299) & (x < 0.754925373))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5025.9741199434275, -41.14769963332496,
                            -0.27973592745384324, 3.741588 ], x - 0.751940299))

        conditions.append((x >= 0.754925373) & (x < 0.757910448))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3275.389529590394, 3.8610143770220127,
                            -0.3910394421586139, 3.74052 ], x - 0.754925373))

        conditions.append((x >= 0.757910448) & (x < 0.760895522))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8000.3679568408725, -25.470835823104117,
                            -0.45554637991177765, 3.7393 ], x - 0.757910448))

        conditions.append((x >= 0.760895522) & (x < 0.763880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11432.168272803627, 46.17423531209312,
                            -0.3937452003855827, 3.737926 ], x - 0.760895522))

        conditions.append((x >= 0.763880597) & (x < 0.766865672))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14757.522496994176, -56.20340380872483,
                            -0.4236830205356657, 3.736858 ], x - 0.763880597))

        conditions.append((x >= 0.766865672) & (x < 0.769850746))
        eval_funcs.append(lambda x: np.polyval(
                          [ -18950.27241739533, 75.95353059441999,
                            -0.3647274108208567, 3.735485 ], x - 0.766865672))

        conditions.append((x >= 0.769850746) & (x < 0.772835821))
        eval_funcs.append(lambda x: np.polyval(
                          [ 20929.30842178209, -93.7503658638337,
                            -0.417852281065867, 3.734569 ], x - 0.769850746))

        conditions.append((x >= 0.772835821) & (x < 0.775820896))
        eval_funcs.append(lambda x: np.polyval(
                          [ -18863.039320767097, 93.67630014762,
                            -0.4180733727836934, 3.733043 ], x - 0.772835821))

        conditions.append((x >= 0.775820896) & (x < 0.77880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14333.391923979367, -75.2464611536966,
                            -0.3630589211489078, 3.732128 ], x - 0.775820896))

        conditions.append((x >= 0.77880597) & (x < 0.781791045))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9822.804775224797, 53.112245538542396,
                            -0.429131192692097, 3.730755 ], x - 0.77880597))

        conditions.append((x >= 0.781791045) & (x < 0.784776119))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7851.910301433685, -34.853181354670284,
                            -0.3746265166734246, 3.729686 ], x - 0.781791045))

        conditions.append((x >= 0.784776119) & (x < 0.787761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10193.398194876783, 35.46241851875613,
                            -0.37280789865507846, 3.728466 ], x - 0.784776119))

        conditions.append((x >= 0.787761194) & (x < 0.790746269))
        eval_funcs.append(lambda x: np.polyval(
                          [ 9988.498243796039, -55.821755830959496,
                            -0.43358204748230383, 3.727398 ], x - 0.787761194))

        conditions.append((x >= 0.790746269) & (x < 0.793731343))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1075.3467367210662, 33.627493354338974,
                            -0.49983358554470186, 3.725872 ], x - 0.790746269))

        conditions.append((x >= 0.793731343) & (x < 0.796716418))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5687.082342508064, 23.99752460002619,
                            -0.3278186426995915, 3.724651 ], x - 0.793731343))

        conditions.append((x >= 0.796716418) & (x < 0.799701493))
        eval_funcs.append(lambda x: np.polyval(
                          [ 928.082526192252, -26.931577370660733,
                            -0.33657701027389325, 3.723735 ], x - 0.796716418))

        conditions.append((x >= 0.799701493) & (x < 0.802686567))
        eval_funcs.append(lambda x: np.polyval(
                          [ 7614.042079947063, -18.62038953004055,
                            -0.4725530478700048, 3.722515 ], x - 0.799701493))

        conditions.append((x >= 0.802686567) & (x < 0.805671642))
        eval_funcs.append(lambda x: np.polyval(
                          [ -14090.337468084816, 49.56504761322536,
                            -0.3801809535870025, 3.721141 ], x - 0.802686567))

        conditions.append((x >= 0.805671642) & (x < 0.808656716))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14309.919296366303, -76.61709473940451,
                            -0.46093334316218215, 3.720073 ], x - 0.805671642))

        conditions.append((x >= 0.808656716) & (x < 0.811641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3035.0233991345335, 51.531409361640904,
                            -0.5358159703555252, 3.718395 ], x - 0.808656716))

        conditions.append((x >= 0.811641791) & (x < 0.814626866))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7884.337885841697, 24.352091942126318,
                            -0.309298027701182, 3.717174 ], x - 0.811641791))

        conditions.append((x >= 0.814626866) & (x < 0.81761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5924.691421915478, -46.25392780161038,
                            -0.37467665037943176, 3.716258 ], x - 0.814626866))

        conditions.append((x >= 0.81761194) & (x < 0.820597015))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1366.6966384485345, 6.802999163139008,
                            -0.49244059173398946, 3.714885 ], x - 0.81761194))

        conditions.append((x >= 0.820597015) & (x < 0.82358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5676.99720473314, 19.042075067189188,
                            -0.41529110677589215, 3.713512 ], x - 0.820597015))

        conditions.append((x >= 0.82358209) & (x < 0.826567164))
        eval_funcs.append(lambda x: np.polyval(
                          [ 9912.288247028098, -31.796712225567344,
                            -0.4533646552914375, 3.712291 ], x - 0.82358209))

        conditions.append((x >= 0.826567164) & (x < 0.829552239))
        eval_funcs.append(lambda x: np.polyval(
                          [ -16791.02466384812, 56.97002955455776,
                            -0.3782204402389207, 3.710918 ], x - 0.826567164))

        conditions.append((x >= 0.829552239) & (x < 0.832537313))
        eval_funcs.append(lambda x: np.polyval(
                          [ 17099.925835909038, -93.39737429075157,
                            -0.4869587963273149, 3.70985 ], x - 0.829552239))

        conditions.append((x >= 0.832537313) & (x < 0.835522388))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5742.263616462703, 59.73625775335111,
                            -0.5874397201140802, 3.708019 ], x - 0.832537313))

        conditions.append((x >= 0.835522388) & (x < 0.838507463))
        eval_funcs.append(lambda x: np.polyval(
                          [ 192.21082641628635, 8.312995058613831,
                            -0.3843075967764037, 3.706645 ], x - 0.835522388))

        conditions.append((x >= 0.838507463) & (x < 0.841492537))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12282.872824337192, 10.03428625660763,
                            -0.3295395860043689, 3.705577 ], x - 0.838507463))

        conditions.append((x >= 0.841492537) & (x < 0.844477612))
        eval_funcs.append(lambda x: np.polyval(
                          [ 26043.715863581154, -99.9615666831002,
                            -0.5979791726962034, 3.704356 ], x - 0.841492537))

        conditions.append((x >= 0.844477612) & (x < 0.847462687))
        eval_funcs.append(lambda x: np.polyval(
                          [ -23130.087929857316, 133.26576871133852,
                            -0.4985636318267595, 3.702373 ], x - 0.844477612))

        conditions.append((x >= 0.847462687) & (x < 0.850447761))
        eval_funcs.append(lambda x: np.polyval(
                          [ 9181.286410897446, -73.86937297031825,
                            -0.3212609358101332, 3.701457 ], x - 0.847462687))

        conditions.append((x >= 0.850447761) & (x < 0.853432836))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2166.018036597907, 8.351085084849657,
                            -0.516837873501556, 3.700084 ], x - 0.850447761))

        conditions.append((x >= 0.853432836) & (x < 0.85641791))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5234.846924131092, -11.046093786942855,
                            -0.5248826766029567, 3.698558 ], x - 0.853432836))

        conditions.append((x >= 0.85641791) & (x < 0.859402985))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7306.737378307295, 35.83312255466875,
                            -0.45089156149116544, 3.697032 ], x - 0.85641791))

        conditions.append((x >= 0.859402985) & (x < 0.86238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1058.912705902023, -29.60035468398314,
                            -0.43228628193957896, 3.695811 ], x - 0.859402985))

        conditions.append((x >= 0.86238806) & (x < 0.865373134))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8823.159293023344, -20.117553147271707,
                            -0.580697965658962, 3.694285 ], x - 0.86238806))

        conditions.append((x >= 0.865373134) & (x < 0.868358209))
        eval_funcs.append(lambda x: np.polyval(
                          [ -13455.918043119089, 58.89579706311635,
                            -0.46494203798011513, 3.692607 ], x - 0.865373134))

        conditions.append((x >= 0.868358209) & (x < 0.871343284))
        eval_funcs.append(lambda x: np.polyval(
                          [ 10638.331183738397, -61.604976594574886,
                            -0.4730291420699839, 3.691386 ], x - 0.868358209))

        conditions.append((x >= 0.871343284) & (x < 0.874328358))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6201.851499146273, 33.663672780319196,
                            -0.556436029553324, 3.689708 ], x - 0.871343284))

        conditions.append((x >= 0.874328358) & (x < 0.877313433))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2702.517787993327, -21.87528420556699,
                            -0.5212468173169352, 3.688182 ], x - 0.874328358))

        conditions.append((x >= 0.877313433) & (x < 0.880298507))
        eval_funcs.append(lambda x: np.polyval(
                          [ 1181.4249257458926, 2.326370652415817,
                            -0.5796017904416089, 3.686503 ], x - 0.877313433))

        conditions.append((x >= 0.880298507) & (x < 0.883283582))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1751.2737489512617, 12.906293138803704,
                            -0.5341311618076969, 3.684825 ], x - 0.880298507))

        conditions.append((x >= 0.883283582) & (x < 0.886268657))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5605.341434401712, -2.7767573196484534,
                            -0.5038937376723319, 3.683299 ], x - 0.883283582))

        conditions.append((x >= 0.886268657) & (x < 0.889253731))
        eval_funcs.append(lambda x: np.polyval(
                          [ 12630.839293065927, -52.973851066538714,
                            -0.6703134850007295, 3.681621 ], x - 0.886268657))

        conditions.append((x >= 0.889253731) & (x < 0.892238806))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4690.90458770018, 60.138118849191045,
                            -0.6489276155136965, 3.679484 ], x - 0.889253731))

        conditions.append((x >= 0.892238806) & (x < 0.895223881))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11085.921435905046, 18.13001281280369,
                            -0.41529137239276737, 3.677958 ], x - 0.892238806))

        conditions.append((x >= 0.895223881) & (x < 0.898208955))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14597.213516966345, -81.14690797804906,
                            -0.6034015307281627, 3.676585 ], x - 0.895223881))

        conditions.append((x >= 0.898208955) & (x < 0.90119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7151.008935066567, 49.57437964778193,
                            -0.6976478641611041, 3.674449 ], x - 0.898208955))

        conditions.append((x >= 0.90119403) & (x < 0.904179104))
        eval_funcs.append(lambda x: np.polyval(
                          [ 2540.19148138613, -14.464514342749464,
                            -0.5928422829856848, 3.672618 ], x - 0.90119403))

        conditions.append((x >= 0.904179104) & (x < 0.907164179))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3009.7066147722444, 8.283464295572466,
                            -0.6112931747742119, 3.670787 ], x - 0.904179104))

        conditions.append((x >= 0.907164179) & (x < 0.910149254))
        eval_funcs.append(lambda x: np.polyval(
                          [ 3746.5224715824897, -18.669135623701226,
                            -0.6422951826140262, 3.668956 ], x - 0.907164179))

        conditions.append((x >= 0.910149254) & (x < 0.913134328))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6224.320541123043, 14.88181607687574,
                            -0.653600615510265, 3.666972 ], x - 0.910149254))

        conditions.append((x >= 0.913134328) & (x < 0.916119403))
        eval_funcs.append(lambda x: np.polyval(
                          [ 9684.21416688224, -40.85835616804177,
                            -0.7311425099463635, 3.664988 ], x - 0.913134328))

        conditions.append((x >= 0.916119403) & (x < 0.919104478))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9541.795187775857, 45.86596064457575,
                            -0.7161944350135722, 3.662699 ], x - 0.916119403))

        conditions.append((x >= 0.919104478) & (x < 0.922089552))
        eval_funcs.append(lambda x: np.polyval(
                          [ 5474.576109444154, -39.582962165871464,
                            -0.697439213329754, 3.660716 ], x - 0.919104478))

        conditions.append((x >= 0.922089552) & (x < 0.925074627))
        eval_funcs.append(lambda x: np.polyval(
                          [ -852.261986803286, 9.443082250097724,
                            -0.7874089852294531, 3.658427 ], x - 0.922089552))

        conditions.append((x >= 0.925074627) & (x < 0.928059701))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2065.584392080834, 1.8108843993271515,
                            -0.7538150507334208, 3.656138 ], x - 0.925074627))

        conditions.append((x >= 0.928059701) & (x < 0.931044776))
        eval_funcs.append(lambda x: np.polyval(
                          [ 3400.1588955625157, -16.686882391492052,
                            -0.7982210055638846, 3.653849 ], x - 0.928059701))

        conditions.append((x >= 0.931044776) & (x < 0.934029851))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5858.180033929496, 13.762305554022003,
                            -0.8069510867669961, 3.651408 ], x - 0.931044776))

        conditions.append((x >= 0.934029851) & (x < 0.937014925))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8641.133749113367, -38.69901474032463,
                            -0.8813890339412975, 3.648966 ], x - 0.934029851))

        conditions.append((x >= 0.937014925) & (x < 0.94))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11562.786748721248, 38.68425631467848,
                            -0.8814330889339742, 3.64622 ], x - 0.937014925))

        conditions.append((x >= 0.94) & (x < 0.942985075))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14676.809959597205, -64.86310064713548,
                            -0.9595789026796793, 3.643626 ], x - 0.94))

        conditions.append((x >= 0.942985075) & (x < 0.945970149))
        eval_funcs.append(lambda x: np.polyval(
                          [ -18459.232215357733, 66.57103482329839,
                            -0.9544805910687693, 3.640574 ], x - 0.942985075))

        conditions.append((x >= 0.945970149) & (x < 0.948955224))
        eval_funcs.append(lambda x: np.polyval(
                          [ 19045.907094356877, -98.73548761478398,
                            -1.0504938628208607, 3.637827 ], x - 0.945970149))

        conditions.append((x >= 0.948955224) & (x < 0.951940299))
        eval_funcs.append(lambda x: np.polyval(
                          [ -17647.777324835533, 71.82489574427831,
                            -1.1308239978487105, 3.634318 ], x - 0.948955224))

        conditions.append((x >= 0.951940299) & (x < 0.954925373))
        eval_funcs.append(lambda x: np.polyval(
                          [ 11468.495529970663, -86.21492094952195,
                            -1.1737793023382537, 3.631113 ], x - 0.951940299))

        conditions.append((x >= 0.954925373) & (x < 0.957910448))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5330.51087183316, 16.48800252737417,
                            -1.3819193136203305, 3.627146 ], x - 0.954925373))

        conditions.append((x >= 0.957910448) & (x < 0.960895522))
        eval_funcs.append(lambda x: np.polyval(
                          [ 4138.952507542343, -31.247921694837338,
                            -1.4259787793291474, 3.623026 ], x - 0.957910448))

        conditions.append((x >= 0.960895522) & (x < 0.963880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5510.685907667969, 5.817316857659868,
                            -1.5018910166328776, 3.618601 ], x - 0.960895522))

        conditions.append((x >= 0.963880597) & (x < 0.966865672))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5029.432637120602, -43.532115349835934,
                            -1.6144725187419107, 3.614023 ], x - 0.963880597))

        conditions.append((x >= 0.966865672) & (x < 0.969850746))
        eval_funcs.append(lambda x: np.polyval(
                          [ 8447.249611570098, -88.57181623759361,
                            -2.0088126623252593, 3.608682 ], x - 0.966865672))

        conditions.append((x >= 0.969850746) & (x < 0.972835821))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5863.811606501247, -12.924820676568523,
                            -2.311787634265169, 3.602121 ], x - 0.969850746))

        conditions.append((x >= 0.972835821) & (x < 0.975820896))
        eval_funcs.append(lambda x: np.polyval(
                          [ -19354.238571813446, -65.43657297039812,
                            -2.5457022714058897, 3.594949 ], x - 0.972835821))

        conditions.append((x >= 0.975820896) & (x < 0.97880597))
        eval_funcs.append(lambda x: np.polyval(
                          [ 14368.300366948208, -238.75813408466593,
                            -3.4537462865682875, 3.586252 ], x - 0.975820896))

        conditions.append((x >= 0.97880597) & (x < 0.981791045))
        eval_funcs.append(lambda x: np.polyval(
                          [ -32328.94959806069, -110.0868145359612,
                            -4.49507427272707, 3.574197 ], x - 0.97880597))

        conditions.append((x >= 0.981791045) & (x < 0.984776119))
        eval_funcs.append(lambda x: np.polyval(
                          [ -17012.681566062525, -399.59983220025407,
                            -6.016527139733181, 3.558938 ], x - 0.981791045))

        conditions.append((x >= 0.984776119) & (x < 0.987761194))
        eval_funcs.append(lambda x: np.polyval(
                          [ -192186.89769555847, -551.9521724396507,
                            -8.856980288431556, 3.536965 ], x - 0.984776119))

        conditions.append((x >= 0.987761194) & (x < 0.990746269))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11637.558571088133, -2273.0290833553604,
                            -17.28976121057386, 3.500496 ], x - 0.987761194))

        conditions.append((x >= 0.990746269) & (x < 0.993731343))
        eval_funcs.append(lambda x: np.polyval(
                          [ -237447.6533006623, -2377.2460388101313,
                            -31.17118122087204, 3.428321 ], x - 0.990746269))

        conditions.append((x >= 0.993731343) & (x <= 0.996716418))
        eval_funcs.append(lambda x: np.polyval(
                          [ -237447.65330065906, -4503.642487496604,
                            -51.71114265764886, 3.307774 ], x - 0.993731343))

        return np.piecewise(x, conditions, eval_funcs)
