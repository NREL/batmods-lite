from numpy import ndarray as _ndarray


class GraphiteSlow(object):

    def __init__(self, alpha_a: float, alpha_c: float, Li_max: float) -> None:
        """
        Computationally slow Graphite kinetic and transport properties.

        Differs from ``GraphiteFast`` because the equilibrium potential is
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

        Ds = 3e-14 * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15))

        if np.atleast_1d(x).size > 1:
            Ds = Ds * np.ones_like(x)

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

        i0 = 2.5 * 0.27 * np.exp(-30e6 / c.R * (1 / T - 1 / 303.15)) \
           * C_Li**self.alpha_a * (self.Li_max * x)**self.alpha_c \
           * (self.Li_max - self.Li_max * x)**self.alpha_a

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

        conditions.append((x >= 0.) & (x < 0.0037))
        eval_funcs.append(lambda x: np.polyval(
                          [  181010.33572359575, -2672.363009950459,
                             -18.69569376464478, 1.170454  ], x - 0.0))

        conditions.append((x >= 0.0037) & (x < 0.0074))
        eval_funcs.append(lambda x: np.polyval(
                          [  181010.33572359628, -663.148283418546,
                             -31.037085550110096, 1.073864  ], x - 0.0037))

        conditions.append((x >= 0.0074) & (x < 0.0111))
        eval_funcs.append(lambda x: np.polyval(
                          [ -148966.15554926347, 1346.0664431133728,
                             -28.51028835923924, 0.959117  ], x - 0.0074))

        conditions.append((x >= 0.0111) & (x < 0.0148))
        eval_funcs.append(lambda x: np.polyval(
                          [  134712.9325556226, -307.4578834834498,
                             -24.66743668860853, 0.864511  ], x - 0.0111))

        conditions.append((x >= 0.0148) & (x < 0.0185))
        eval_funcs.append(lambda x: np.polyval(
                          [ -100702.30813422689, 1187.8556678839682,
                             -21.409964886326637, 0.775856  ], x - 0.0148))

        conditions.append((x >= 0.0185) & (x < 0.0222))
        eval_funcs.append(lambda x: np.polyval(
                          [  30143.9575731349, 70.06004759405101,
                             -16.755676739057968, 0.7078  ], x - 0.0185))

        conditions.append((x >= 0.0222) & (x < 0.0259))
        eval_funcs.append(lambda x: np.polyval(
                          [ -16872.712729461055, 404.65797665585376,
                             -14.999220049333337, 0.64829  ], x - 0.0222))

        conditions.append((x >= 0.0259) & (x < 0.0296))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4842.118164967984, 217.3708653588359,
                             -12.697713333878985, 0.597478  ], x - 0.0259))

        conditions.append((x >= 0.0296) & (x < 0.0333))
        eval_funcs.append(lambda x: np.polyval(
                          [  18157.360146991865, 163.62335372769078,
                             -11.288034723258832, 0.553227  ], x - 0.0296))

        conditions.append((x >= 0.0333) & (x < 0.037))
        eval_funcs.append(lambda x: np.polyval(
                          [ -40641.8423922001, 365.1700513593008,
                             -9.331499124436963, 0.514621  ], x - 0.0333))

        conditions.append((x >= 0.037) & (x < 0.0407))
        eval_funcs.append(lambda x: np.polyval(
                          [  23884.078085068348, -85.95439919412051,
                             -8.298401211425796, 0.483035  ], x - 0.037))

        conditions.append((x >= 0.0407) & (x < 0.0444))
        eval_funcs.append(lambda x: np.polyval(
                          [ -15706.267867250335, 179.15886755013975,
                             -7.953544678508529, 0.452364  ], x - 0.0407))

        conditions.append((x >= 0.0444) & (x < 0.0481))
        eval_funcs.append(lambda x: np.polyval(
                          [  23857.977570454488, 4.8192942236626095,
                             -7.272825479945466, 0.424593  ], x - 0.0444))

        conditions.append((x >= 0.0481) & (x < 0.0518))
        eval_funcs.append(lambda x: np.polyval(
                          [ -34535.82147602317, 269.64284525570946,
                             -6.257315563871799, 0.398958  ], x - 0.0481))

        conditions.append(( x >= 0.0518) & (x < 0.0555))
        eval_funcs.append(lambda x: np.polyval(
                          [  20904.857027693295, -113.70477312814793,
                             -5.680344696999822, 0.377748  ], x - 0.0518))

        conditions.append((x >= 0.0555) & (x < 0.0592))
        eval_funcs.append(lambda x: np.polyval(
                          [ -9934.888888519139, 118.33913987924774,
                             -5.663197540020752, 0.356233  ], x - 0.0555))

        conditions.append((x >= 0.0592) & (x < 0.0629))
        eval_funcs.append(lambda x: np.polyval(
                          [  15833.889097527042, 8.061873216685747,
                             -5.195513791565799, 0.336396  ], x - 0.0592))

        conditions.append((x >= 0.0629) & (x < 0.0666))
        eval_funcs.append(lambda x: np.polyval(
                          [ -26294.671805381935, 183.81804219923714,
                             -4.4855581045268895, 0.318085  ], x - 0.0629))

        conditions.append((x >= 0.0666) & (x < 0.0703))
        eval_funcs.append(lambda x: np.polyval(
                          [  20069.533085401796, -108.05281484050391,
                             -4.205226763299572, 0.302673  ], x - 0.0666))

        conditions.append((x >= 0.0703) & (x < 0.074))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2772.2785726656616, 114.71900240745559,
                             -4.18056186930185, 0.286651  ], x - 0.0703))

        conditions.append((x >= 0.074) & (x < 0.115))
        eval_funcs.append(lambda x: np.polyval(
                          [ -708.4919240956782, 83.94671025086652,
                             -3.445498732466058, 0.272613  ], x - 0.074))

        conditions.append((x >= 0.115) & (x < 0.118))
        eval_funcs.append(lambda x: np.polyval(
                          [ -920.0013687165061, -3.1977964129019405,
                             -0.13479326510950837, 0.223632  ], x - 0.115))

        conditions.append((x >= 0.118) & (x < 0.122))
        eval_funcs.append(lambda x: np.polyval(
                          [  2123.8322167289843, -11.477808731350406,
                             -0.1788200805422654, 0.223174  ], x - 0.118))

        conditions.append((x >= 0.122) & (x < 0.126))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4880.256717969286, 14.008177869397421,
                             -0.16869860399007733, 0.222411  ], x - 0.122))

        conditions.append((x >= 0.126) & (x < 0.129))
        eval_funcs.append(lambda x: np.polyval(
                          [  7616.690192902428, -44.554902746234006,
                             -0.290885503497424, 0.221648  ], x - 0.126))

        conditions.append((x >= 0.129) & (x < 0.133))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3041.684449568075, 23.995308989887963,
                             -0.35256428476646234, 0.22058  ], x - 0.129))

        conditions.append((x >= 0.133) & (x < 0.137))
        eval_funcs.append(lambda x: np.polyval(
                          [  835.7677528966614, -12.504904404929063,
                             -0.30660266642662637, 0.219359  ], x - 0.133))

        conditions.append((x >= 0.137) & (x < 0.141))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2707.6365620194147, -2.4756913701691907,
                             -0.3665250495270192, 0.217986  ], x - 0.137))

        conditions.append((x >= 0.141) & (x < 0.144))
        eval_funcs.append(lambda x: np.polyval(
                          [  6873.976941685764, -34.96733011440202,
                             -0.516297135465303, 0.216307  ], x - 0.141))

        conditions.append((x >= 0.144) & (x < 0.148))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3943.1319198047695, 26.8984623607699,
                             -0.5405037387261993, 0.214629  ], x - 0.144))

        conditions.append((x >= 0.148) & (x < 0.152))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3264.8215813615097, -20.41912067688729,
                             -0.5145863719906691, 0.212645  ], x - 0.148))

        conditions.append((x >= 0.152) & (x < 0.155))
        eval_funcs.append(lambda x: np.polyval(
                          [  10863.893956014545, -59.5969796532254,
                             -0.8346507733111204, 0.210051  ], x - 0.152))

        conditions.append((x >= 0.155) & (x < 0.159))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1034.6718365963468, 38.178065950905605,
                             -0.8989075144180798, 0.207304  ], x - 0.155))

        conditions.append((x >= 0.159) & (x < 0.163))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11556.298792470894, 25.762003911749492,
                             -0.6431472349674593, 0.204253  ], x - 0.159))

        conditions.append((x >= 0.163) & (x < 0.166))
        eval_funcs.append(lambda x: np.polyval(
                          [  23499.365611752546, -112.91358159790153,
                             -0.9917535457120671, 0.201353  ], x - 0.163))

        conditions.append((x >= 0.166) & (x < 0.17))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10035.666990583222, 98.58070890787154,
                             -1.034752163782157, 0.197996  ], x - 0.166))

        conditions.append((x >= 0.17) & (x < 0.174))
        eval_funcs.append(lambda x: np.polyval(
                          [  3262.9804989802706, -21.84729497912729,
                             -0.7278185080671795, 0.194792  ], x - 0.17))

        conditions.append((x >= 0.174) & (x < 0.178))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5391.255005337232, 17.308471008635596,
                             -0.7459738039491457, 0.19174  ], x - 0.174))

        conditions.append((x >= 0.178) & (x < 0.181))
        eval_funcs.append(lambda x: np.polyval(
                          [  10345.856663238894, -47.38658905541122,
                             -0.8662862761362484, 0.188688  ], x - 0.178))

        conditions.append(( x >= 0.181) & (x < 0.185))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4664.800193355612, 45.726120913738995,
                             -0.8712676805612654, 0.185942  ], x - 0.181))

        conditions.append((x >= 0.185) & (x < 0.189))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1914.0594900911544, -10.251481406528502,
                             -0.7293691225324229, 0.18289  ], x - 0.185))

        conditions.append((x >= 0.189) & (x < 0.192))
        eval_funcs.append(lambda x: np.polyval(
                          [  4027.7498339137455, -33.22019528762244,
                             -0.9032558293090267, 0.179686  ], x - 0.189))

        conditions.append((x >= 0.192) & (x < 0.196))
        eval_funcs.append(lambda x: np.polyval(
                          [  4138.096415542864, 3.029553217601343,
                             -0.9938277555190902, 0.176786  ], x - 0.192))

        conditions.append((x >= 0.196) & (x < 0.2))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7892.758686515039, 52.68671020411581,
                             -0.7709627018322216, 0.173124  ], x - 0.196))

        conditions.append((x >= 0.2) & (x < 0.203))
        eval_funcs.append(lambda x: np.polyval(
                          [  10155.62436157858, -42.02639403406461,
                             -0.7283214371520172, 0.170378  ], x - 0.2))

        conditions.append((x >= 0.203) & (x < 0.207))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1576.1848304243954, 49.374225220142755,
                             -0.7062779435937828, 0.168089  ], x - 0.203))

        conditions.append((x >= 0.207) & (x < 0.211))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12056.18845794916, 30.46000725505008,
                             -0.3869410136930132, 0.165953  ], x - 0.207))

        conditions.append((x >= 0.211) & (x < 0.214))
        eval_funcs.append(lambda x: np.polyval(
                          [  22214.899743170252, -114.21425424034004,
                             -0.7219580016341731, 0.164121  ], x - 0.211))

        conditions.append((x >= 0.214) & (x < 0.218))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11496.133736384785, 85.71984344819253,
                             -0.8074412340106162, 0.161527  ], x - 0.214))

        conditions.append((x >= 0.218) & (x < 0.222))
        eval_funcs.append(lambda x: np.polyval(
                          [  7480.121957827935, -52.23376138842496,
                             -0.673496905771546, 0.158933  ], x - 0.218))

        conditions.append((x >= 0.222) & (x < 0.226))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11299.354094927237, 37.527702105510336,
                             -0.7323211429032045, 0.155882  ], x - 0.222))

        conditions.append((x >= 0.226) & (x < 0.229))
        eval_funcs.append(lambda x: np.polyval(
                          [  22258.758931461252, -98.06454703361639,
                             -0.9744685226156298, 0.15283  ], x - 0.226))

        conditions.append((x >= 0.229) & (x < 0.233))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5980.488733141848, 102.26428334953513,
                             -0.9618693136678738, 0.149625  ], x - 0.229))

        conditions.append((x >= 0.233) & (x < 0.237))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2151.5730090580137, 30.498418551832863,
                             -0.4308185060624013, 0.147031  ], x - 0.233))

        conditions.append((x >= 0.237) & (x < 0.24))
        eval_funcs.append(lambda x: np.polyval(
                          [  2414.966824419423, 4.679542443136912,
                             -0.2901066620825233, 0.145658  ], x - 0.237))

        conditions.append((x >= 0.24) & (x < 0.244))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3833.229517954336, 26.41424386291173,
                             -0.19682530316437724, 0.144895  ], x - 0.24))

        conditions.append((x >= 0.244) & (x < 0.248))
        eval_funcs.append(lambda x: np.polyval(
                          [  1177.7756583158327, -19.584510352540356,
                             -0.16950636912289171, 0.144285  ], x - 0.244))

        conditions.append((x >= 0.248) & (x < 0.251))
        eval_funcs.append(lambda x: np.polyval(
                          [  3518.8327076641463, -5.451202452750396,
                             -0.2696492203440546, 0.143369  ], x - 0.248))

        conditions.append((x >= 0.251) & (x < 0.255))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3126.5759819551536, 26.21829191622691,
                             -0.2073479519536249, 0.142606  ], x - 0.251))

        conditions.append((x >= 0.255) & (x < 0.259))
        eval_funcs.append(lambda x: np.polyval(
                          [  133.1089516621029, -11.300619867234952,
                             -0.14767726375765705, 0.141996  ], x - 0.255))

        conditions.append((x >= 0.259) & (x < 0.263))
        eval_funcs.append(lambda x: np.polyval(
                          [  219.1401753074067, -9.703312447289715,
                             -0.2316929930157558, 0.141233  ], x - 0.259))

        conditions.append((x >= 0.263) & (x < 0.266))
        eval_funcs.append(lambda x: np.polyval(
                          [  7298.70243075365, -7.073630343600879,
                             -0.298800764179318, 0.140165  ], x - 0.263))

        conditions.append((x >= 0.266) & (x < 0.27))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8033.199095134383, 58.614691533182,
                             -0.14417758061057442, 0.139402  ], x - 0.266))

        conditions.append((x >= 0.27) & (x < 0.274))
        eval_funcs.append(lambda x: np.polyval(
                          [  3718.0247090807597, -37.78369760843068,
                             -0.060853604911569124, 0.139249  ], x - 0.27))

        conditions.append((x >= 0.274) & (x < 0.277))
        eval_funcs.append(lambda x: np.polyval(
                          [  1277.059597205721, 6.832598900538448,
                             -0.18465799974313804, 0.138639  ], x - 0.274))

        conditions.append((x >= 0.277) & (x < 0.281))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2523.2964928879237, 18.32613527538993,
                             -0.1091817972153528, 0.138181  ], x - 0.277))

        conditions.append((x >= 0.281) & (x < 0.285))
        eval_funcs.append(lambda x: np.polyval(
                          [  1062.7898267447642, -11.953422639265167,
                             -0.08369094667085378, 0.137876  ], x - 0.281))

        conditions.append((x >= 0.285) & (x < 0.288))
        eval_funcs.append(lambda x: np.polyval(
                          [  2693.064843282718, 0.80005528167184,
                             -0.1283044161012265, 0.137418  ], x - 0.285))

        conditions.append((x >= 0.288) & (x < 0.292))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3084.9513651439543, 25.037638871216306,
                             -0.050791333642561956, 0.137113  ], x - 0.288))

        conditions.append((x >= 0.292) & (x < 0.296))
        eval_funcs.append(lambda x: np.polyval(
                          [  515.3123901114132, -11.981777510511177,
                             0.0014321118002585999, 0.137113  ], x - 0.292))

        conditions.append((x >= 0.296) & (x < 0.3))
        eval_funcs.append(lambda x: np.polyval(
                          [  1039.3268046987487, -5.798028829174207,
                             -0.06968711355848302, 0.13696  ], x - 0.296))

        conditions.append((x >= 0.3) & (x < 0.303))
        eval_funcs.append(lambda x: np.polyval(
                          [ -500.5208424397201, 6.67389282721079,
                             -0.06618365756633668, 0.136655  ], x - 0.3))

        conditions.append((x >= 0.303) & (x < 0.307))
        eval_funcs.append(lambda x: np.polyval(
                          [ -454.5286020045236, 2.1692052452533055,
                             -0.039654363348944364, 0.136503  ], x - 0.303))

        conditions.append((x >= 0.307) & (x < 0.311))
        eval_funcs.append(lambda x: np.polyval(
                          [  1188.0403873964, -3.2851379788009822,
                             -0.04411809428313508, 0.13635  ], x - 0.307))

        conditions.append((x >= 0.311) & (x < 0.314))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2171.1978323724275, 10.971346669955835,
                             -0.013373259518515652, 0.136197  ], x - 0.311))

        conditions.append((x >= 0.314) & (x < 0.318))
        eval_funcs.append(lambda x: np.polyval(
                          [  152.82851615106762, -8.569433821396027,
                             -0.006167520972836227, 0.136197  ], x - 0.314))

        conditions.append((x >= 0.318) & (x < 0.322))
        eval_funcs.append(lambda x: np.polyval(
                          [  1129.9493299428868, -6.735491627583205,
                             -0.06738722276875324, 0.136045  ], x - 0.318))

        conditions.append((x >= 0.322) & (x < 0.325))
        eval_funcs.append(lambda x: np.polyval(
                          [ -493.1236714477341, 6.8239003317314655,
                             -0.06703358795216025, 0.13574  ], x - 0.322))

        conditions.append((x >= 0.325) & (x < 0.329))
        eval_funcs.append(lambda x: np.polyval(
                          [ -508.66400399689894, 2.385787288701845,
                             -0.039404525090860267, 0.135587  ], x - 0.325))

        conditions.append((x >= 0.329) & (x < 0.333))
        eval_funcs.append(lambda x: np.polyval(
                          [  1334.8013756339924, -3.7181807592609504,
                             -0.044734098973096686, 0.135435  ], x - 0.329))

        conditions.append((x >= 0.333) & (x < 0.337))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2424.291498539706, 12.299435748346973,
                             -0.010409079016752566, 0.135282  ], x - 0.333))

        conditions.append((x >= 0.337) & (x < 0.34))
        eval_funcs.append(lambda x: np.polyval(
                          [  3083.974629140745, -16.792062234129535,
                             -0.028379584959882808, 0.135282  ], x - 0.337,))

        conditions.append((x >= 0.34) & (x < 0.344))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2249.3871459178845, 10.963709428137191,
                             -0.045864643377859844, 0.135129  ], x - 0.34))

        conditions.append((x >= 0.344) & (x < 0.348))
        eval_funcs.append(lambda x: np.polyval(
                          [  3374.456015520031, -16.02893632287706,
                             -0.0661255509568191, 0.134977  ], x - 0.344))

        conditions.append((x >= 0.348) & (x < 0.351))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4556.71719946846, 24.464535863363352,
                             -0.0323831527948739, 0.134672  ], x - 0.348))

        conditions.append((x >= 0.351) & (x < 0.355))
        eval_funcs.append(lambda x: np.polyval(
                          [  2285.0611079848145, -16.545918931852828,
                             -0.008627302000342311, 0.134672  ], x - 0.351))

        conditions.append((x >= 0.355) & (x < 0.359))
        eval_funcs.append(lambda x: np.polyval(
                          [ -761.7210739978832, 10.874814363964969,
                             -0.031311720271893766, 0.134519  ], x - 0.355))

        conditions.append((x >= 0.359) & (x < 0.362))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2702.9630573220747, 1.7341614759903656,
                             0.019124183087927613, 0.134519  ], x - 0.359))

        conditions.append((x >= 0.362) & (x < 0.366))
        eval_funcs.append(lambda x: np.polyval(
                          [  3598.1796727162355, -22.59250603990833,
                             -0.04345085060382633, 0.134519  ], x - 0.362))

        conditions.append((x >= 0.366) & (x < 0.37))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1929.0203436270383, 20.58565003268653,
                             -0.051478274632713525, 0.134214  ], x - 0.366))

        conditions.append((x >= 0.37) & (x < 0.374))
        eval_funcs.append(lambda x: np.polyval(
                          [ -647.7232982080634, -2.5625940908379534,
                             0.020613949134680856, 0.134214  ], x - 0.37))

        conditions.append((x >= 0.374) & (x < 0.377))
        eval_funcs.append(lambda x: np.polyval(
                          [  1220.3714348909568, -10.33527366933472,
                             -0.03097752190600989, 0.134214  ], x - 0.374))

        conditions.append((x >= 0.377) & (x < 0.381))
        eval_funcs.append(lambda x: np.polyval(
                          [  1215.4286375764655, 0.6480692446839038,
                             -0.06003913517996237, 0.134061  ], x - 0.377))

        conditions.append((x >= 0.381) & (x < 0.385))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4026.177810224078, 15.2332128956015,
                             0.0034859933811793005, 0.133909  ], x - 0.381))

        conditions.append((x >= 0.385) & (x < 0.388))
        eval_funcs.append(lambda x: np.polyval(
                          [  7238.622314002987, -33.08092082708747,
                             -0.06790483834476466, 0.133909  ], x - 0.385))

        conditions.append((x >= 0.388) & (x < 0.392))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3582.447447909331, 32.06667999893948,
                             -0.07094756082920863, 0.133603  ], x - 0.388))

        conditions.append((x >= 0.392) & (x < 0.396))
        eval_funcs.append(lambda x: np.polyval(
                          [ -496.1027599232692, -10.922689375972535,
                             0.013628401662659223, 0.133603  ], x - 0.392))

        conditions.append((x >= 0.396) & (x < 0.399))
        eval_funcs.append(lambda x: np.polyval(
                          [  5169.682959991873, -16.87592249505176,
                             -0.09756604582143809, 0.133451  ], x - 0.396))

        conditions.append((x >= 0.399) & (x < 0.403))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3710.297231720785, 29.65122414487512,
                             -0.05924014087196791, 0.133146  ], x - 0.399))

        conditions.append((x >= 0.403) & (x < 0.407))
        eval_funcs.append(lambda x: np.polyval(
                          [  1335.2490861665503, -14.872342635774345,
                             -0.00012461483556475847, 0.133146  ], x - 0.403))

        conditions.append((x >= 0.407) & (x < 0.41))
        eval_funcs.append(lambda x: np.polyval(
                          [  62.162287899734665, 1.1506463982240576,
                             -0.05501139978576523, 0.132993  ], x - 0.407))

        conditions.append((x >= 0.41) & (x < 0.414))
        eval_funcs.append(lambda x: np.polyval(
                          [  99.29447911487931, 1.7101069893216727,
                             -0.04642913962312804, 0.13284  ], x - 0.41))

        conditions.append((x >= 0.414) & (x < 0.418))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1367.1508902352487, 2.901640738700228,
                             -0.02798214871104043, 0.132688  ], x - 0.414))

        conditions.append((x >= 0.418) & (x < 0.422))
        eval_funcs.append(lambda x: np.polyval(
                          [  3009.9340818263636, -13.504169944122772,
                             -0.07039226553273063, 0.132535  ], x - 0.418))

        conditions.append((x >= 0.422) & (x < 0.425))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3766.2586617037464, 22.615039037793633,
                             -0.03394878915804718, 0.13223  ], x - 0.422))

        conditions.append((x >= 0.425) & (x < 0.429))
        eval_funcs.append(lambda x: np.polyval(
                          [  442.04340421567105, -11.281288917540113,
                             5.246120271340376e-05, 0.13223  ], x - 0.425))

        conditions.append((x >= 0.429) & (x < 0.433))
        eval_funcs.append(lambda x: np.polyval(
                          [  1024.177437691464, -5.976768066952054,
                             -0.06897976673525533, 0.132078  ], x - 0.429))

        conditions.append((x >= 0.433) & (x < 0.436))
        eval_funcs.append(lambda x: np.polyval(
                          [ -219.26177344735353, 6.313361185345528,
                             -0.06763339426168145, 0.131772  ], x - 0.433))

        conditions.append((x >= 0.436) & (x < 0.44))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1246.045366536695, 4.340005224319338,
                             -0.0356732950326868, 0.13162  ], x - 0.436))

        conditions.append((x >= 0.44) & (x < 0.444))
        eval_funcs.append(lambda x: np.polyval(
                          [  1685.2242205236093, -10.612539174121022,
                             -0.06076343083189353, 0.131467  ], x - 0.44))

        conditions.append((x >= 0.444) & (x < 0.447))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1673.052530751471, 9.61015147216231,
                             -0.06477298163972839, 0.131162  ], x - 0.444))

        conditions.append((x >= 0.447) & (x < 0.451))
        eval_funcs.append(lambda x: np.polyval(
                          [ -136.01397778448415, -5.447321304600952,
                             -0.05228449113704428, 0.131009  ], x - 0.447))

        conditions.append((x >= 0.451) & (x < 0.455))
        eval_funcs.append(lambda x: np.polyval(
                          [  1013.1055412231151, -7.079489038014774,
                             -0.1023917325075072, 0.130704  ], x - 0.451))

        conditions.append((x >= 0.455) & (x < 0.459))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1510.158187108603, 5.077777456662591,
                             -0.11039857883291583, 0.130246  ], x - 0.455))

        conditions.append((x >= 0.459) & (x < 0.462))
        eval_funcs.append(lambda x: np.polyval(
                          [  3192.1830955653672, -13.044120788640697,
                             -0.14226395216082816, 0.129789  ], x - 0.459))

        conditions.append((x >= 0.462) & (x < 0.466))
        eval_funcs.append(lambda x: np.polyval(
                          [ -290.7734358364457, 15.68552707144764,
                             -0.13433973331240734, 0.129331  ], x - 0.462))

        conditions.append((x >= 0.466) & (x < 0.47))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6388.896356541788, 12.196245841410299,
                             -0.022812641660975522, 0.129026  ], x - 0.466))

        conditions.append((x >= 0.47) & (x < 0.473))
        eval_funcs.append(lambda x: np.polyval(
                          [  13331.98866906993, -64.4705104370901,
                             -0.23190970004369224, 0.128721  ], x - 0.47))

        conditions.append((x >= 0.473) & (x < 0.477))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4862.530108551009, 55.51738758453929,
                             -0.25876906860134447, 0.127805  ], x - 0.473))

        conditions.append((x >= 0.477) & (x < 0.481))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3446.0432495141367, -2.8329737180728873,
                             -0.0480314131354786, 0.127347  ], x - 0.477))

        conditions.append((x >= 0.481) & (x < 0.484))
        eval_funcs.append(lambda x: np.polyval(
                          [  7073.528554829887, -44.185492712242606,
                             -0.2361052788567406, 0.126889  ], x - 0.481))

        conditions.append((x >= 0.484) & (x < 0.488))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2167.005810944752, 19.476264281226406,
                             -0.31023296414978924, 0.125974  ], x - 0.484))

        conditions.append((x >= 0.488) & (x < 0.492))
        eval_funcs.append(lambda x: np.polyval(
                          [  1096.89691411057, -6.527805450110662,
                             -0.2584391288253261, 0.124906  ], x - 0.488))

        conditions.append((x >= 0.492) & (x < 0.496))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4611.206845497518, 6.634957519216176,
                             -0.258010520548904, 0.123838  ], x - 0.492))

        conditions.append((x >= 0.496) & (x < 0.499))
        eval_funcs.append(lambda x: np.polyval(
                          [  7077.85513251692, -48.699524626754155,
                             -0.42626878897905585, 0.122617  ], x - 0.496))

        conditions.append((x >= 0.499) & (x < 0.503))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1790.0523813730345, 15.001171565898225,
                             -0.5273638481616238, 0.121091  ], x - 0.499))

        conditions.append((x >= 0.503) & (x < 0.507))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3315.948876083893, -6.479457010578303,
                             -0.49327698994034375, 0.119107  ], x - 0.503))

        conditions.append((x >= 0.507) & (x < 0.51))
        eval_funcs.append(lambda x: np.polyval(
                          [  3269.339553454129, -46.27084352358504,
                             -0.7042781920769975, 0.116818  ], x - 0.507))

        conditions.append((x >= 0.51) & (x < 0.514))
        eval_funcs.append(lambda x: np.polyval(
                          [  5219.889715327224, -16.84678754249792,
                             -0.8936310852752465, 0.114377  ], x - 0.51))

        conditions.append((x >= 0.514) & (x < 0.518))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8129.179805386951, 45.79188904142879,
                             -0.7778506792795229, 0.110867  ], x - 0.514))

        conditions.append((x >= 0.518) & (x < 0.521))
        eval_funcs.append(lambda x: np.polyval(
                          [  15924.926312182808, -51.75826862321478,
                             -0.8017161976066666, 0.107968  ], x - 0.518))

        conditions.append(( x >= 0.521) & (x < 0.525))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8857.592114293848, 91.5660681864303,
                             -0.6822927989170191, 0.105527  ], x - 0.521))

        conditions.append((x >= 0.525) & (x < 0.529))
        eval_funcs.append(lambda x: np.polyval(
                          [  5645.551478254127, -14.725037185096062,
                             -0.37492867491168164, 0.103696  ], x - 0.525))

        conditions.append((x >= 0.529) & (x < 0.533))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4162.113798722643, 53.021580553953505,
                             -0.22174250143625165, 0.102322  ], x - 0.529))

        conditions.append((x >= 0.533) & (x < 0.536))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1319.7028405038614, 3.076214969281744,
                             0.002648680656689534, 0.102017  ], x - 0.533))

        conditions.append((x >= 0.536) & (x < 0.54))
        eval_funcs.append(lambda x: np.polyval(
                          [  733.1530376397869, -8.80111059525302,
                             -0.014526006221224307, 0.102017  ], x - 0.536))

        conditions.append((x >= 0.54) & (x < 0.544))
        eval_funcs.append(lambda x: np.polyval(
                          [  719.1651094275551, -0.0032741435755676562,
                             -0.04974354517653869, 0.101865  ], x - 0.54))

        conditions.append((x >= 0.544) & (x < 0.547))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1181.1453817827548, 8.626707169555104,
                             -0.015249813072620522, 0.101712  ], x - 0.544))

        conditions.append((x >= 0.547) & (x < 0.551))
        eval_funcs.append(lambda x: np.polyval(
                          [  212.18127683644317, -2.0036012664897,
                             0.004619504636575709, 0.101712  ], x - 0.547))

        conditions.append((x >= 0.551) & (x < 0.555))
        eval_funcs.append(lambda x: np.polyval(
                          [ -59.10575093736702, 0.5425740555476208,
                             -0.0012246042071926116, 0.101712  ], x - 0.551))

        conditions.append((x >= 0.555) & (x < 0.558))
        eval_funcs.append(lambda x: np.polyval(
                          [  24.574741656401592, -0.16669495570078396,
                             0.0002789121921947376, 0.101712  ], x - 0.555))

        conditions.append((x >= 0.558) & (x < 0.562))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10.010709971262454, 0.054477719206830554,
                             -5.7739517287122876e-05, 0.101712  ], x - 0.558))

        conditions.append((x >= 0.562) & (x < 0.566))
        eval_funcs.append(lambda x: np.polyval(
                          [  22.814690252897822, -0.06565080044831899,
                             -0.00010243184225307665, 0.101712  ], x - 0.562))

        conditions.append((x >= 0.566) & (x < 0.57))
        eval_funcs.append(lambda x: np.polyval(
                          [ -81.24805104032566, 0.20812548258644747,
                             0.0004674668862994222, 0.101712  ], x - 0.566))

        conditions.append((x >= 0.57) & (x < 0.573))
        eval_funcs.append(lambda x: np.polyval(
                          [  451.99878807077937, -0.7668511298974615,
                             -0.0017674357029446352, 0.101712  ], x - 0.57))

        conditions.append((x >= 0.573) & (x < 0.577))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1189.9985404087415, 3.301137962739557,
                             0.005835424795581657, 0.101712  ], x - 0.573))

        conditions.append((x >= 0.577) & (x < 0.581))
        eval_funcs.append(lambda x: np.polyval(
                          [  1908.7987206739301, -10.97884452216535,
                             -0.024875401442121553, 0.101712  ], x - 0.577))

        conditions.append((x >= 0.581) & (x < 0.584))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1632.9334834077656, 11.926740125921834,
                             -0.021083819027095618, 0.101559  ], x - 0.581))

        conditions.append((x >= 0.584) & (x < 0.588))
        eval_funcs.append(lambda x: np.polyval(
                          [  293.20170141041064, -2.7696612247480683,
                             0.006387417676425704, 0.101559  ], x - 0.584))

        conditions.append((x >= 0.588) & (x < 0.592))
        eval_funcs.append(lambda x: np.polyval(
                          [ -81.17789467802031, 0.7487591921768625,
                             -0.0016961904538591254, 0.101559  ], x - 0.588))

        conditions.append((x >= 0.592) & (x < 0.595))
        eval_funcs.append(lambda x: np.polyval(
                          [  30.975832540816405, -0.22537554395938197,
                             0.0003973441390107983, 0.101559  ], x - 0.592))

        conditions.append((x >= 0.595) & (x < 0.599))
        eval_funcs.append(lambda x: np.polyval(
                          [ -5.941634343025867, 0.05340694890796599,
                             -0.00011856164614345012, 0.101559  ], x - 0.595))

        conditions.append((x >= 0.599) & (x < 0.603))
        eval_funcs.append(lambda x: np.polyval(
                          [  3.004697261146363, -0.017892663208344477,
                             2.349549665503608e-05, 0.101559  ], x - 0.599))

        conditions.append((x >= 0.603) & (x < 0.606))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8.785641255504606, 0.01816370392541191,
                             2.457965952330581e-05, 0.101559  ], x - 0.603))

        conditions.append((x >= 0.606) & (x < 0.61))
        eval_funcs.append(lambda x: np.polyval(
                          [  21.70491876996034, -0.06090706737412961,
                             -0.0001036504308228474, 0.101559  ], x - 0.606))

        conditions.append((x >= 0.61) & (x < 0.614))
        eval_funcs.append(lambda x: np.polyval(
                          [ -78.07106016273691, 0.1995519578653947,
                             0.00045092913114221344, 0.101559  ], x - 0.61))

        conditions.append((x >= 0.614) & (x < 0.618))
        eval_funcs.append(lambda x: np.polyval(
                          [  290.57932188098744, -0.7373007640874494,
                             -0.0017000660937460065, 0.101559  ], x - 0.614))

        conditions.append((x >= 0.618) & (x < 0.621))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1622.032059921667, 2.749651098484403,
                             0.006349335243841813, 0.101559  ], x - 0.618))

        conditions.append((x >= 0.621) & (x < 0.625))
        eval_funcs.append(lambda x: np.polyval(
                          [  1896.385846648717, -11.848637440810611,
                             -0.020947623783136842, 0.101559  ], x - 0.621))

        conditions.append((x >= 0.625) & (x < 0.629))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1182.6105128383015, 10.907992718974015,
                             -0.02471020267048324, 0.101407  ], x - 0.625))

        conditions.append((x >= 0.629) & (x < 0.632))
        eval_funcs.append(lambda x: np.polyval(
                          [  451.2850933540516, -3.2833334350856154,
                             0.005788434465070384, 0.101407  ], x - 0.629))

        conditions.append((x >= 0.632) & (x < 0.636))
        eval_funcs.append(lambda x: np.polyval(
                          [ -86.62881221996847, 0.778232405100852,
                             -0.001726868624883913, 0.101407  ], x - 0.632))

        conditions.append((x >= 0.636) & (x < 0.64))
        eval_funcs.append(lambda x: np.polyval(
                          [  44.02785854941666, -0.2613133415387705,
                             0.00034080762936441496, 0.101407  ], x - 0.636))

        conditions.append((x >= 0.64) & (x < 0.643))
        eval_funcs.append(lambda x: np.polyval(
                          [ -129.41122117654908, 0.26702096105423,
                             0.0003636381074262529, 0.101407  ], x - 0.64))

        conditions.append((x >= 0.643) & (x < 0.647))
        eval_funcs.append(lambda x: np.polyval(
                          [  319.9412010096276, -0.8976800295347127,
                             -0.0015283390980151965, 0.101407  ], x - 0.643))

        conditions.append((x >= 0.647) & (x < 0.651))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1150.865990280782, 2.9416143825808216,
                             0.006647398314169247, 0.101407  ], x - 0.647))

        conditions.append((x >= 0.651) & (x < 0.655))
        eval_funcs.append(lambda x: np.polyval(
                          [  1892.8977601135002, -10.868777500788571,
                             -0.02506125415866179, 0.101407  ], x - 0.651))

        conditions.append((x >= 0.655) & (x < 0.658))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1598.4005757997904, 11.845995620573456,
                             -0.02115238167952226, 0.101254  ], x - 0.655))

        conditions.append((x >= 0.658) & (x < 0.662))
        eval_funcs.append(lambda x: np.polyval(
                          [  211.97885932341057, -2.539609561624671,
                             0.006766776497324119, 0.101254  ], x - 0.658))

        conditions.append((x >= 0.662) & (x < 0.666))
        eval_funcs.append(lambda x: np.polyval(
                          [  209.9104841952816, 0.004136750256258037,
                             -0.003375114748149544, 0.101254  ], x - 0.662))

        conditions.append((x >= 0.666) & (x < 0.669))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1589.2077974525507, 2.5230625605996395,
                             0.006733682495274056, 0.101254  ], x - 0.666))

        conditions.append((x >= 0.669) & (x < 0.673))
        eval_funcs.append(lambda x: np.polyval(
                          [  1884.7364461400336, -11.779807616473331,
                             -0.021036552672347036, 0.101254  ], x - 0.669))

        conditions.append((x >= 0.673) & (x < 0.677))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1158.7784224635257, 10.837029737207097,
                             -0.02480766418941198, 0.101102  ], x - 0.673))

        conditions.append((x >= 0.677) & (x < 0.68))
        eval_funcs.append(lambda x: np.polyval(
                          [  326.4138407855707, -3.0683113323552234,
                             0.006267209429995537, 0.101102  ], x - 0.677))

        conditions.append((x >= 0.68) & (x < 0.684))
        eval_funcs.append(lambda x: np.polyval(
                          [  240.73949525410794, -0.13058676528508462,
                             -0.003329484862925395, 0.101102  ], x - 0.68))

        conditions.append((x >= 0.684) & (x < 0.688))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1138.4040936280414, 2.7582871777642133,
                             0.007181316786991129, 0.101102  ], x - 0.684))

        conditions.append((x >= 0.688) & (x < 0.692))
        eval_funcs.append(lambda x: np.polyval(
                          [  1922.2518792579017, -10.902561945771922,
                             -0.025395782285038814, 0.101102  ], x - 0.688))

        conditions.append((x >= 0.692) & (x < 0.695))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1793.9104632371034, 12.164460605322915,
                             -0.02034818764683482, 0.100949  ], x - 0.692))

        conditions.append((x >= 0.695) & (x < 0.699))
        eval_funcs.append(lambda x: np.polyval(
                          [  732.4962985964537, -3.9807335638110315,
                             0.004202993477700858, 0.100949  ], x - 0.695))

        conditions.append((x >= 0.699) & (x < 0.703))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1672.1147110767547, 4.8092220193464215,
                             0.007516947299842421, 0.100949  ], x - 0.699))

        conditions.append((x >= 0.703) & (x < 0.706))
        eval_funcs.append(lambda x: np.polyval(
                          [  3263.61995012534, -15.256154513574645,
                             -0.034270782677070545, 0.100949  ], x - 0.703))

        conditions.append((x >= 0.706) & (x < 0.71))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1173.4830653174774, 14.116425037553444,
                             -0.037689971105134154, 0.100797  ], x - 0.706))

        conditions.append((x >= 0.71) & (x < 0.714))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1190.7971921893286, 0.03462825374370304,
                             0.01891424206005448, 0.100797  ], x - 0.71))

        conditions.append((x >= 0.714) & (x < 0.717))
        eval_funcs.append(lambda x: np.polyval(
                          [  3303.534588074269, -14.254938052528244,
                             -0.03796699713508377, 0.100797  ], x - 0.714))

        conditions.append((x >= 0.717) & (x < 0.721))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1725.393836769556, 15.476873240140197,
                             -0.034301191572247904, 0.100644  ], x - 0.717))

        conditions.append((x >= 0.721) & (x < 0.725))
        eval_funcs.append(lambda x: np.polyval(
                          [  888.5325637776872, -5.227852801094491,
                             0.00669489018393496, 0.100644  ], x - 0.721))

        conditions.append((x >= 0.725) & (x < 0.729))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1828.736418341193, 5.434537964237766,
                             0.0075216308365080595, 0.100644  ], x - 0.725))

        conditions.append((x >= 0.729) & (x < 0.732))
        eval_funcs.append(lambda x: np.polyval(
                          [  3923.5900775040905, -16.510299055856567,
                             -0.0367814135299672, 0.100644  ], x - 0.729))

        conditions.append((x >= 0.732) & (x < 0.736))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2831.360674639062, 18.80201164168027,
                             -0.029906275772496077, 0.100491  ], x - 0.732))

        conditions.append((x >= 0.736) & (x < 0.74))
        eval_funcs.append(lambda x: np.polyval(
                          [  2380.797552355198, -15.1743164539885,
                             -0.015395495021728989, 0.100491  ], x - 0.736))

        conditions.append((x >= 0.74) & (x < 0.743))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1963.7798202482595, 13.395254174273905,
                             -0.022511744140587386, 0.100339  ], x - 0.74))

        conditions.append((x >= 0.743) & (x < 0.747))
        eval_funcs.append(lambda x: np.polyval(
                          [  767.333192093047, -4.278764207960444,
                             0.004837725758353019, 0.100339  ], x - 0.743))

        conditions.append((x >= 0.747) & (x < 0.751))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1697.2838564850151, 4.929234097156129,
                             0.007439605315135756, 0.100339  ], x - 0.747))

        conditions.append((x >= 0.751) & (x < 0.754))
        eval_funcs.append(lambda x: np.polyval(
                          [  3323.40706232091, -15.438172180664077,
                             -0.034596147018896045, 0.100339  ], x - 0.751))

        conditions.append((x >= 0.754) & (x < 0.758))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1274.7985062925447, 14.472491380224142,
                             -0.03749318942021586, 0.100186  ], x - 0.754))

        conditions.append((x >= 0.758) & (x < 0.762))
        eval_funcs.append(lambda x: np.polyval(
                          [ -862.2531586493396, -0.8250906952864105,
                             0.017096413319535105, 0.100186  ], x - 0.758))

        conditions.append((x >= 0.762) & (x < 0.766))
        eval_funcs.append(lambda x: np.polyval(
                          [  2348.81114088992, -11.172128599078496,
                             -0.030892463857924563, 0.100186  ], x - 0.762))

        conditions.append((x >= 0.766) & (x < 0.769))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4834.917487440602, 17.013605091600567,
                             -0.007526557887836259, 0.100034  ], x - 0.766))

        conditions.append((x >= 0.769) & (x < 0.773))
        eval_funcs.append(lambda x: np.polyval(
                          [  4093.144292536792, -26.50065229536489,
                             -0.035987699499129236, 0.100034  ], x - 0.769))

        conditions.append((x >= 0.773) & (x < 0.777))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2434.145315001528, 22.617079215076657,
                             -0.0515219918202822, 0.099728  ], x - 0.773))

        conditions.append((x >= 0.777) & (x < 0.78))
        eval_funcs.append(lambda x: np.polyval(
                          [  800.2585460630472, -6.5926645649417015,
                             0.012575666780257686, 0.099728  ], x - 0.777))

        conditions.append((x >= 0.78) & (x < 0.784))
        eval_funcs.append(lambda x: np.polyval(
                          [  183.41815419920795, 0.6096623496257284,
                             -0.005373339865690248, 0.099728  ], x - 0.78))

        conditions.append((x >= 0.784) & (x < 0.788))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1221.921945808904, 2.8106802000162268,
                             0.00830803033287758, 0.099728  ], x - 0.784))

        conditions.append((x >= 0.788) & (x < 0.791))
        eval_funcs.append(lambda x: np.polyval(
                          [  1416.5849164695105, -11.852383149690638,
                             -0.02785878146582008, 0.099728  ], x - 0.788))

        conditions.append((x >= 0.791) & (x < 0.795))
        eval_funcs.append(lambda x: np.polyval(
                          [  1180.4852015717026, 0.8968810985349469,
                             -0.060725287619287124, 0.099576  ], x - 0.791))

        conditions.append((x >= 0.795) & (x < 0.799))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3960.2415571259885, 15.062703517395395,
                             0.003113050844434294, 0.099423  ], x - 0.795))

        conditions.append((x >= 0.799) & (x < 0.802))
        eval_funcs.append(lambda x: np.polyval(
                          [  6910.092732903701, -32.46019516811653,
                             -0.06647691575845022, 0.099423  ], x - 0.799))

        conditions.append((x >= 0.802) & (x < 0.806))
        eval_funcs.append(lambda x: np.polyval(
                          [ -2766.060920832377, 29.73063942801683,
                             -0.07466558297874931, 0.099118  ], x - 0.802))

        conditions.append((x >= 0.806) & (x < 0.81))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1035.015109846516, -3.462091621971725,
                             0.0304086082454312, 0.099118  ], x - 0.806))

        conditions.append((x >= 0.81) & (x < 0.814))
        eval_funcs.append(lambda x: np.polyval(
                          [  2140.49636021833, -15.882272940129912,
                             -0.04696885000297549, 0.099118  ], x - 0.81))

        conditions.append((x >= 0.814) & (x < 0.817))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1014.2046571034608, 9.803683382489368,
                             -0.07128320823353704, 0.098813  ], x - 0.814))

        conditions.append((x >= 0.817) & (x < 0.821))
        eval_funcs.append(lambda x: np.polyval(
                          [ -53.67076211490019, 0.6758414685582176,
                             -0.03984463368039427, 0.09866  ], x - 0.817))

        conditions.append((x >= 0.821) & (x < 0.825))
        eval_funcs.append(lambda x: np.polyval(
                          [ -85.19192370462227, 0.031792323179412996,
                             -0.03701409851344375, 0.098508  ], x - 0.821))

        conditions.append((x >= 0.825) & (x < 0.828))
        eval_funcs.append(lambda x: np.polyval(
                          [ -797.721716704623, -0.9905107612760462,
                             -0.04084897226583031, 0.098355  ], x - 0.825))

        conditions.append((x >= 0.828) & (x < 0.832))
        eval_funcs.append(lambda x: np.polyval(
                          [  1547.5342519363887, -8.170006211617661,
                             -0.06833052318451145, 0.098202  ], x - 0.828))

        conditions.append((x >= 0.832) & (x < 0.836))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1262.0431538731202, 10.400404811619028,
                             -0.059408928784505995, 0.097897  ], x - 0.832))

        conditions.append((x >= 0.836) & (x < 0.839))
        eval_funcs.append(lambda x: np.polyval(
                          [  1.7889757820873253, -4.744113034858421,
                             -0.03678376167746358, 0.097745  ], x - 0.836))

        conditions.append((x >= 0.839) & (x < 0.843))
        eval_funcs.append(lambda x: np.polyval(
                          [  491.3866594860325, -4.72801225281965,
                             -0.06520013754049778, 0.097592  ], x - 0.839,))

        conditions.append((x >= 0.843) & (x < 0.847))
        eval_funcs.append(lambda x: np.polyval(
                          [ -92.92717102033794, 1.168627661012742,
                             -0.07943767590772542, 0.097287  ], x - 0.843,))

        conditions.append((x >= 0.847) & (x < 0.851))
        eval_funcs.append(lambda x: np.polyval(
                          [ -119.67797540468248, 0.05350160876868909,
                             -0.0745491588285997, 0.096982  ], x - 0.847))

        conditions.append((x >= 0.851) & (x < 0.854))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1998.4898815402985, -1.382634096087492,
                             -0.07986568877787495, 0.096677  ], x - 0.851))

        conditions.append((x >= 0.854) & (x < 0.858))
        eval_funcs.append(lambda x: np.polyval(
                          [  4193.555767236822, -19.369043029950216,
                             -0.14212072015598806, 0.096371  ], x - 0.854))

        conditions.append((x >= 0.858) & (x < 0.862))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6517.632321209027, 30.95362617689169,
                             -0.0957823875682221, 0.095761  ], x - 0.858,))

        conditions.append((x >= 0.862) & (x < 0.865))
        eval_funcs.append(lambda x: np.polyval(
                          [  11048.920141182143, -47.257961677616684,
                             -0.16099972957112224, 0.095456  ], x - 0.862,))

        conditions.append((x >= 0.865) & (x < 0.869))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6297.038909199155, 52.182319593022655,
                             -0.1462266558249042, 0.094846  ], x - 0.865))

        conditions.append((x >= 0.869) & (x < 0.873))
        eval_funcs.append(lambda x: np.polyval(
                          [  628.4097494844855, -23.38214731736727,
                             -0.031025966722282552, 0.094693  ], x - 0.869,))

        conditions.append((x >= 0.873) & (x < 0.876))
        eval_funcs.append(lambda x: np.polyval(
                          [  3567.7594359214195, -15.841230323553413,
                             -0.1879194772859655, 0.094235  ], x - 0.873))

        conditions.append((x >= 0.876) & (x < 0.88))
        eval_funcs.append(lambda x: np.polyval(
                          [ -1949.1914963468762, 16.2686045997394,
                             -0.18663735445740756, 0.093625  ], x - 0.876))

        conditions.append((x >= 0.88) & (x < 0.884))
        eval_funcs.append(lambda x: np.polyval(
                          [  1627.2801818646924, -7.1216933564230915,
                             -0.15004970948414245, 0.093014  ], x - 0.88))

        conditions.append((x >= 0.884) & (x < 0.888))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4575.554231111926, 12.405668825953253,
                             -0.12891380760602184, 0.092404  ], x - 0.884))

        conditions.append((x >= 0.888) & (x < 0.891))
        eval_funcs.append(lambda x: np.polyval(
                          [  7940.519177845019, -42.50098194738993,
                             -0.24929506009176852, 0.091794  ], x - 0.888))

        conditions.append((x >= 0.891) & (x < 0.895))
        eval_funcs.append(lambda x: np.polyval(
                          [ -3418.614289910524, 28.963690653215295,
                             -0.2899069339742924, 0.090878  ], x - 0.891))

        conditions.append((x >= 0.895) & (x < 0.899))
        eval_funcs.append(lambda x: np.polyval(
                          [  2595.601122944965, -12.05968082571102,
                             -0.22229089466427523, 0.089963  ], x - 0.895))

        conditions.append((x >= 0.899) & (x < 0.902))
        eval_funcs.append(lambda x: np.polyval(
                          [ -7379.604879290401, 19.087532649628557,
                             -0.19417948736860494, 0.089047  ], x - 0.899))

        conditions.append((x >= 0.902) & (x < 0.906))
        eval_funcs.append(lambda x: np.polyval(
                          [  7794.954266725975, -47.32891126398517,
                             -0.2789036232116747, 0.088437  ], x - 0.902))

        conditions.append((x >= 0.906) & (x < 0.91))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10529.065701637304, 46.21053993672657,
                             -0.2833771085207089, 0.087063  ], x - 0.906))

        conditions.append((x >= 0.91) & (x < 0.913))
        eval_funcs.append(lambda x: np.polyval(
                          [  16759.55794306497, -80.13824848292111,
                             -0.4190879427054875, 0.085995  ], x - 0.91))

        conditions.append((x >= 0.913) & (x < 0.917))
        eval_funcs.append(lambda x: np.polyval(
                          [ -8773.857679899693, 70.69777300466382,
                             -0.4474093691402596, 0.084469  ], x - 0.913))

        conditions.append((x >= 0.917) & (x < 0.921))
        eval_funcs.append(lambda x: np.polyval(
                          [  3739.1518971666014, -34.58851915413262,
                             -0.3029723537381347, 0.083249  ], x - 0.917))

        conditions.append((x >= 0.921) & (x < 0.925))
        eval_funcs.append(lambda x: np.polyval(
                          [ -6167.124908766686, 10.281303611866596,
                             -0.4002012159071987, 0.081723  ], x - 0.921))

        conditions.append((x >= 0.925) & (x < 0.928))
        eval_funcs.append(lambda x: np.polyval(
                          [  15979.115019970644, -63.72419529333374,
                             -0.6139727826330673, 0.079892  ], x - 0.925))

        conditions.append((x >= 0.928) & (x < 0.932))
        eval_funcs.append(lambda x: np.polyval(
                          [ -10951.219418234163, 80.08783988640221,
                             -0.5648818488538618, 0.077908  ], x - 0.928))

        conditions.append((x >= 0.932) & (x < 0.936))
        eval_funcs.append(lambda x: np.polyval(
                          [  2805.9271479697836, -51.326793132407836,
                             -0.44983766183788426, 0.076229  ], x - 0.932))

        conditions.append((x >= 0.936) & (x < 0.939))
        eval_funcs.append(lambda x: np.polyval(
                          [  1748.278429430659, -17.655667356770387,
                             -0.7257675037945975, 0.073788  ], x - 0.936))

        conditions.append((x >= 0.939) & (x < 0.943))
        eval_funcs.append(lambda x: np.polyval(
                          [  4214.5397692609, -1.9211614918954432,
                             -0.7844979903405918, 0.071499  ], x - 0.939))

        conditions.append((x >= 0.943) & (x < 0.947))
        eval_funcs.append(lambda x: np.polyval(
                          [ -15346.493100356824, 48.653315739235495,
                             -0.5975693733512317, 0.0686  ], x - 0.943))

        conditions.append((x >= 0.947) & (x < 0.95))
        eval_funcs.append(lambda x: np.polyval(
                          [  25832.035627735237, -135.50460146504656,
                             -0.9449745162544761, 0.066006  ], x - 0.947))

        conditions.append((x >= 0.95) & (x < 0.954))
        eval_funcs.append(lambda x: np.polyval(
                          [ -12806.107102648655, 96.98371918457075,
                             -1.0605371630959035, 0.062649  ], x - 0.95,))

        conditions.append((x >= 0.954) & (x < 0.958))
        eval_funcs.append(lambda x: np.polyval(
                          [  8398.050920957965, -56.689566047213205,
                             -0.8993605505464733, 0.059139  ], x - 0.954))

        conditions.append((x >= 0.958) & (x < 0.962))
        eval_funcs.append(lambda x: np.polyval(
                          [ -13645.471581183196, 44.087045004282274,
                             -0.9497706347181964, 0.055172  ], x - 0.958))

        conditions.append((x >= 0.962) & (x < 0.965))
        eval_funcs.append(lambda x: np.polyval(
                          [  20744.379906349895, -119.6586139699161,
                             -1.2520569105807324, 0.051205  ], x - 0.962))

        conditions.append((x >= 0.965) & (x < 0.969))
        eval_funcs.append(lambda x: np.polyval(
                          [ -4937.6802387593925, 67.04080518723337,
                             -1.4099103369287818, 0.046932  ], x - 0.965))

        conditions.append((x >= 0.969) & (x < 0.973))
        eval_funcs.append(lambda x: np.polyval(
                          [ -11207.001399819874, 7.788642322120893,
                             -1.1105925468913653, 0.042049  ], x - 0.969))

        conditions.append((x >= 0.973) & (x < 0.976))
        eval_funcs.append(lambda x: np.polyval(
                          [  20663.585066619387, -126.69537447571766,
                             -1.5862194755057528, 0.037014  ], x - 0.973))

        conditions.append((x >= 0.976) & (x <= 0.98))
        eval_funcs.append(lambda x: np.polyval(
                          [  20663.585066619267, 59.276891123856856,
                             -1.788474925561335, 0.031673  ], x - 0.976))

        return np.piecewise(x, conditions, eval_funcs)
