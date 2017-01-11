"""Module exchange_model contains code for setting up different three-site
exchange models."""

from scipy import constants as cst

from chemex import parameters


def update_params(params=None,
                  map_names=None,
                  model=None,
                  temperature=None,
                  p_total=None,
                  l_total=None):
    """Update the experimental and fitting parameters depending on the
    model."""

    if model not in {'3st.pb_kex', '3st.temperature'}:
        print("Warning: The 'model' option should be either '3st.pb_kex'")
        print("or '3st.temperature'.")
        print("\nSet it to the default model: '3st.pb_kex'.")
        model = '3st.pb_kex'

    if model == '3st.pb_kex':
        pass

    elif model == '3st.temperature':
        # temperature-dependent population and kinetics:
        #
        # G ^
        #   |         -TSAB-
        #   |        /      \       -TSBC-
        #   |        |      |      /      \      -TSAC-
        #   |        |      \     /       \     /      \
        #   |       /        --B--         --C--       \
        #   |  --A--                                    --A--
        #   |
        #
        # The exchange rate is calculated relative to that at a reference temperature, T0,
        # using transition state theory, assuming (a) no temperature dependence of the
        # pre-exponential factor, and (b) no significant effect of changes in heat capacity over
        # the observed temperature range. All energies are defined relative to state A.
        #
        # pB(T1) = KAB(T1) / (1 + KAB(T1) + KAC(T1))
        # pC(T1) = KAC(T1) / (1 + KAB(T1) + KAC(T1))
        #     where KAB(T1) = ((pB(T0) / (1-pB(T0)-pC(T0))) * exp(-dH_AB/R * (1/T1 - 1/T0)))
        #       and KAC(T1) = ((pC(T0) / (1-pB(T0)-pC(T0))) * exp(-dH_AC/R * (1/T1 - 1/T0)))
        #
        # kex_AB(T1) = kAB(T1) + kBA(T1)
        #         = kex_AB(T0) * exp(-dH_TSAB/R * (1/T1 - 1/T0)) *
        #           ( (pB(T0)) / (1-pC(T0)) +
        #             (1-pB(T0)-pC(T0)) / (1-pC(T0)) * exp(dH_AB/R * (1/T1 - 1/T0)) )
        #
        # kex_AC(T1) = kAC(T1) + kCA(T1)
        #         = kex_AC(T0) * exp(-dH_TSAC/R * (1/T1 - 1/T0)) *
        #           ( (pC(T0)) / (1-pB(T0)) +
        #             (1-pB(T0)-pC(T0)) / (1-pB(T0)) * exp(dH_AC/R * (1/T1 - 1/T0)) )
        #
        # kex_BC(T1) = kBC(T1) + kCB(T1)
        #         = kex_BC(T0) * exp(-dH_TSBC/R * (1/T1 - 1/T0)) *
        #           ( (pC(T0)) / (pB(T0)+pC(T0)) * exp(dH_AB/R * (1/T1 - 1/T0)) +
        #             (pB(T0)) / (pB(T0)+pC(T0)) * exp(dH_AC/R * (1/T1 - 1/T0)) )
        #
        # Parameters:
        #     T0: reference temperature (in degrees C)
        #     pb0: population B at reference temperature
        #     pc0: population C at reference temperature
        #     kexAB0: exchange rate AB at reference temperature (in s-1)
        #     kexBC0: exchange rate BC at reference temperature (in s-1)
        #     kexAC0: exchange rate AC at reference temperature (in s-1)
        #     dH_AB: enthalpy of reaction, deltaH_(A-B) = H_B - H_A (in kJ mol-1)
        #     dH_AC: enthalpy of reaction, deltaH_(A-C) = H_C - H_A (in kJ mol-1)
        #     dH_TSAB: enthalpy of activation, deltaH_(A-TSAB) = H_TSAB - H_A (in kJ mol-1)
        #     dH_TSBC: enthalpy of activation, deltaH_(A-TSBC) = H_TSBC - H_A (in kJ mol-1)
        #     dH_TSAC: enthalpy of activation, deltaH_(A-TSAC) = H_TSAC - H_A (in kJ mol-1)

        map_names.update({
            'T0': parameters.ParameterName('T0').to_full_name(),
            'pb0': parameters.ParameterName('pb0').to_full_name(),
            'pc0': parameters.ParameterName('pc0').to_full_name(),
            'kexAB0': parameters.ParameterName('kexAB0').to_full_name(),
            'kexBC0': parameters.ParameterName('kexBC0').to_full_name(),
            'kexAC0': parameters.ParameterName('kexAC0').to_full_name(),
            'dH_AB': parameters.ParameterName('dH_AB').to_full_name(),
            'dH_AC': parameters.ParameterName('dH_AC').to_full_name(),
            'dH_TSAB': parameters.ParameterName('dH_TSAB').to_full_name(),
            'dH_TSBC': parameters.ParameterName('dH_TSBC').to_full_name(),
            'dH_TSAC': parameters.ParameterName('dH_TSAC').to_full_name(),
        })

        R = cst.R * 0.001 # convert to kJ mol-1 K-1

        # KAB= (({pb0} / (1-{pb0}-{pc0})) * exp(-{dH_AB}/{R} * (1/({T}+273.15) - 1/({T0}+273.15))))
        # KAC= (({pc0} / (1-{pb0}-{pc0})) * exp(-{dH_AC}/{R} * (1/({T}+273.15) - 1/({T0}+273.15))))
        pb = ('(({pb0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AB}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15))))'
            ' / (1.0 + '
            '(({pb0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AB}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15)))) +'
            '(({pc0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AC}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15)))) )'
            .format(R=R, T=temperature, **map_names))

        pc = ('(({pc0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AC}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15))))'
            ' / (1.0 + '
            '(({pb0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AB}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15)))) +'
            '(({pc0} / (1.0-{pb0}-{pc0})) * exp(-{dH_AC}/{R} * (1.0/({T}+273.15) - 1.0/({T0}+273.15)))) )'
            .format(R=R, T=temperature, **map_names))

        kex_ab = ('{kexAB0} * exp(-{dH_TSAB}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) * '
            '( ({pb0}) / (1.0-{pc0}) +'
            '(1.0-{pb0}-{pc0}) / (1.0-{pc0}) * exp({dH_AB}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) )'
            .format(R=R, T=temperature, **map_names))

        kex_bc = ('{kexBC0} * exp(-{dH_TSBC}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) * '
            '( ({pc0}) / ({pb0}+{pc0}) * exp({dH_AB}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) + '
            '({pb0}) / ({pb0}+{pc0}) * exp({dH_AC}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) )'
            .format(R=R, T=temperature, **map_names))

        kex_ac = ('{kexAC0} * exp(-{dH_TSAC}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) * '
            '( ({pc0}) / (1.0-{pb0}) +'
            '(1.0-{pb0}-{pc0}) / (1.0-{pb0}) * exp({dH_AC}/{R} * (1.0/({T}+273.15) - 1/({T0}+273.15))) )'
            .format(R=R, T=temperature, **map_names))

        params.add_many(  # Name, Value, Vary, Min, Max, Expr
            (map_names['T0'], 25, False, None, None, None),
            (map_names['pb0'], 0.05, True, 0.0, 1.0, None),
            (map_names['pc0'], 0.05, True, 0.0, 1.0, None),
            (map_names['kexAB0'], 100.0, True, 0.0, None, None),
            (map_names['kexBC0'], 100.0, True, 0.0, None, None),
            (map_names['kexAC0'], 0.0, False, 0.0, None, None),
            (map_names['dH_AB'], 0.0, True, -200, 200, None),
            (map_names['dH_AC'], 0.0, True, -200, 200, None),
            (map_names['dH_TSAB'], 10.0, True, 0.0, 200, None),
            (map_names['dH_TSBC'], 10.0, True, 0.0, 200, None),
            (map_names['dH_TSAC'], 10.0, False, 0.0, 200, None),
            (map_names['pb'], 0.0, None, 0.0, 1.0, pb),
            (map_names['pc'], 0.0, None, 0.0, 1.0, pc),
            (map_names['kex_ab'], 0.0, None, 0.0, None, kex_ab),
            (map_names['kex_bc'], 0.0, None, 0.0, None, kex_bc),
            (map_names['kex_ac'], 0.0, None, 0.0, None, kex_ac),
            )

    return map_names, params
