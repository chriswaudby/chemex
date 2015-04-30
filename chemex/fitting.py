from __future__ import print_function

import copy
import itertools

import lmfit as lf

from chemex import parameters
from chemex import utils
from chemex import chi2


product = itertools.product


def run_fit(fit_filename, params, data):

    utils.header1("Fit")

    fit_config = utils.read_cfg_file(fit_filename)

    if not fit_config.sections():
        fit_config.add_section('Standard Calculation')

    for section in fit_config.sections():

        utils.header2(section)

        items = fit_config.items(section)

        parameters.set_param_status(params, items)

        independent_clusters = find_independent_clusters(data, params)
        independent_clusters_no = len(independent_clusters)

        if independent_clusters_no > 1:

            for index, independent_cluster in enumerate(independent_clusters, 1):

                print('\nChi2 / Reduced Chi2 (cluster {}/{}):'.format(index, independent_clusters_no))

                c_data, c_params = independent_cluster

                func = chi2.make_calc_residuals(verbose=True)
                args = (c_data, )

                lf.minimize(func, c_params, args=args)

                for name, param in c_params.items():
                    params[name] = param

            print("\nFinalizing...")

        else:
            print("\nChi2 / Reduced Chi2:")

        func = chi2.make_calc_residuals(verbose=True)
        args = (data, )

        out = lf.minimize(func, params, args=args)

        print("\nFinal Chi2        : {:.3e}".format(out.chisqr))
        print("Final Reduced Chi2: {:.3e}".format(out.redchi))

    return params


def find_independent_clusters(data, params):
    """
    Finds clusters of data points that depend on independent sets of variables.
    For example, if the population of the minor state and the exchange rate are
    set to 'fix', chances are that the fit can be decomposed
    residue-specifically.
    """

    clusters = []

    for profile in data:

        params_profile = lf.Parameters(
            params[name] for name in profile.map_names.values()
        )

        for data_cluster, params_cluster in clusters:

            params_profile_fit = set([
                name for name, param in params_profile if param.vary
            ])

            params_cluster_fit = set([
                name for name, param in params_cluster if param.vary
            ])

            if params_profile_fit.intersection(params_cluster_fit):
                data_cluster.append(profile)
                params_cluster.update(params_profile)
                break

        else:
            data_cluster = [profile]
            params_cluster = copy.deepcopy(params_profile)
            clusters.append((data_cluster, params_cluster))

    return clusters

