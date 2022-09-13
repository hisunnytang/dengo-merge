import logging
import sys

import numpy as np
import sympy
import yt
from sympy import lambdify

import dengo.primordial_cooling
import dengo.primordial_rates
from dengo.chemical_network import (
    ChemicalNetwork,
    cooling_registry,
    reaction_registry,
    species_registry,
)

grackle_species_field = [
    "HI_Density",
    "HII_Density",
    "HeI_Density",
    "HeII_Density",
    "HeIII_Density",
    "Electron_Density",
    "HM_Density",
    "H2I_Density",
    "H2II_Density",
]
dengo_species_field = [
    "H_1Density",
    "H_2Density",
    "He_1Density",
    "He_2Density",
    "He_3Density",
    "deDensity",
    "H_m0Density",
    "H2_1Density",
    "H2_2Density",
]

field2mass = {
    "H_1Density": 1.00794,
    "H_2Density": 1.00794,
    "H_m0Density": 1.00894,
    "H2_1Density": 2.01588,
    "H2_2Density": 2.01588,
    "He_1Density": 4.0022602,
    "He_2Density": 4.0022602,
    "He_3Density": 4.0022602,
    "deDensity": 1.00794,
}


dengo.primordial_rates.setup_primordial()

from dengo.chemistry_constants import tevk
from dengo.reaction_classes import ReactionCoefficient


class temperature_data:
    def __init__(self, T):
        self.T = T
        self.logT = np.log(self.T)
        self.tev = self.T / tevk
        self.logtev = np.log(self.tev)
        self.z = 16.5
        self.threebody = 4


def evaluate_cooling_terms(cooling_action, init_values):
    eval_rates = {}
    _names = {}
    init_values["z"] = np.ones_like(init_values["T"]) * 20.0
    for tab_name in cooling_action.tables:
        symbol_name = f"{cooling_action.name}_{tab_name}[i]"
        _names[symbol_name] = f"{cooling_action.name}_{tab_name}"
        eval_rates[f"{cooling_action.name}_{tab_name}"] = cooling_action.tables[
            tab_name
        ](temperature_data(init_values["T"]))

    sym_names = sorted(list(init_values.keys()) + list(eval_rates.keys()))
    sym_names += [
        "h2_optical_depth_approx",
    ]
    symbols = sympy.symbols(sym_names)

    eq_ = cooling_action.equation
    for i, j in _names.items():

        tmp = ReactionCoefficient(i)
        eq_ = eq_.subs(tmp, j)
    inputs = []
    for s in sym_names:
        if s in init_values:
            inputs.append(init_values[s])
        elif s in eval_rates:
            inputs.append(eval_rates[s])
        else:
            inputs.append(np.ones_like(inputs[-1]))
    func1 = lambdify(symbols, eq_, "numpy")
    return func1(*inputs)


def prepare_dengo_inputs(data):
    inits = {}
    inits["H2_1"] = (data["H2I_Density"] / yt.units.mh / 2.01588).to("1/cm**3").v
    inits["H2_2"] = (data["H2II_Density"] / yt.units.mh / 2.01588).to("1/cm**3").v
    inits["H_1"] = (data["HI_Density"] / yt.units.mh / 1.00794).to("1/cm**3").v
    inits["H_2"] = (data["HII_Density"] / yt.units.mh / 1.00794).to("1/cm**3").v
    inits["H_m0"] = (data["HM_Density"] / yt.units.mh / 1.00794).to("1/cm**3").v
    inits["He_1"] = (data["HeI_Density"] / yt.units.mh / 4).to("1/cm**3").v
    inits["He_2"] = (data["HeII_Density"] / yt.units.mh / 4).to("1/cm**3").v
    inits["He_3"] = (data["HeIII_Density"] / yt.units.mh / 4).to("1/cm**3").v
    inits["de"] = (data["Electron_Density"] / yt.units.mh / 1.00794).to("1/cm**3").v
    inits["ge"] = data["thermal_energy"]
    inits["mdensity"] = data["density"].to("g/cm**3").v
    inits["T"] = data["temperature"].v
    return inits


def add_dengo_thermal_fields(thermal_process_name):
    def add_terms(field, data):
        inits = prepare_dengo_inputs(data)
        tdyn = data["dynamical_time"].to("s").v
        thermal_process = (
            evaluate_cooling_terms(cooling_registry[thermal_process_name], inits)
            * inits["mdensity"] ** (-1)
            * yt.units.erg
            / yt.units.g
            * tdyn
        )
        if thermal_process_name == "gloverabel08":
            h2optical = np.minimum(
                1.0, pow((0.76 * inits["mdensity"] / (1.34e-14)), -0.45)
            )
            return thermal_process * h2optical
        if thermal_process_name == "cie_cooling":
            mdensity = inits["mdensity"]
            tau = (mdensity / 3.3e-8) ** 2.8
            tau = np.maximum(tau, 1e-4)

            cie_fudge = np.minimum(1, (1 - np.exp(-tau)) / tau)
            tau = (mdensity / 3.3e-6) ** 8
            tau = np.maximum(tau, 1e-5)
            cie_fudge *= np.minimum(1, (1 - np.exp(-tau)) / tau)
            return thermal_process * cie_fudge
        return thermal_process

    return add_terms


def calculate_temperature_from_dict(data):
    species = ["H_1", "H_2", "H_m0", "He_1", "He_2", "He_3"]
    species_h2 = ["H2_1", "H2_2"]

    epsilon = data["ge"]

    tmp = 0
    gamma = 5 / 3
    for i, s in enumerate(species):
        tmp += data[s] * (1 / (gamma - 1)) / yt.units.atomic_mass_unit_cgs.v
    gammaH2 = 7.0 / 5.0
    for s in species_h2:
        tmp += data[s] * (1 / (gammaH2 - 1)) / yt.units.atomic_mass_unit_cgs.v
    temperature = epsilon / tmp / yt.units.boltzmann_constant_cgs.v * data["density"]

    # use the approximate temperature to calculate the gammaH2
    gammaH2 = get_gammaH2(temperature)
    for s in species_h2:
        tmp += data[s] * (1 / (gammaH2 - 1)) / yt.units.atomic_mass_unit_cgs.v
        tmp -= data[s] * (1 / (7.0 / 5.0 - 1)) / yt.units.atomic_mass_unit_cgs.v
    temperature = epsilon / tmp / yt.units.boltzmann_constant_cgs.v * data["density"]

    return temperature


def calculate_instant_cooling_rates(traj, density, tidx=0):
    tmp_dict = {}
    for k, v in traj.items():
        if k not in ["successful", "t", "dt"]:
            if v.ndim == 2:
                tmp_dict[k] = np.copy(v[:, tidx])
            else:
                tmp_dict[k] = np.copy(v)
    tmp_dict["density"] = density
    tmp_dict["mdensity"] = density * 1.67e-24
    tmp_dict["z"] = np.ones_like(tmp_dict["ge"]) * 20

    tmp_dict["H2_1"] /= 2.01588
    tmp_dict["H2_2"] /= 2.01588

    tmp_dict["He_1"] /= 4.0
    tmp_dict["He_2"] /= 4.0
    tmp_dict["He_3"] /= 4.0

    T = calculate_temperature_from_dict(tmp_dict)
    tmp_dict["T"] = T

    mdensity = tmp_dict["mdensity"]
    tau = (mdensity / 3.3e-8) ** 2.8
    tau = np.maximum(tau, 1e-4)

    cie_fudge = np.minimum(1, (1 - np.exp(-tau)) / tau)
    tau = (mdensity / 3.3e-6) ** 8
    tau = np.maximum(tau, 1e-5)
    cie_fudge *= np.minimum(1, (1 - np.exp(-tau)) / tau)

    all_rates = {}

    combined = np.zeros_like(mdensity) * yt.units.g * yt.units.cm**-3

    for ca in cooling_registry.items():
        all_rates[ca[0]] = evaluate_cooling_terms(ca[1], tmp_dict)
        if ca[0] in ["gloverabel08"]:
            #         print("here")
            all_rates[ca[0]] *= np.minimum(
                1.0, pow((0.76 * mdensity / (1.34e-14)), -0.45)
            )
        if ca[0] not in ["h2formation", "h2formation_extra"]:
            all_rates[ca[0]] *= cie_fudge

        if ca[0] not in ["h2formation_extra", "gammah"]:
            combined += all_rates[ca[0]]

    all_rates["chemical"] = combined
    all_rates["T"] = T
    return all_rates


def add_net_heating(field, data):
    return (
        data[("gas", "h2formation")]
        + data[("gas", "gloverabel08")]
        + data[("gas", "cie_cooling")]
    )


def add_thermal_timescale(field, data):
    return np.abs(
        data[("gas", "specific_thermal_energy")]
        / (data[("gas", "net_heating_over_tdyn")] / data[("gas", "dynamical_time")])
    )


def get_gammaH2(T):
    """Get Temperature dependent gamma for molecular hydrogen"""
    if isinstance(T, yt.units.unyt_array):
        T0 = T.v ** (1 / 6.5)
    else:
        T0 = T ** (1 / 6.5)
    a0 = 64.2416
    a1 = -9.36334
    a2 = -0.377266
    a3 = 69.8091
    a4 = 0.0493452
    a5 = 2.28021
    a6 = 0.115357
    a7 = 0.114188
    a8 = 2.90778
    a9 = 0.689708

    gammaH2_expr = (
        np.exp(-a0 * T0**a1) * (a2 + T0**-a3)
        + a4 * np.exp(-((T0 - a5) ** 2) / a6)
        + a7 * np.exp(-((T0 - a8) ** 2) / a9)
        + 5.0 / 3.0
    )
    return gammaH2_expr


def get_new_field_function(field_name):
    def assign_dengo_field_unit(field, data):
        return data[field_name] * data.ds.mass_unit / data.ds.length_unit**3

    return assign_dengo_field_unit


def get_new_fraction_function(field_name):
    def assign_dengo_fraction(field, data):
        return (
            data[field_name]
            * data.ds.mass_unit
            / data.ds.length_unit**3
            / data["density"]
        )

    return assign_dengo_fraction


def get_grackle_fraction(field_name):
    def assign_grackle_fraction(field, data):
        return data[field_name] / data["density"]

    return assign_grackle_fraction


def calculate_temperature(field, data):
    """Evaluate temperature in the presence of molecular hydrogen"""
    epsilon = data["specific_thermal_energy"]

    species = [
        "HI_Density",
        "HII_Density",
        "HeI_Density",
        "HeII_Density",
        "HeIII_Density",
        "Electron_Density",
        "HM_Density",
    ]
    m_species = [1.00794, 1.00794, 4, 4, 4, 1.00794, 1.00794]
    species_h2 = ["H2I_Density", "H2II_Density"]

    # calculate sum of n / (gamma - 1)
    tmp = 0
    gamma = 5 / 3
    for i, s in enumerate(species):
        tmp += data[s] * (1 / (gamma - 1)) / m_species[i] / yt.units.atomic_mass_unit
    gammaH2 = 7.0 / 5.0
    for s in species_h2:
        tmp += data[s] * (1 / (gammaH2 - 1)) / 2.01588 / yt.units.atomic_mass_unit
    temperature = epsilon / tmp / yt.units.boltzmann_constant * data["density"]

    # use the approximate temperature to calculate the gammaH2
    gammaH2 = get_gammaH2(temperature)
    for s in species_h2:
        tmp += data[s] * (1 / (gammaH2 - 1)) / 2.01588 / yt.units.atomic_mass_unit
        tmp -= data[s] * (1 / (7.0 / 5.0 - 1)) / 2.01588 / yt.units.atomic_mass_unit
    temperature = epsilon / tmp / yt.units.boltzmann_constant * data["density"]

    # use the approximate temperature to calculate the gammaH2
    gammaH2 = get_gammaH2(temperature)
    for s in species_h2:
        tmp += data[s] * (1 / (gammaH2 - 1)) / 2.01588 / yt.units.atomic_mass_unit
        tmp -= data[s] * (1 / (7.0 / 5.0 - 1)) / 2.01588 / yt.units.atomic_mass_unit
    temperature = epsilon / tmp / yt.units.boltzmann_constant * data["density"]
    return temperature


def add_dengo_fields(dataset):
    """add chemical abundance/density as derived fields, correct for temperature"""
    for old_name, new_name in zip(dengo_species_field, grackle_species_field):
        dataset.add_field(
            ("gas", new_name),
            function=get_new_field_function(old_name),
            sampling_type="cell",
            units="g/ cm**3",
            force_override=True,
        )

        fraction_name = new_name.split("_")[0] + "_fraction"

        dataset.add_field(
            ("gas", fraction_name),
            function=get_new_fraction_function(old_name),
            sampling_type="cell",
            units="dimensionless",
            force_override=True,
        )


def add_grackle_fields(dataset):
    """add abundance fraction as derived fields, correct for temperature"""
    for old_name, new_name in zip(dengo_species_field, grackle_species_field):
        fraction_name = new_name.split("_")[0] + "_fraction"

        dataset.add_field(
            ("gas", fraction_name),
            function=get_grackle_fraction(new_name),
            sampling_type="cell",
            units="dimensionless",
            force_override=True,
        )


def load_dengo_dataset(filename, thermal_process=cooling_registry.keys()):
    """add additional fields and correct for temperature in dengo
    Parameters
    ----------
        filename: input file
    Returns
        yt dataset
    -------
    """
    dengo_ds = yt.load(filename)

    add_dengo_fields(dengo_ds)

    dengo_ds.add_field(
        ("gas", "temperature"),
        function=calculate_temperature,
        sampling_type="cell",
        units="K",
        force_override=True,
    )

    for ca in thermal_process:
        dengo_ds.add_field(
            ("gas", ca),
            function=add_dengo_thermal_fields(ca),
            sampling_type="cell",
            units="erg/g",
        )

    dengo_ds.add_field(
        ("gas", "net_heating_over_tdyn"),
        function=add_net_heating,
        sampling_type="cell",
        units="erg/g",
    )

    dengo_ds.add_field(
        ("gas", "thermal_timescale"),
        function=add_thermal_timescale,
        sampling_type="cell",
        units="s",
    )
    return dengo_ds


def load_grackle_dataset(filename, thermal_process=cooling_registry.keys()):

    grackle_ds = yt.load(filename)
    add_grackle_fields(grackle_ds)

    grackle_ds.add_field(
        ("gas", "temperature"),
        function=calculate_temperature,
        sampling_type="cell",
        units="K",
        force_override=True,
    )

    for ca in thermal_process:
        grackle_ds.add_field(
            ("gas", ca),
            function=add_dengo_thermal_fields(ca),
            sampling_type="cell",
            units="erg/g",
        )

    grackle_ds.add_field(
        ("gas", "net_heating_over_tdyn"),
        function=add_net_heating,
        sampling_type="cell",
        units="erg/g",
    )

    grackle_ds.add_field(
        ("gas", "thermal_timescale"),
        function=add_thermal_timescale,
        sampling_type="cell",
        units="s",
    )
    return grackle_ds
