from itertools import product
from math import sqrt
from scipy.optimize import brentq

from . import ideal
from ..common import *
from ..gas import GasMix

# Van der Waals equation for real gases:
# (P + a * n^2 / V^2) * (V - n * b) = n * R * T

VDW_PARAMETERS = [
    #a,      b,       MW
    (0.1382,  3.186e-5, 31.9988),     # O2
    (0.00346, 2.380e-5, 4.0020602),   # He
    (0.1370,  3.870e-5, 28.01348)     # N2
]

def mixing_rules(fractions: tuple[float]) -> tuple[float]:
    '''Apply mixing rules for VDW parameters.

    Calculate the VDW parameters (a, b + molar weight) for the mix,
    in SI units, using the general mixing rules (see wikipedia).
    
    Args:
        fractions (tuple): the gas mix for which VDW parameters should be calculated,
                           provided as a tuple[float] containing the molar fractions
                           of resp. O2, He and N2
    
    Returns:
        a 3-values tuple, containing:
        - a (float): the interaction parameter [m^6 * Pa * mol^-2]
        - b (float): the molar covolume parameter [m^3 * mol^-1]
        - mw (float): the averaged molar weight of the mix [g * mol^-1]
    '''
    a_mix = sum(fractions[i] * fractions[j] * sqrt(VDW_PARAMETERS[i][0] * VDW_PARAMETERS[j][0]) for i,j in product(range(3), repeat=2))
    b_mix = sum(fractions[i] * VDW_PARAMETERS[i][1] for i in range(3))
    mw_mix = sum(fractions[i] * VDW_PARAMETERS[i][2] for i in range(3))
    return a_mix, b_mix, mw_mix

def get_pressure(n:float, V:float, T:float, fractions:tuple[float], use_SI:bool = False) -> float:
    '''Wrapper function for unit management.'''
    if use_SI:
        return _get_pressure_SI(n, V, T, fractions)
    return _get_pressure_SI(n,
                            V * LITER_TO_CUBIC_METER,
                            T + CELSIUS_TO_KELVIN,
                            fractions) / BAR_TO_PASCAL

def _get_pressure_SI(n:float, V:float, T:float, fractions:tuple[float]) -> float:
    '''Use VDW equation to calculate the pressure of a gas.

    Calculate the pressure of a gas mix given the number of moles, the volume
    and the temperature, using the Van Der Waals equation of state.

    Note:
        Should not be called directly, but by get_pressure(...)

    Args:
        n (float): nomber of moles
        V (float): volume, in [m^3]
        T (float): temperature, in [K]
        fractions (tuple): the gas mix for which the pressure should be calculated,
                           provided as a tuple[float] containing the molar fractions
                           of resp. O2, He and N2

    Returns:
        pressure (float), in [Pa]
    '''
    a, b, _ = mixing_rules(fractions)
    # VDW is explicit in pressure -> easy
    return ((n * R_IDEAL_GASES * T) / (V - n * b)) - (a * (n / V)**2.0)

def get_moles(p:float, V:float, T:float, fractions:tuple[float], SI:bool = False) -> float:
    '''Wrapper function for unit management.'''
    if SI:
        return _get_moles_SI(p, V, T)
    return _get_moles_SI(p * BAR_TO_PASCAL,
                         V * LITER_TO_CUBIC_METER,
                         T + CELSIUS_TO_KELVIN,
                         fractions)

def _get_moles_SI(p:float, V:float, T:float, fractions:tuple[float]) -> float:
    '''Use VDW equation to calculate the number of moles of a gas.

    Calculate the number of moles of a gas mix given the pressure, the volume
    and the temperature, using a root-finding algorithms applied to the
    Van Der Waals equation of state.

    Note:
        Should not be called directly, but by get_moles(...)

    Args:
        p (float): pressure, in [Pa]
        V (float): volume, in [m^3]
        T (float): temperature, in [K]
        fractions (tuple): the gas mix for which the pressure should be calculated,
                           provided as a tuple[float] containing the molar fractions
                           of resp. O2, He and N2

    Returns:
        number of moles (float)
    '''
    a, b, _ = mixing_rules(fractions)
    # not explicit for n, so we'll need some scipy magic here
    def vdw_solve_n(n: float) -> float:
        # TODO: maybe use a partial() here...
        return (p + a * (n / V)**2.0) * (V - n * b) - n * R_IDEAL_GASES * T
    # use ideal gas law to get a first estimate
    root_estimate = ideal.get_moles(p, V, T, SI=True)
    # Brent's method requires you to bound the root search interval
    root = brentq(vdw_solve_n, root_estimate / 2.0, root_estimate * 2.0)
    return root

