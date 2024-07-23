from itertools import product
from math import sqrt
from scipy.optimize import brentq
from . import ideal

from ..common import *

# Van der Waals equation for real gases:
# (P + a * n^2 / V^2) * (V - n * b) = n * R * T

VDW_PARAMETERS = [
    #a,      b,       MW
    (0.1382,  3.186e-5, 31.9988),     # O2
    (0.00346, 2.380e-5, 4.0020602),   # He
    (0.1370,  3.870e-5, 28.01348)     # N2
]

def mixing_rules(fractions: tuple[float]) -> tuple[float]:
    # return the VDW parameters (a, b + molar weight) for the mix
    # using the general mixing rules (see wikipedia)
    a_mix = sum(fractions[i] * fractions[j] * sqrt(VDW_PARAMETERS[i][0] * VDW_PARAMETERS[j][0]) for i,j in product(range(3), repeat=2))
    b_mix = sum(fractions[i] * VDW_PARAMETERS[i][1] for i in range(3))
    mw_mix = sum(fractions[i] * VDW_PARAMETERS[i][2] for i in range(3))
    return a_mix, b_mix, mw_mix

def get_pressure(n:float, V:float, T:float, fractions:tuple[float], use_SI:bool = False) -> float:
    if use_SI:
        return _get_pressure_SI(n, V, T, fractions)
    return _get_pressure_SI(n,
                            V * LITER_TO_CUBIC_METER,
                            T + CELSIUS_TO_KELVIN,
                            fractions) / BAR_TO_PASCAL

def _get_pressure_SI(n:float, V:float, T:float, fractions:tuple[float]) -> float:
    a, b, _ = mixing_rules(fractions)
    # VDW is explicit in pressure
    return ((n * R_IDEAL_GASES * T) / (V - n * b)) - (a * (n / V)**2.0)

def get_moles(p:float, V:float, T:float, fractions:tuple[float], SI:bool = False) -> float:
    if SI:
        return _get_moles_SI(p, V, T)
    return _get_moles_SI(p * BAR_TO_PASCAL,
                         V * LITER_TO_CUBIC_METER,
                         T + CELSIUS_TO_KELVIN,
                         fractions)

def _get_moles_SI(p:float, V:float, T:float, fractions:tuple[float]) -> float:
    a, b, _ = mixing_rules(fractions)
    # we'll need some magic here
    def vdw_solve_n(n: float) -> float:
        # maybe use a partial() here...
        return (p + a * (n / V)**2.0) * (V - n * b) - n * R_IDEAL_GASES * T
    root_estimate = ideal.get_moles(p, V, T, SI=True)
    root = brentq(vdw_solve_n, root_estimate / 2.0, root_estimate * 2.0)
    return root

if __name__ == "__main__":
    AIR = (0.2095, 0.0, 0.7905)
    print(mixing_rules(AIR))

