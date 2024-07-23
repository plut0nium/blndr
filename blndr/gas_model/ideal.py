from ..common import *

# ideal gas law : p * V = n * R * T

def get_moles(p:float, V:float, T:float, fractions:tuple[float]|None = None, SI:bool = False) -> float:
    if SI:
        return _get_moles_SI(p, V, T)
    return _get_moles_SI(p * BAR_TO_PASCAL,
                         V * LITER_TO_CUBIC_METER,
                         T + CELSIUS_TO_KELVIN)

def _get_moles_SI(p:float, V:float, T:float) -> float:
    return (p * V) / (R_IDEAL_GASES * T)

def get_pressure(n:float, V:float, T:float, fractions:tuple[float]|None = None, SI:bool = False) -> float:
    if SI:
        return _get_pressure_SI(n, V, T)
    return _get_pressure_SI(n,
                            V * LITER_TO_CUBIC_METER,
                            T + CELSIUS_TO_KELVIN) / BAR_TO_PASCAL

def _get_pressure_SI(n:float, V:float, T:float) -> float:
    return (n * R_IDEAL_GASES * T) / V

