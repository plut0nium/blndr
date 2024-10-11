from . import gas
from .gas_model import ideal, vanderwaals
from enum import Enum


class FillStep(Enum):
    BLEED = "BLEED"
    FILL = "FILL"
    TOP_UP = "TOP-UP"


def make_fill_plan(gas_target:gas.GasMix, pressure_target:float = 200.0,
                   gas_initial:gas.GasMix|None = None, pressure_initial:float = 0.0,
                   volume:float = 10.0, temperature:float = 20.0,
                   gas_top_up:gas.GasMix|None = None) -> list[str]:
    # default to AIR
    if gas_initial is None:
        gas_initial = gas.AIR
    if gas_top_up is None:
        gas_top_up = gas.AIR

    steps = []

    moles_total_initial = ideal.get_moles(pressure_initial, volume, temperature)
    moles_initial = tuple(x * moles_total_initial for x in gas_initial.fractions)

    moles_total_target = ideal.get_moles(pressure_target, volume, temperature)
    moles_target = tuple(x * moles_total_target for x in gas_target.fractions)

    moles_current = list(moles_initial)

    # check if bleeding is required
    if any(moles_current[i] > moles_target[i] for i in range(3)):
        if any(moles_current[i] > 0.0 and moles_target[i] == 0.0 for i in range(3)):
            # we have an unwanted gas -> empty completely
            moles_current = [0.0 for _ in range(3)]
        else:    
            bleed_ratio = max(moles_current[i] / moles_target[i] for i in range(3) if moles_target[i] > 0.0)
            for i in range(3):
                moles_current[i] *= 1.0 / bleed_ratio
        #print("[BLEED] to", ideal.get_pressure(sum(moles_current), volume, temperature), "bar")

    if (moles_target[2] - moles_current[2]) > 0 and gas_top_up.N2 == 0.0:
        # we need N2 but none is available
        raise ValueError("A source of N2 (top-up gas) is required for this Mix")
    if gas_top_up.N2 > 0.0:
        moles_total_top_up = (moles_target[2] - moles_current[2]) / gas_top_up.N2
    else:
        # no top-up if N2 fraction = 0
        moles_total_top_up = 0.0
    moles_top_up = tuple(x * moles_total_top_up for x in gas_top_up.fractions)

    # TODO: check that the fraction of N2 is high enough to reach target without overfilling
    # ideas: gas_top_up.N2 > gas_target.N2

    moles_O2_fill = moles_target[0] - moles_current[0] - moles_top_up[0]
    moles_He_fill = moles_target[1] - moles_current[1] - moles_top_up[1]

    if moles_O2_fill < 0.0 or moles_He_fill < 0.0:
        # quantities of O2 and He would be exceeded -> bleeding
        moles_fill = [moles_O2_fill, moles_He_fill]
        bleed_ratio = max(-moles_fill[i]/moles_current[i] for i in range(2) if moles_current[i] > 0)
        assert(bleed_ratio <= 1.0)
        for i in range(3):
            moles_current[i] *= 1.0 - bleed_ratio
        # print("[BLEED] (2) to", ideal.get_pressure(sum(moles_current), volume, temperature), "bar")

        # recalculate top-up & fill
        if gas_top_up.N2 > 0.0:
            moles_total_top_up = (moles_target[2] - moles_current[2]) / gas_top_up.N2
        else:
            # no top-up if N2 fraction = 0
            moles_total_top_up = 0.0
        moles_top_up = tuple(x * moles_total_top_up for x in gas_top_up.fractions)

        moles_O2_fill = moles_target[0] - moles_current[0] - moles_top_up[0]
        moles_He_fill = moles_target[1] - moles_current[1] - moles_top_up[1]

    if moles_O2_fill < 0.0 or moles_He_fill < 0.0:
        raise ValueError("Unable to blend this mix using the available gases")

    if ideal.get_pressure(sum(moles_current), volume, temperature) < pressure_initial:
        steps.append(f"[BLEED] down to {ideal.get_pressure(sum(moles_current), volume, temperature):.2f} bar")

    if moles_O2_fill > 0.0:
        moles_current[0] += moles_O2_fill
        gas_current = gas.GasMix.from_moles(*moles_current)
        # print("O2 FILL")
        # print(gas_current)
        # print(ideal.get_pressure(sum(moles_current), volume, temperature), "bar")
        steps.append(f"[FILL] with O2 up to {ideal.get_pressure(sum(moles_current), volume, temperature):.2f} bar")

    if moles_He_fill > 0.0:
        moles_current[1] += moles_He_fill
        gas_current = gas.GasMix.from_moles(*moles_current)
        # print("He FILL")
        # print(gas_current)
        # print(ideal.get_pressure(sum(moles_current), volume, temperature), "bar")
        steps.append(f"[FILL] with He up to {ideal.get_pressure(sum(moles_current), volume, temperature):.2f} bar")

    moles_current = [c + t for c, t in zip(moles_current, moles_top_up)]
    gas_current = gas.GasMix.from_moles(*moles_current)
    # print("TOP-UP")
    # print(gas_current)
    # print(ideal.get_pressure(sum(moles_current), volume, temperature), "bar")
    # print(vanderwaals.get_pressure(sum(moles_current), volume, temperature, gas_current.fractions), "bar")
    steps.append(f"[TOP-UP] to {ideal.get_pressure(sum(moles_current), volume, temperature):.2f}bar")

    return steps