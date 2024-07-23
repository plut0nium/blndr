#! /bin/env python

from blndr import gas
from blndr.blndr import make_fill_plan

if __name__ == "__main__":
    volume = 10.0
    temperature = 20.0

    pressure_initial = 50.0
    pressure_target = 200.0

    gas_initial = gas.O2
    gas_target = gas.GasMix(0.2, 0.15)
    gas_top_up = gas.AIR

    print(f"[START] {gas_initial} @ {pressure_initial} bar")
    print(f"[TARGET] {gas_target} @ {pressure_target} bar")

    steps = make_fill_plan(gas_target, pressure_target,
                           gas_initial, pressure_initial,
                           volume, temperature, gas_top_up)

    print("[STEPS]")
    for i, s in enumerate(steps):
        print(f" {i+1}. {s}")


