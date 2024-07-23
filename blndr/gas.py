from dataclasses import dataclass
from typing import Self
from math import floor

@dataclass
class GasMix:
    O2: float
    He: float

    @classmethod
    def from_moles(cls, O2:float = 0.0, He:float = 0.0, N2:float = 0.0) -> Self:
        m_tot = O2 + He + N2
        if m_tot == 0.0:
            return cls(0.0, 0.0)
        return cls(O2/m_tot, He/m_tot)

    @property
    def N2(self) -> float:
        return 1.0 - self.O2 - self.He

    def MOD(self, ppO2:float = 1.4, round:bool = True) -> float:
        assert(ppO2 < 2.0) # ppO2 > 2 bar is (very) dangerous
        mod = 10 * ((ppO2 / self.O2) - 1)
        if round:
            return floor(mod) * 1.0
        return mod
    
    def EAD(self, depth: float|None = None, round:bool = True) -> float:
        if depth is None:
            depth = self.MOD()
        ead = (depth + 10.0) * self.N2 / 0.7905 - 10
        if round:
            return floor(ead) * 1.0
        return ead

    @property
    def fractions(self) -> tuple[float]:
        return (self.O2, self.He, self.N2)

    def __post_init__(self) -> None:
        if self.O2 < 0.0 or self.O2 > 1.0:
            raise ValueError(f"incorrect O2 fraction: {self.O2}")
        if self.He < 0.0 or self.He > 1.0:
            raise ValueError(f"incorrect He fraction: {self.He}")
        if self.O2 + self.He > 1.0:
            raise ValueError(f"sum of O2 ({self.O2}) and He ({self.He}) fractions cannot exceed 100%")
    
    def __str__(self) -> str:
        return f"Gas mix: O2 {self.O2 * 100:.2f}%, He {self.He * 100:.2f}%, N2 {self.N2 * 100:.2f}%"


AIR = GasMix(0.2095, 0.0)
EAN32 = GasMix(0.32, 0.0)
EAN36 = GasMix(0.36, 0.0)
EAN40 = GasMix(0.40, 0.0)
O2 = GasMix(1.0, 0.0)
He = GasMix(0.0, 1.0)
N2 = GasMix(0.0, 0.0)


if __name__ == "__main__":
    print(AIR, AIR.MOD())

