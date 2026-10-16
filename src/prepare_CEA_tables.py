import cea
import cea.matlab as ceam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# CEA configuration
Pc = cea.units.psi_to_bar(340)
expansion_ratio = 4
T_reac = 298.15

# Mixture ratio range
OF = np.linspace(2, 12, 101)

paraffin = cea.Reactant(
    name="paraffin",
    formula={"C": 32, "H": 66},
    molecular_weight=450.8664,
    enthalpy=-967.8,
    enthalpy_units="kJ/mol",
    temperature=298.15
)

# Reactants
reactants = [paraffin, "N2O"]
fuel = [1.0, 0.0]
oxid = [0.0, 1.0]

# Generate isobaric table
results = []

for of in OF:
    sol = ceam.rocket_solve(
        reactants,
        pc=Pc,
        supar=expansion_ratio,
        T_reac=T_reac,
        fuel_amounts=fuel,
        oxid_amounts=oxid,
        of_ratio=float(of),
        iac=True
    )

    if not sol.converged:
        raise RuntimeError(f"CEA failed at O/F = {of:.3f}")

    # Extract chamber and nozzle exit properties
    ve = float(sol.Isp[-1])        # Effective exhaust velocity [m/s]
    T0 = float(sol.T[0])           # Chamber temperature [K]
    cstar = float(sol.c_star[0])   # Characteristic velocity [m/s]
    pe = float(sol.P[-1])          # Exit pressure [bar]

    results.append([of, ve, T0, cstar, pe])

# Export table
table = pd.DataFrame(
    results,
    columns=["OF", "ve", "T0", "cstar", "pe"]
)

table.to_csv("CEA_Isobaric.csv", index=False)
print(table.to_string(index=False))

# Plot exhaust velocity
plt.figure()
plt.plot(table["OF"], table["ve"])
plt.xlabel("O/F")
plt.ylabel("Effective Exhaust Velocity [m/s]")
plt.title("CEA Isobaric Performance - 340 psia")
plt.grid(True)
plt.tight_layout()
plt.show()