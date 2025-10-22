"""
Run tests for 1D linear advection:
- Compare FTCS and Upwind schemes
- Check convergence for Upwind
- Investigate CFL effects
"""

import os
import numpy as np
import advection_template as adv

def main():
    # Ensure results folder exists
    if not os.path.exists("results"):
        os.makedirs("results")

    testproblem = adv.get_testproblem()

    # --- Part D: FTCS and Upwind comparison ---
    print("=== FTCS Scheme ===")
    parameters = adv.define_default_parameters()
    parameters["Nrefine"] = 2  # N = 40, 80, 160
    adv.my_driver("FTCS", testproblem, parameters, parameters["N"])

    print("=== Upwind Scheme ===")
    adv.my_driver("upwind", testproblem, parameters, parameters["N"])

    # --- Part E: CFL influence for Upwind ---
    CFL_values = [0.8, 1.0, 1.2]
    N = 40

    for CFL in CFL_values:
        print(f"\n--- Upwind Scheme with CFL = {CFL} ---")
        parameters = adv.define_default_parameters()
        parameters["CFL"] = CFL
        parameters["Nrefine"] = 0  # single run
        adv.my_driver("upwind", testproblem, parameters, N)

    # --- Part 3: Show convergence table ---
    print("\nConvergence table stored in 'results/error.txt':")
    if os.path.exists("results/error.txt"):
        with open("results/error.txt") as f:
            print(f.read())
    else:
        print("No convergence table found. Run with Nrefine > 0.")

if __name__ == "__main__":
    main()
