import math

# --- Task 2: IBC intracellular compartment, closed form -------------------
psi_max_in = 1.3        # /h  (PNAS 2004 median of 1.19-1.39)
k_A_in     = 0.0007     # /h
k_exf      = 0.0        # /h
k_emerge   = 0.065      # /h
k_invade   = 0.01       # /h
a_bar      = -0.283     # /h  (Table 9.4.23)

d = psi_max_in - k_A_in - k_exf - k_emerge
tr  = a_bar + d
det = a_bar*d - k_emerge*k_invade
disc = tr*tr - 4*det
Lam = (tr + math.sqrt(disc))/2
Lam_minus = (tr - math.sqrt(disc))/2

print("=== Task 2.1  net intracellular growth / dominant eigenvalue ===")
print(f"d        = {d:+.4f} /h")
print(f"tr       = {tr:+.4f} /h   det = {det:+.6f} /h^2   disc = {disc:.6f}")
print(f"Lambda+  = {Lam:+.4f} /h   (doubling time ln2/L = {math.log(2)/Lam:.3f} h)")
print(f"Lambda-  = {Lam_minus:+.4f} /h")
print(f"Lambda+ - d = {Lam-d:+.6f} /h  -> fate set by intracellular compartment")
print()

# --- Table 9.4.50 : single-lever value required for d = 0 -----------------
print("=== Table 9.4.50  single-lever requirement for d = 0 ===")
print(f"psi_max_in must FALL to <= {k_A_in + k_exf + k_emerge:.4f} /h"
      f"   (IBC band floor 0.69 -> unreachable)")
print(f"k_exf      must RISE to >= {psi_max_in - k_A_in - k_emerge:.4f} /h")
print(f"k_emerge   must RISE to >= {psi_max_in - k_A_in - k_exf:.4f} /h")
print(f"k_A_in     must RISE to >= {psi_max_in - k_exf - k_emerge:.4f} /h")
span = 1.70   # psi_max - psi_min at rho_in = 1
need = psi_max_in - k_exf - k_emerge
print(f"  k_A_in ceiling (psi_max - psi_min, rho_in=1) = {span:.2f} /h"
      f"  -> need {100*need/span:.1f}% of full span")
kappa = 2.95
frac = need/span
x = (frac/(1-frac))**(1/kappa)
print(f"  Hill inversion (kappa={kappa}): normalised exposure x = eta*A/MIC_R >= {x:.4f}")
print(f"  at MIC_R(0)=163.2 ug/mL, A_peak=300 ug/mL -> eta >= {x*163.2/300:.4f}")
print("  [Hill form assumed psi = psi_max - (psi_max-psi_min)*x^k/(1+x^k);")
print("   must be confirmed against sunrise_model.py before publishing]")
print()

# --- Task 8: mechanical-efficiency axis ----------------------------------
print("=== Task 8  eta_mech scaling (closed form, no ODE) ===")
print("N_0 = eta_mech * phi * m * P  =>  m_inf, D, particles/course all ~ 1/eta_mech")
m_inf_ref = 411.0   # particles/CFU at eta_mech = 1, Ver.22 baseline D = 3.0e11
m_1h_ref  = 86.0    # particles/CFU at tau_P = 1 h, eta_mech = 1
for em in (1.0, 0.5, 0.2):
    print(f"  eta_mech = {em:>4}:  m_inf = {m_inf_ref/em:7.0f} particles/CFU"
          f"   | tau_P=1h requirement = {m_1h_ref/em:5.0f} particles/CFU"
          f"   {'(within 10-100)' if m_1h_ref/em <= 100 else '(OUTSIDE 10-100)'}")
print("  NOTE: m_inf column is anchored to the Ver.22 D = 3.0e11 baseline and")
print("        must be rescaled once Task 1 fixes the new D.")
