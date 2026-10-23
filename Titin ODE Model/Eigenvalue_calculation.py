import numpy as np
#Parameters and optimized interaction coefficients from Table4.
alpha = np.array([0.95, 0.60, 0.50, 0.50, 0.70, 0.35, 0.30, 0.40])
K = np.array([0.45, 0.25, 0.20, 0.20, 0.30, 0.12, 0.03, 0.08])
beta_12 = 0.15
beta_14 = 0.10
beta_25 = 0.08
eps_34 = 0.20
gamma_31 = 0.05
gamma_13 = 0.03
eta_45 = 0.10
theta_56 = 0.08
mu_67 = 0.04
mu_76 = 0.02
nu_78 = 0.03
nu_87 = 0.02
rho_81 = 0.06
#Stable equilibrium pointfrom Table 5.
x_star = np.array([0.462, 0.268, 0.225, 0.253, 0.339, 0.165, 0.046, 0.123])
x1, x2, x3, x4, x5, x6, x7, x8 = x_star
#Constructing the Jacobian matrix at the equilibrium point(X*)
#Matrix elements based on partial derivatives of differential equations
J = np.array([
    [alpha[0]*(1 - 2*x1/K[0]) + beta_12*x2 - beta_14*x4 - gamma_13,
     beta_12*x1, gamma_31, -beta_14*x1, 0, 0, 0, rho_81],
    [beta_12*x2, alpha[1]*(1 - 2*x2/K[1]) + beta_12*x1 - \
     beta_25*x5, 0, 0, -beta_25*x2, 0, 0, 0],
    [gamma_13, 0, alpha[2]*(1 - 2*x3/K[2]) - gamma_31 + \
     eps_34*x4, eps_34*x3, 0, 0, 0, 0],
    [-beta_14*x4, 0, eps_34*x4, alpha[3] * \
        (1 - 2*x4/K[3]) - beta_14*x1 + eps_34*x3, eta_45, 0, 0, 0],
    [0, -beta_25*x5, 0, eta_45, alpha[4] * \
        (1 - 2*x5/K[4]) - beta_25*x2, theta_56, 0, 0],
    [0, 0, 0, 0, theta_56, alpha[5]*(1 - 2*x6/K[5]) - mu_67, mu_76, 0],
    [0, 0, 0, 0, 0, mu_67, alpha[6]*(1 - 2*x7/K[6]) - nu_78, nu_87],
    [rho_81, 0, 0, 0, 0, 0, nu_78, alpha[7]*(1 - 2*x8/K[7]) - nu_87]])
#Eigenvalue calculation
eigenvalues = np.linalg.eigvals(J)
idx = np.argsort(np.real(eigenvalues))
sorted_eigs = eigenvalues[idx]
for i, eig in enumerate(sorted_eigs):
    real_part = np.real(eig)
    imag_part = np.imag(eig)
    print(f"eigenvalue{i+1}:{real_part:8.4f}+{imag_part:8.4f}i")
print("-" * 30)
print("Trace:",np.trace(J))