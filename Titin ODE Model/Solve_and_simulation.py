import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#Parameters and optimized interaction coefficients from Table4.
alpha = [0.95, 0.6, 0.5, 0.5, 0.7, 0.35, 0.3, 0.4]
K = [0.45, 0.25, 0.2, 0.2, 0.3, 0.12, 0.03, 0.08]
beta12 = 0.15
beta14 = 0.1
beta25 = 0.08
gamma31 = 0.05
gamma13 = 0.03
epsilon34 = 0.2
eta45 = 0.1
theta56 = 0.08
mu67 = 0.04
mu76 = 0.02
nu78 = 0.03
nu87 = 0.02
rho81 = 0.06
#Definition DE set
def equations(s, x):
    x1, x2, x3, x4, x5, x6, x7, x8 = x
    dx1 = alpha[0]*x1*(1-x1/K[0]) + beta12*x1*x2 - beta14*x1*x4 + gamma31*x3 - gamma13*x1 + rho81*x8
    dx2 = alpha[1]*x2*(1-x2/K[1]) + beta12*x1*x2 - beta25*x2*x5
    dx3 = alpha[2]*x3*(1-x3/K[2]) - gamma31*x3 + gamma13*x1 + epsilon34*x3*x4
    dx4 = alpha[3]*x4*(1-x4/K[3]) - beta14*x1*x4 + epsilon34*x3*x4 + eta45*x5
    dx5 = alpha[4]*x5*(1-x5/K[4]) - beta25*x2*x5 + eta45*x4 + theta56*x6
    dx6 = alpha[5]*x6*(1-x6/K[5]) + theta56*x5 - mu67*x6 + mu76*x7
    dx7 = alpha[6]*x7*(1-x7/K[6]) + mu67*x6 - nu78*x7 + nu87*x8
    dx8 = alpha[7]*x8*(1-x8/K[7]) + nu78*x7 - nu87*x8 + rho81*x1
    return [dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8]

#Initial conditions written in the caption of Figure1.
x0 = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.02, 0.05]
#Numerical solution with RK45 algorithm.
s_span = (0, 100)
s_eval = np.linspace(0, 100, 1000)
sol = solve_ivp(equations, s_span, x0, t_eval=s_eval, method='RK45')
#Figure1
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
#Panel a
labels = ['Hyd+', 'Hyd-', 'Pos+', 'Neg-', 'Polar', 'Gly', 'Cys', 'Aro']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf']
for i in range(8):
    ax1.plot(sol.t, sol.y[i], label=labels[i], color=colors[i], linewidth=2)
ax1.set_xlabel('Model coordinate s', fontsize=12)
ax1.set_ylabel('Model abundance', fontsize=12)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
ax1.set_xlim(0, 100)
ax1.text(0.02, 0.95, '(a)', transform=ax1.transAxes, fontsize=14, fontweight='bold')
#Panel b
ax2.plot(sol.y[0], sol.y[5], 'b-', linewidth=2)
ax2.plot(x0[0], x0[5], 'ro', markersize=10)
ax2.plot(sol.y[0,-1], sol.y[5,-1], 'gs', markersize=10)
ax2.set_xlabel('x₁ (Hyd+)', fontsize=12)
ax2.set_ylabel('x₆ (Gly)', fontsize=12)
ax2.legend(['Trajectory', 'Initial condition', 'Equilibrium'])
ax2.text(0.02, 0.95, '(b)', transform=ax2.transAxes, fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figure1.png', dpi=300, bbox_inches='tight')
plt.show()
#equilibrium values
print("equilibrium values:")
for i in range(8):
    print(f"x{i+1} = {sol.y[i,-1]:.3f}")