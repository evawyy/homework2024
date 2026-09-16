#program for problem 5
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye

# =====================================================
# Parameters
# =====================================================

n = 30
h = 1.0 / (n + 1)
N = n * n
maxit = 300

# =====================================================
# Construct 2D Poisson matrix
# =====================================================

e = np.ones(n)

T = diags(
    [-e, 4*e, -e],
    [-1, 0, 1],
    shape=(n, n)
)

I = eye(n)

A = kron(I, T) + kron(
    diags([-e, -e], [-1, 1], shape=(n, n)),
    -I
)

A = A.toarray() / h**2

# =====================================================
# Right-hand side
# =====================================================

np.random.seed(0)

x_exact = np.random.randn(N)

b = A @ x_exact

# =====================================================
# Young theory
# =====================================================

rho_J = np.cos(np.pi/(n+1))

rho_GS = rho_J**2

omega_opt = 2/(1 + np.sqrt(1-rho_J**2))

rho_SOR = omega_opt - 1

print("Theoretical quantities")
print("----------------------")
print(f"rho(Jacobi) = {rho_J:.8f}")
print(f"rho(GS)     = {rho_GS:.8f}")
print(f"omega_opt   = {omega_opt:.8f}")
print(f"rho(SOR)    = {rho_SOR:.8f}")

# =====================================================
# Splitting
# =====================================================

D = np.diag(np.diag(A))
L = -np.tril(A, -1)
U = -np.triu(A, 1)

# =====================================================
# Relative residual
# =====================================================

r0_norm = np.linalg.norm(b)

def relres(x):
    return np.linalg.norm(b - A @ x) / r0_norm

# =====================================================
# Jacobi
# =====================================================

resJ = []

x = np.zeros(N)

for k in range(maxit):

    resJ.append(relres(x))

    x = np.linalg.solve(
        D,
        (L + U) @ x + b
    )

# =====================================================
# Gauss-Seidel
# =====================================================

resGS = []

x = np.zeros(N)

M = D - L
Nmat = U

for k in range(maxit):

    resGS.append(relres(x))

    x = np.linalg.solve(
        M,
        Nmat @ x + b
    )

# =====================================================
# SOR
# =====================================================

resSOR = []

x = np.zeros(N)

M = (1/omega_opt)*D - L

Nmat = ((1/omega_opt)-1)*D + U

for k in range(maxit):

    resSOR.append(relres(x))

    x = np.linalg.solve(
        M,
        Nmat @ x + b
    )

# =====================================================
# Figure 1
# Residual curves
# =====================================================

plt.figure(figsize=(8,6))

plt.semilogy(
    resJ,
    linewidth=2,
    label='Jacobi'
)

plt.semilogy(
    resGS,
    linewidth=2,
    label='Gauss-Seidel'
)

plt.semilogy(
    resSOR,
    linewidth=2,
    label=rf'SOR ($\omega={omega_opt:.4f}$)'
)

plt.xlabel('Iteration')
plt.ylabel(r'$\|r_k\|_2/\|r_0\|_2$')

plt.title('Residual Convergence')

plt.grid(True)

plt.legend()

plt.tight_layout()

plt.savefig(
    "figure1_residual_convergence.png",
    dpi=300,
    bbox_inches="tight"
)

print("Saved: figure1_residual_convergence.png")
# =====================================================
# Figure 2
# Verify Young theory
# =====================================================

qJ = np.array(resJ[1:]) / np.array(resJ[:-1])

qGS = np.array(resGS[1:]) / np.array(resGS[:-1])

qSOR = np.array(resSOR[1:]) / np.array(resSOR[:-1])

plt.figure(figsize=(8,6))

plt.plot(
    qJ,
    linewidth=2,
    label='Jacobi'
)

plt.plot(
    qGS,
    linewidth=2,
    label='Gauss-Seidel'
)

plt.plot(
    qSOR,
    linewidth=2,
    label='SOR'
)

# Young theory lines

plt.axhline(
    rho_J,
    color='C0',
    linestyle='--',
    linewidth=2,
    label=rf'$\rho_J={rho_J:.6f}$'
)

plt.axhline(
    rho_GS,
    color='C1',
    linestyle='--',
    linewidth=2,
    label=rf'$\rho_{{GS}}={rho_GS:.6f}$'
)

plt.axhline(
    rho_SOR,
    color='C2',
    linestyle='--',
    linewidth=2,
    label=rf'$\rho_{{SOR}}={rho_SOR:.6f}$'
)

plt.xlabel('Iteration')

plt.ylabel(
    r'$\|r_{k+1}\|_2/\|r_k\|_2$'
)

plt.title(
    "Verification of Young's Theory"
)

plt.grid(True)

plt.legend()

plt.tight_layout()

plt.savefig(
    "figure2_verify_young.png",
    dpi=300,
    bbox_inches="tight"
)

print("Saved: figure2_verify_young.png")


# =====================================================
# Numerical comparison
# =====================================================

print("\nObserved asymptotic factors")
print("---------------------------")

print(
    f"Jacobi : {np.mean(qJ[-20:]):.8f}"
)

print(
    f"GS     : {np.mean(qGS[-20:]):.8f}"
)

print(
    f"SOR    : {np.mean(qSOR[-20:]):.8f}"
)
