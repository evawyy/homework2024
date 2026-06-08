import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt


# ============================================================
# Basic iterative methods
# ============================================================

def bicg(A, b, x0=None, rhat0=None, maxit=200, xtrue=None):
    n = b.size
    x = np.zeros(n) if x0 is None else x0.copy()
    r = b - A @ x
    rt = r.copy() if rhat0 is None else rhat0.copy()

    p = r.copy()
    pt = rt.copy()

    res, xnorm = [], []

    rho_old = np.dot(rt, r)

    for k in range(maxit):
        Ap = A @ p
        ATpt = A.T @ pt

        alpha = rho_old / np.dot(pt, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rt = rt - alpha * ATpt

        res.append(np.linalg.norm(b - A @ x))
        if xtrue is not None:
            xnorm.append(np.linalg.norm(x) / np.linalg.norm(xtrue))

        rho_new = np.dot(rt, r)
        if abs(rho_new) < 1e-30:
            break

        beta = rho_new / rho_old
        p = r + beta * p
        pt = rt + beta * pt
        rho_old = rho_new

    return x, np.array(res), np.array(xnorm)


def cgs(A, b, x0=None, rhat0=None, maxit=200, xtrue=None):
    n = b.size
    x = np.zeros(n) if x0 is None else x0.copy()
    r = b - A @ x
    rt = r.copy() if rhat0 is None else rhat0.copy()

    p = r.copy()
    u = r.copy()
    rho_old = np.dot(rt, r)

    res, xnorm = [], []

    for k in range(maxit):
        Ap = A @ p
        alpha = rho_old / np.dot(rt, Ap)

        q = u - alpha * Ap
        x = x + alpha * (u + q)

        Auq = A @ (u + q)
        r = r - alpha * Auq

        res.append(np.linalg.norm(b - A @ x))
        if xtrue is not None:
            xnorm.append(np.linalg.norm(x) / np.linalg.norm(xtrue))

        rho_new = np.dot(rt, r)
        if abs(rho_new) < 1e-30:
            break

        beta = rho_new / rho_old
        u = r + beta * q
        p = u + beta * (q + beta * p)
        rho_old = rho_new

    return x, np.array(res), np.array(xnorm)


def bicgstab(A, b, x0=None, maxit=200, xtrue=None):
    n = b.size
    x = np.zeros(n) if x0 is None else x0.copy()
    r = b - A @ x
    rt = r.copy()

    rho_old = alpha = omega = 1.0
    v = np.zeros(n)
    p = np.zeros(n)

    res = []

    for k in range(maxit):
        rho_new = np.dot(rt, r)
        if abs(rho_new) < 1e-30:
            break

        beta = (rho_new / rho_old) * (alpha / omega)
        p = r + beta * (p - omega * v)

        v = A @ p
        alpha = rho_new / np.dot(rt, v)

        s = r - alpha * v
        t = A @ s

        omega = np.dot(t, s) / np.dot(t, t)

        x = x + alpha * p + omega * s
        r = s - omega * t

        res.append(np.linalg.norm(b - A @ x))

        if np.linalg.norm(r) < 1e-14:
            break

        rho_old = rho_new

    return x, np.array(res)


def cgne(A, b, x0=None, maxit=200):
    n = b.size
    x = np.zeros(n) if x0 is None else x0.copy()

    r = b - A @ x
    z = A.T @ r
    p = z.copy()

    res = []

    for k in range(maxit):
        Ap = A @ p
        alpha = np.dot(z, z) / np.dot(Ap, Ap)

        x = x + alpha * p
        r = r - alpha * Ap
        z_new = A.T @ r

        res.append(np.linalg.norm(b - A @ x))

        beta = np.dot(z_new, z_new) / np.dot(z, z)
        p = z_new + beta * p
        z = z_new

    return x, np.array(res)


def gmres_hist(A, b, restart=None, maxit=200):
    residuals = []

    def callback(rk):
        residuals.append(rk)

    x, info = spla.gmres(
        A, b,
        restart=restart,
        maxiter=maxit,
        callback=callback,
        callback_type="pr_norm",
        rtol=1e-14,
        atol=0.0,
    )

    return x, np.array(residuals) * np.linalg.norm(b)


# ============================================================
# Simple QMR by bi-Lanczos projection
# ============================================================

def qmr_bilanczos(A, b, maxit=200):
    n = b.size
    x0 = np.zeros(n)
    r0 = b - A @ x0

    beta = np.linalg.norm(r0)
    v = r0 / beta
    w = v.copy()

    V = []
    alpha_list = []
    beta_list = []
    gamma_list = []

    v_old = np.zeros(n)
    w_old = np.zeros(n)
    beta_old = 0.0
    gamma_old = 0.0

    res_true = []
    res_quasi = []

    for k in range(maxit):
        V.append(v.copy())

        Av = A @ v
        ATw = A.T @ w

        alpha = np.dot(w, Av)
        alpha_list.append(alpha)

        v_new = Av - alpha * v - gamma_old * v_old
        w_new = ATw - alpha * w - beta_old * w_old

        delta = np.dot(w_new, v_new)
        if abs(delta) < 1e-30:
            break

        beta_k = np.sqrt(abs(delta))
        gamma_k = delta / beta_k

        beta_list.append(beta_k)
        gamma_list.append(gamma_k)

        m = k + 1
        Tbar = np.zeros((m + 1, m))

        for j in range(m):
            Tbar[j, j] = alpha_list[j]
            if j > 0:
                Tbar[j, j - 1] = beta_list[j - 1]
            if j < m - 1:
                Tbar[j, j + 1] = gamma_list[j]
            else:
                Tbar[j + 1, j] = beta_list[j]

        rhs = np.zeros(m + 1)
        rhs[0] = beta

        y, *_ = np.linalg.lstsq(Tbar, rhs, rcond=None)

        Vk = np.column_stack(V)
        xk = Vk @ y

        rq = np.linalg.norm(rhs - Tbar @ y)
        rt = np.linalg.norm(b - A @ xk)

        res_quasi.append(rq)
        res_true.append(rt)

        v_old, w_old = v, w
        v = v_new / beta_k
        w = w_new / gamma_k

        beta_old = beta_k
        gamma_old = gamma_k

    return np.array(res_true), np.array(res_quasi)


# ============================================================
# Part (a): convection-diffusion matrix
# ============================================================

def convection_diffusion_matrix(n=32):
    h = 1.0 / (n + 1)
    N = n * n

    rows, cols, data = [], [], []

    def idx(i, j):
        return i * n + j

    for i in range(n):
        x = (i + 1) * h
        for j in range(n):
            y = (j + 1) * h
            k = idx(i, j)

            center = 4.0 / h**2 - 100.0
            rows.append(k); cols.append(k); data.append(center)

            # x direction
            if i > 0:
                rows.append(k); cols.append(idx(i - 1, j))
                data.append(-1.0 / h**2 - 40.0 * x / (2 * h))
            if i < n - 1:
                rows.append(k); cols.append(idx(i + 1, j))
                data.append(-1.0 / h**2 + 40.0 * x / (2 * h))

            # y direction
            if j > 0:
                rows.append(k); cols.append(idx(i, j - 1))
                data.append(-1.0 / h**2 - 40.0 * y / (2 * h))
            if j < n - 1:
                rows.append(k); cols.append(idx(i, j + 1))
                data.append(-1.0 / h**2 + 40.0 * y / (2 * h))

    return sp.csr_matrix((data, (rows, cols)), shape=(N, N))


def exact_solution_grid(n=32):
    h = 1.0 / (n + 1)
    u = np.zeros(n * n)

    for i in range(n):
        x = (i + 1) * h
        for j in range(n):
            y = (j + 1) * h
            u[i * n + j] = x * (x - 1)**2 * y**2 * (y - 1)**2

    return u


def experiment_a():
    n = 32
    A = convection_diffusion_matrix(n)
    xtrue = exact_solution_grid(n)
    b = A @ xtrue

    Anorm = spla.svds(A, k=1, return_singular_vectors=False)[0]
    denom = Anorm * np.linalg.norm(xtrue)

    _, res_bicg, xnorm_bicg = bicg(A, b, maxit=200, xtrue=xtrue)
    _, res_cgs, xnorm_cgs = cgs(A, b, maxit=200, xtrue=xtrue)

    rel_bicg = res_bicg / denom
    rel_cgs = res_cgs / denom

    plt.figure()
    plt.semilogy(rel_bicg, label="BCG")
    plt.semilogy(rel_cgs, label="CGS")

    k1 = np.argmax(xnorm_bicg)
    k2 = np.argmax(xnorm_cgs)

    plt.scatter(k1, rel_bicg[k1], marker="o")
    plt.scatter(k2, rel_cgs[k2], marker="s")

    plt.xlabel("iteration n")
    plt.ylabel(r"$\|b-Ax_n\| / (\|A\|\|x^*\|)$")
    plt.title("Part (a): BCG vs CGS")
    plt.legend()
    plt.grid(True)
    plt.savefig("part_a_bcg_cgs.png", dpi=200)

    print("Part (a)")
    print("BCG max ||x_n||/||x*|| =", xnorm_bicg[k1], "at n =", k1)
    print("CGS max ||x_n||/||x*|| =", xnorm_cgs[k2], "at n =", k2)
    print("BCG attainable accuracy ≈", np.min(rel_bicg))
    print("CGS attainable accuracy ≈", np.min(rel_cgs))


# ============================================================
# Part (b)(c): random diagonalizable real matrix
# ============================================================

def random_real_nonnormal_matrix(seed=1):
    rng = np.random.default_rng(seed)

    blocks = []

    # real eigenvalues
    blocks.append(np.array([[4.0]]))
    blocks.append(np.array([[0.5]]))
    blocks.append(np.array([[-1.0]]))

    # 50 conjugate pairs
    for _ in range(50):
        a = rng.uniform(1.0, 2.0)
        b = rng.uniform(-1.0, 1.0)
        blocks.append(np.array([[a, b], [-b, a]]))

    D = sp.block_diag(blocks).toarray()

    V = rng.normal(size=(103, 103))
    while np.linalg.cond(V) > 1e5:
        V = rng.normal(size=(103, 103))

    A = V @ D @ np.linalg.inv(V)
    return A


def experiment_b():
    rng = np.random.default_rng(2)

    A = random_real_nonnormal_matrix()
    n = A.shape[0]

    xtrue = rng.normal(size=n)
    b = A @ xtrue

    _, res_bicg, _ = bicg(A, b, maxit=200)
    res_qmr, qres_qmr = qmr_bilanczos(A, b, maxit=200)

    plt.figure()
    plt.semilogy(res_bicg, label="BCG residual")
    plt.semilogy(res_qmr, label="QMR residual")
    plt.semilogy(qres_qmr, "--", label="QMR quasi-residual")

    plt.xlabel("iteration n")
    plt.ylabel("norm")
    plt.title("Part (b): BCG, QMR residuals")
    plt.legend()
    plt.grid(True)
    plt.savefig("part_b_qmr.png", dpi=200)


def experiment_c():
    rng = np.random.default_rng(3)

    A = random_real_nonnormal_matrix()
    n = A.shape[0]

    xtrue = rng.normal(size=n)
    b = A @ xtrue

    maxit = 200

    _, res_gmres = gmres_hist(A, b, restart=None, maxit=maxit)
    _, res_gmres10 = gmres_hist(A, b, restart=10, maxit=maxit)

    _, res_cgs, _ = cgs(A, b, maxit=maxit)
    _, res_bicgstab = bicgstab(A, b, maxit=maxit)
    res_qmr, _ = qmr_bilanczos(A, b, maxit=maxit)
    _, res_cgne = cgne(A, b, maxit=maxit)

    methods = {
        "GMRES": (res_gmres, 1),
        "GMRES(10)": (res_gmres10, 1),
        "CGS": (res_cgs, 2),
        "BCGSTAB": (res_bicgstab, 2),
        "QMR": (res_qmr, 2),
        "CGNE": (res_cgne, 2),
    }

    plt.figure()
    for name, (res, mv_per_iter) in methods.items():
        matvecs = mv_per_iter * np.arange(1, len(res) + 1)
        plt.semilogy(matvecs, res / np.linalg.norm(b), label=name)

    plt.xlabel("number of matrix-vector multiplications")
    plt.ylabel(r"$\|r_n\|/\|b\|$")
    plt.title("Part (c): residual vs matvecs")
    plt.legend()
    plt.grid(True)
    plt.savefig("part_c_matvecs.png", dpi=200)

    plt.figure()
    for name, (res, mv_per_iter) in methods.items():
        matvecs = mv_per_iter * np.arange(1, len(res) + 1)
        flops = 9 * n * matvecs
        plt.semilogy(flops, res / np.linalg.norm(b), label=name)

    plt.xlabel("estimated floating point operations")
    plt.ylabel(r"$\|r_n\|/\|b\|$")
    plt.title("Part (c): residual vs flops")
    plt.legend()
    plt.grid(True)
    plt.savefig("part_c_flops.png", dpi=200)


if __name__ == "__main__":
    experiment_a()
    experiment_b()
    experiment_c()
    plt.show()
