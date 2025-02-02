import itertools as it
from collections import namedtuple
from pathlib import Path
from typing import TypeAlias

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy import Symbol, integrate, pi, sqrt

GSVec = namedtuple(
    "GSVec",
    [
        "f",
        "norm",
        "norm_sq",
        "e",
    ],
)
GSVecs: TypeAlias = dict[int, GSVec]


# Polynomial basis
x = Symbol("x")
v1 = sp.Integer("1")
v2 = x
v3 = x**2
v4 = x**3
v5 = x**4
v6 = x**5

v = sp.sin(x)  # Function we approximate on [-π, π]


def comp_int(f):
    """Compute integral of f on [-pi, pi]"""
    return integrate(f, (x, -pi, pi))


def inp(f, g):
    r"""Compute inner product of f and g

    The inner product is defined as:
      \begin{align*}
          \langle f, g \rangle = \int_{-\pi}^{\pi} fg.
      \end{align*}
    i.e. the integral of the product of the functions from -π to π.
    """
    return comp_int(f * g)


def norm_sq(f):
    return inp(f, f)


def norm(f):
    return sqrt(norm_sq(f))


def gs_step(v, gs_vecs: list[GSVec]):
    """Compute Gram-Schmidt step for vector v in basis.

    gs_vecs are the Gram-Schmidt vectors computed so far
    i.e. vectors for j = 1, 2, ..., k - 1
    when v = v_k.
    """
    return v - sum([inp(v, gs.f) / gs.norm_sq * gs.f for gs in gs_vecs])


def gram_schmidt(basis):
    gs_vecs: GSVecs = dict()
    f1 = basis[1]
    ns_f1 = norm_sq(f1)
    n_f1 = sqrt(ns_f1)
    e1 = f1 / n_f1
    gs_vecs[1] = GSVec(f1, n_f1, ns_f1, e1)
    for k in range(2, 6 + 1):
        v_k = basis[k]
        f_k = gs_step(v_k, list(gs_vecs.values()))
        ns_fk = norm_sq(f_k)
        n_fk = sqrt(ns_fk)
        gs_vecs[k] = GSVec(f_k, n_fk, ns_fk, f_k / n_fk)
    return gs_vecs


def projection_gs(v, gs_vecs: GSVecs):
    r"""Projects v onto the subspace spanned by the orthonormal basis.

    If $e_{1}, \ldots, e_{m}$ be an orthonormal basis of $U$
    then the projection of $v$ onto $U$ is given by
    $$
        P_{U}v =
        \langle v, e_{1} \rangle e_{1}
        + \cdots
        + \langle v, e_{m} \rangle e_{m}
    $$
    """
    return sum([inp(v, gs.e) * gs.e for gs in gs_vecs.values()])


def taylor_sin_5(x):
    return x - x**3 / 6 + x**5 / 120  # 3! = 6;  5! = 120


def plot_projection(proj):
    f_proj_np = sp.lambdify(x, proj, "numpy")

    def plot(xs):

        fs_proj = f_proj_np(xs)
        fs_sin = np.sin(xs)
        fs_taylor = taylor_sin_5(xs)

        # Create plot
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(xs, fs_sin, label=r"$v(x) = \sin(x)$", linestyle="dashed", color="blue")
        ax.plot(xs, fs_proj, label=r"Projection $P_{U}v$", color="red")
        ax.plot(
            xs,
            fs_taylor,
            label=r"$p(x) = x - \frac{x^{3}}{3!} + \frac{x^{5}}{5!}$",
            color="green",
        )

        # Labels and legend
        ax.set_xlabel("x")
        ax.set_ylabel("f(x)")
        ax.set_title("Projection and Taylor approximation of sin(x)")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return fig

    (fp := Path.cwd().joinpath("plots")).mkdir(exist_ok=True)
    plot(xs=np.linspace(-np.pi, np.pi, 100)).savefig(fp / "proj.png")
    plot(xs=np.linspace(2, np.pi, 500)).savefig(fp / "proj2.png")
    plot(xs=np.linspace(3, np.pi, 1000)).savefig(fp / "proj3.png")

    # Relative approximation error in tails
    xs = np.linspace(3, np.pi, 1000)
    fs_proj = f_proj_np(xs)
    fs_sin = np.sin(xs)
    fs_taylor = taylor_sin_5(xs)
    rel_error = (fs_taylor - fs_sin) / (fs_proj - fs_sin)
    print(
        f"Largest relative error for x ∈ [3, π]:"
        f" (x={xs[np.argmax(rel_error)]}, error= {np.max(rel_error):.2f})",
    )


def test_inp_orth(gs_vecs: GSVecs):
    for i, j in list(it.combinations(range(1, 6 + 1), 2)):
        print(
            f"Inner product of {i} and {j} = {inp(gs_vecs[i].e, gs_vecs[j].e).evalf()}"
        )


def main():
    basis = dict(zip(range(1, 6 + 1), [v1, v2, v3, v4, v5, v6]))
    print("Creating orthonormal basis...")
    gs_vecs = gram_schmidt(basis)
    print("Projecting sin(x) onto orthonormal basis...")
    proj = projection_gs(v, gs_vecs)
    proj_a = proj.evalf(strict=True)
    print("Projection: ", proj_a)
    plot_projection(proj_a)
    print("Computing pairwise inner products (might take some time)...")
    test_inp_orth(gs_vecs)
    print(f"Projection exact (latex):\n\t{sp.latex(proj)}")


if __name__ == "__main__":
    main()
