"""Companion calculations for "Multiperiod Attribution: Who Needs Linking Functions?"""

import math

CASE = "base"  # "base" or "pure_selection"
START_WEALTH = 100.0
TOL = 1e-12
SECTORS = [[0, 1], [2, 3]]
q = lambda x: round(x, 10)

W_BENCH = [
    [0.20, 0.10, 0.10, 0.60],
    [0.20, 0.20, 0.40, 0.20],
    [0.15, 0.05, 0.10, 0.70],
    [0.10, 0.10, 0.20, 0.60],
]

W_PORT_BASE = [
    [0.05, 0.05, 0.20, 0.70],
    [0.05, 0.05, 0.40, 0.50],
    [0.20, 0.10, 0.20, 0.50],
    [0.20, 0.20, 0.20, 0.40],
]

W_PORT_PURE_SEL = [
    [0.10, 0.00, 0.70, 0.20],
    [0.10, 0.00, 0.75, 0.15],
    [0.30, 0.00, 0.30, 0.40],
    [0.40, 0.00, 0.25, 0.35],
]

W_BENCH_PURE_SEL = W_PORT_BASE

RETURNS = [
    [0.05, 0.05, 0.05, 0.10],
    [-0.10, -0.05, 0.02, 0.08],
    [-0.05, 0.06, -0.02, 0.03],
    [0.05, 0.10, 0.04, 0.08],
]


def period_stats(wp, wb, r):
    ap = [sum(wp[i] for i in s) for s in SECTORS]
    ab = [sum(wb[i] for i in s) for s in SECTORS]
    rp = [sum(wp[i] * r[i] for i in s) / ap[j] for j, s in enumerate(SECTORS)]
    rb = [sum(wb[i] * r[i] for i in s) / ab[j] for j, s in enumerate(SECTORS)]

    RP = sum(a * x for a, x in zip(ap, rp))
    RB = sum(a * x for a, x in zip(ab, rb))
    A = sum((ap[j] - ab[j]) * (rb[j] - RB) for j in range(len(SECTORS)))
    S = sum(ab[j] * (rp[j] - rb[j]) for j in range(len(SECTORS)))
    I = sum((ap[j] - ab[j]) * (rp[j] - rb[j]) for j in range(len(SECTORS)))
    RPB = sum(ap[j] * rb[j] for j in range(len(SECTORS)))
    RBP = sum(ab[j] * rp[j] for j in range(len(SECTORS)))
    return RP, RB, (A, S, I), RPB, RBP


def run(w_port, w_bench):
    n = len(RETURNS)
    stats = [period_stats(w_port[t], w_bench[t], RETURNS[t]) for t in range(n)]
    RP = [x[0] for x in stats]
    RB = [x[1] for x in stats]
    G = [x[2] for x in stats]
    RA = [rp - rb for rp, rb in zip(RP, RB)]

    for t in range(n):
        assert abs(sum(G[t]) - RA[t]) < TOL

    WP = WB = START_WEALTH
    WA = 0.0
    ewa = [0.0, 0.0, 0.0]
    carry = 0.0

    print("EWA account across time (per 100 initial):")
    print(f"{'t':>2} {'A':>9} {'S':>9} {'I':>9} {'Carry':>9} {'sum':>9} {'dW^A':>9}")
    for t in range(n):
        E = [WP * g for g in G[t]]
        C = WA * RB[t]
        dWA = WP * RP[t] - WB * RB[t]
        assert abs(sum(E) + C - dWA) < TOL

        ewa = [x + y for x, y in zip(ewa, E)]
        carry += C
        print(f"{t+1:>2} {q(E[0]):>+9.4f} {q(E[1]):>+9.4f} {q(E[2]):>+9.4f} "
              f"{q(C):>+9.4f} {q(sum(E)+C):>+9.4f} {q(dWA):>+9.4f}")

        WA += dWA
        WP *= 1 + RP[t]
        WB *= 1 + RB[t]

    print("-" * 62)
    print(f"{'':>2} {q(ewa[0]):>+9.4f} {q(ewa[1]):>+9.4f} {q(ewa[2]):>+9.4f} "
          f"{q(carry):>+9.4f} {q(sum(ewa)+carry):>+9.4f} {q(WA):>+9.4f}")

    F = [0.0, 0.0, 0.0]
    WP = START_WEALTH
    for t in range(n):
        F = [(1 + RB[t]) * f + WP * g for f, g in zip(F, G[t])]
        WP *= 1 + RP[t]

    cPP = cBB = cPB = cBP = 0.0
    for t in range(n):
        cPP = (1 + RP[t]) * (1 + cPP) - 1
        cBB = (1 + RB[t]) * (1 + cBB) - 1
        cPB = (1 + stats[t][3]) * (1 + cPB) - 1
        cBP = (1 + stats[t][4]) * (1 + cBP) - 1
    berg = [
        START_WEALTH * (cPB - cBB),
        START_WEALTH * (cBP - cBB),
        START_WEALTH * (cPP - cPB - cBP + cBB),
    ]

    RPc = math.prod(1 + x for x in RP) - 1
    RBc = math.prod(1 + x for x in RB) - 1
    K = (math.log1p(RPc) - math.log1p(RBc)) / (RPc - RBc)

    def kt(t):
        if abs(RA[t]) < 1e-15:
            return 1 / (1 + RP[t])
        return (math.log1p(RP[t]) - math.log1p(RB[t])) / RA[t]

    carino = [
        START_WEALTH * sum(kt(t) / K * G[t][k] for t in range(n))
        for k in range(3)
    ]

    def rev_grap_multiplier(t):
        return (
            math.prod(1 + RP[s] for s in range(t))
            * math.prod(1 + RB[u] for u in range(t + 1, n))
        )

    rgrap = [
        START_WEALTH * sum(G[t][k] * rev_grap_multiplier(t) for t in range(n))
        for k in range(3)
    ]

    rows = [
        ("EWA (explicit Carry)", ewa, carry),
        ("Frongello (2002)", F, None),
        ("Berg (2014)", berg, None),
        ("Carino (1999)", carino, None),
        ("Reverse GRAP (1997)", rgrap, None),
    ]

    print("\nFive methods on the same data (terminal, per 100 initial):")
    print(f"{'Method':<22}{'A':>9}{'S':>9}{'I':>9}{'Carry':>9}{'Total':>9}")
    for name, effects, C in rows:
        total = sum(effects) + (0.0 if C is None else C)
        assert abs(total - WA) < 1e-9
        ctext = "—" if C is None else f"{q(C):+.4f}"
        print(f"{name:<22}{q(effects[0]):>+9.4f}{q(effects[1]):>+9.4f}{q(effects[2]):>+9.4f}"
              f"{ctext:>9}{q(total):>+9.4f}")

    D_path = []
    WP = WB = START_WEALTH
    for t in range(n):
        WP *= 1 + RP[t]
        WB *= 1 + RB[t]
        D_path.append(WP - WB)

    ewa_c, fro_c, berg_c, car_c, rg_c = [], [], [], [], []
    WP = START_WEALTH
    WAc = Fc = cPP = cBB = car_run = rg_run = 0.0
    for t in range(n):
        WAc += WP * RA[t] + WAc * RB[t]
        Fc = (1 + RB[t]) * Fc + WP * RA[t]
        cPP = (1 + RP[t]) * (1 + cPP) - 1
        cBB = (1 + RB[t]) * (1 + cBB) - 1
        car_run += kt(t) / K * RA[t]
        rg_run += RA[t] * rev_grap_multiplier(t)

        ewa_c.append(WAc)
        fro_c.append(Fc)
        berg_c.append(START_WEALTH * (cPP - cBB))
        car_c.append(START_WEALTH * car_run)
        rg_c.append(START_WEALTH * rg_run)
        WP *= 1 + RP[t]

    print("\nCumulative attribution through time vs realised active wealth")
    print(f"(horizon fixed at T = {n}; per 100 initial):")
    print(f"{'t':>2} {'W^A(t)':>9} {'EWA':>9} {'Frongello':>10} {'Berg':>9} {'Carino':>9} {'RevGRAP':>9}")
    for t in range(n):
        print(f"{t+1:>2} {q(D_path[t]):>+9.4f} {q(ewa_c[t]):>+9.4f} {q(fro_c[t]):>+10.4f} "
              f"{q(berg_c[t]):>+9.4f} {q(car_c[t]):>+9.4f} {q(rg_c[t]):>+9.4f}")


if CASE == "base":
    W_PORT, W_B = W_PORT_BASE, W_BENCH
elif CASE == "pure_selection":
    W_PORT, W_B = W_PORT_PURE_SEL, W_BENCH_PURE_SEL
else:
    raise ValueError("CASE must be 'base' or 'pure_selection'")

print(f"Explicit Wealth Attribution — case: {CASE}\n")
run(W_PORT, W_B)
