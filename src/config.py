ETAS: list[float] = [0.5, 1, 10, 30, 100, float("inf")]
CODE_TYPES: list[str] = ["css", "xzzx"]
DISTANCES: list[int] = [3, 5, 7, 9, 11, 13, 15]

P_STEP: float = 0.001
P_WINDOWS: dict[tuple[float, str], tuple[float, float]] = {
    (0.5, "css"): (0.001, 0.02),
    (0.5, "xzzx"): (0.001, 0.02),
    (1, "css"): (0.001, 0.02),
    (1, "xzzx"): (0.001, 0.02),
    (10, "css"): (0.001, 0.02),
    (10, "xzzx"): (0.005, 0.025),
    (30, "css"): (0.001, 0.02),
    (30, "xzzx"): (0.007, 0.03),
    (100, "css"): (0.001, 0.02),
    (100, "xzzx"): (0.007, 0.03),
    (float("inf"), "css"): (0.001, 0.02),
    (float("inf"), "xzzx"): (0.012, 0.035),
}
