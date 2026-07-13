ETAS: list[float] = [0.5, 1, 10, 30, 100, float("inf")]
CODE_TYPES: list[str] = ["css", "xzzx"]
DISTANCES: list[int] = [3, 5, 7, 9, 11, 13, 15, 17, 19]

def build_window(center: float) -> float:
    half_width = 0.003
    return (center - half_width, center + half_width)

P_STEP: float = 0.0001
P_WINDOWS: dict[tuple[float, str], tuple[float, float]] = {
    (0.5, "css"): build_window(0.007),
    (0.5, "xzzx"): build_window(0.007),
    (1, "css"): build_window(0.0065),
    (1, "xzzx"): build_window(0.0078),
    (10, "css"): build_window(0.006),
    (10, "xzzx"): build_window(0.013),
    (30, "css"): build_window(0.006),
    (30, "xzzx"): build_window(0.016),
    (100, "css"): build_window(0.006),
    (100, "xzzx"): build_window(0.019),
    (float("inf"), "css"): build_window(0.006),
    (float("inf"), "xzzx"): build_window(0.022),
}
