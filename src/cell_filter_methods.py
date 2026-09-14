import numpy as np

# Restates a normal's 10-90 range on the scale stats::mad()'s default constant
# (1.4826) puts MAD on, so the fallback is commensurable with the MAD it replaces.
# Same constant and same three branches as qc_cell_spread() in
# SQANTI-sc_multisample_report.R -- one convention across the tool, not two.
P10_P90_TO_MAD = 2.563


def robust_spread(values):
    """Robust scale of one sample's cell values, with the zero-inflation fallback.

    Returns NaN rather than 0 when even the 10-90 range is empty: a criterion that
    cannot be scaled must be skipped, because median + n*0 flags every cell whose
    value is not exactly the median.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return np.nan
    centre = np.median(v)
    s = np.median(np.abs(v - centre)) * 1.4826
    if not np.isfinite(s) or s <= 0:
        q10, q90 = np.quantile(v, [0.10, 0.90])
        s = (q90 - q10) / P10_P90_TO_MAD
    if not np.isfinite(s) or s <= 0:
        return np.nan
    return s


def method_none(summary, params):
    """No automatic criterion. Rules and the manual lists decide everything."""
    return {}


# The seam the automatic mode plugs into. Each entry takes the (already
# sentinel-stripped) cell summary and returns {CB: reason} for the cells it
# discards; anything absent from the mapping is left to the rules.
AUTO_METHODS = {
    'none': method_none,
}
