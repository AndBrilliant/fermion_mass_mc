"""
engine.running — MSbar mass running with flavor-threshold matching.

Primary path: four-loop CRunDec via the rundec Python wrapper.  This
reproduces the paper's m_s^MSbar(mu_star) = 95.075 MeV reference to 4
significant figures, with propagated sigma 0.6919 MeV matching the
paper's 0.69 MeV.

Fallback path: pure-Python one-loop with threshold matching.  Used
only if `rundec` is not importable.  ~0.7% less accurate than four
loops at scales a few GeV from the reference (enough to shift hit
fractions measurably), so rundec should be preferred for production.
"""
import math

# --- rundec availability --------------------------------------------
try:
    import rundec
    _CRD = rundec.CRunDec()
    HAS_RUNDEC = True
except ImportError:
    _CRD = None
    HAS_RUNDEC = False

# --- Constants ------------------------------------------------------
DEFAULT_MC = 1.273       # GeV, m_c(m_c)
DEFAULT_MB = 4.183       # GeV, m_b(m_b)
DEFAULT_MT = 173.0       # GeV, m_t pole (approximate for nf=6 matching)
ALPHA_S_MZ = 0.1180
M_Z = 91.1876
NLOOPS = 4               # rundec loop order
GAMMA_0 = 8.0            # one-loop mass anomalous dimension

# --- Threshold alpha_s cache (rundec path) --------------------------
_threshold_cache = {}


def _get_threshold_alphas(mc, mb, mt, as_mz):
    """Return a dict of alpha_s values at thresholds, cached."""
    key = (mc, mb, mt, as_mz)
    if key in _threshold_cache:
        return _threshold_cache[key]
    if not HAS_RUNDEC:
        _threshold_cache[key] = None
        return None
    as_mb_5 = _CRD.AlphasExact(as_mz, M_Z, mb, 5, NLOOPS)
    as_mb_4 = _CRD.DecAsDownMS(as_mb_5, mb, mb, 4, NLOOPS)
    as_mc_4 = _CRD.AlphasExact(as_mb_4, mb, mc, 4, NLOOPS)
    as_mc_3 = _CRD.DecAsDownMS(as_mc_4, mc, mc, 3, NLOOPS)
    as_mt_5 = _CRD.AlphasExact(as_mz, M_Z, mt, 5, NLOOPS)
    as_mt_6 = _CRD.DecAsUpMS(as_mt_5, mt, mt, 5, NLOOPS)
    d = {
        "mc_3": as_mc_3, "mc_4": as_mc_4,
        "mb_4": as_mb_4, "mb_5": as_mb_5,
        "mt_5": as_mt_5, "mt_6": as_mt_6,
    }
    _threshold_cache[key] = d
    return d


def _nf_of(mu, mc, mb, mt):
    if mu < mc: return 3
    if mu < mb: return 4
    if mu < mt: return 5
    return 6


def alpha_s_at(mu, mc=DEFAULT_MC, mb=DEFAULT_MB, mt=DEFAULT_MT,
               as_mz=ALPHA_S_MZ):
    """alpha_s at scale mu with proper threshold matching."""
    if mu <= 0:
        return float("nan")
    if HAS_RUNDEC:
        thr = _get_threshold_alphas(mc, mb, mt, as_mz)
        if mu >= mt:
            return _CRD.AlphasExact(thr["mt_6"], mt, mu, 6, NLOOPS)
        elif mu >= mb:
            return _CRD.AlphasExact(as_mz, M_Z, mu, 5, NLOOPS)
        elif mu >= mc:
            return _CRD.AlphasExact(thr["mb_4"], mb, mu, 4, NLOOPS)
        else:
            return _CRD.AlphasExact(thr["mc_3"], mc, mu, 3, NLOOPS)
    return _alpha_s_1loop_fallback(mu, mc, mb, mt, as_mz)


def run_mass(m_ref, mu_ref, mu_target, mc=DEFAULT_MC, mb=DEFAULT_MB,
             mt=DEFAULT_MT, as_mz=ALPHA_S_MZ):
    """Run an MSbar mass from mu_ref to mu_target with threshold matching.

    For light quarks (like m_s), mass matching across thresholds is
    continuous to leading order.  This function uses continuous mass
    matching and relies on rundec's properly-matched alpha_s values.
    """
    if m_ref <= 0 or mu_ref <= 0 or mu_target <= 0:
        return float("nan")
    if mu_ref == mu_target:
        return m_ref
    if not HAS_RUNDEC:
        return _run_mass_1loop_fallback(m_ref, mu_ref, mu_target,
                                        mc, mb, mt, as_mz)
    thr = _get_threshold_alphas(mc, mb, mt, as_mz)
    nf_ref = _nf_of(mu_ref, mc, mb, mt)
    nf_tgt = _nf_of(mu_target, mc, mb, mt)
    if nf_ref == nf_tgt:
        as_r = alpha_s_at(mu_ref, mc, mb, mt, as_mz)
        as_t = alpha_s_at(mu_target, mc, mb, mt, as_mz)
        return _CRD.mMS2mMS(m_ref, as_r, as_t, nf_ref, NLOOPS)
    # Multi-region walk (continuous mass matching at thresholds)
    m = m_ref
    mu = mu_ref
    nf = nf_ref
    while nf != nf_tgt:
        if nf < nf_tgt:
            # Going up: next threshold is mb (nf 4 -> 5) or mt (5 -> 6)
            if nf == 3:
                nxt_mu, nxt_as, nxt_nf = mc, thr["mc_4"], 4
                cur_as_at_nxt = thr["mc_3"]
            elif nf == 4:
                nxt_mu, nxt_as, nxt_nf = mb, thr["mb_5"], 5
                cur_as_at_nxt = thr["mb_4"]
            elif nf == 5:
                nxt_mu, nxt_as, nxt_nf = mt, thr["mt_6"], 6
                cur_as_at_nxt = thr["mt_5"]
            else:
                return float("nan")
            # Run mass to threshold in current nf, then "match" (continuous)
            as_cur = alpha_s_at(mu, mc, mb, mt, as_mz)
            m = _CRD.mMS2mMS(m, as_cur, cur_as_at_nxt, nf, NLOOPS)
            mu, nf = nxt_mu, nxt_nf
        else:
            # Going down: next threshold is mb (5 -> 4), mc (4 -> 3)
            if nf == 6:
                nxt_mu, nxt_as, nxt_nf = mt, thr["mt_5"], 5
                cur_as_at_nxt = thr["mt_6"]
            elif nf == 5:
                nxt_mu, nxt_as, nxt_nf = mb, thr["mb_4"], 4
                cur_as_at_nxt = thr["mb_5"]
            elif nf == 4:
                nxt_mu, nxt_as, nxt_nf = mc, thr["mc_3"], 3
                cur_as_at_nxt = thr["mc_4"]
            else:
                return float("nan")
            as_cur = alpha_s_at(mu, mc, mb, mt, as_mz)
            m = _CRD.mMS2mMS(m, as_cur, cur_as_at_nxt, nf, NLOOPS)
            mu, nf = nxt_mu, nxt_nf
    # Final leg within target nf
    as_mu = alpha_s_at(mu, mc, mb, mt, as_mz)
    as_tgt = alpha_s_at(mu_target, mc, mb, mt, as_mz)
    return _CRD.mMS2mMS(m, as_mu, as_tgt, nf_tgt, NLOOPS)


def m_s_at(mu_target, m_s_2gev=0.09344, sigma_m_s_2gev=0.00068,
           mc=DEFAULT_MC, mb=DEFAULT_MB, mt=DEFAULT_MT,
           as_mz=ALPHA_S_MZ):
    """Run m_s from 2 GeV to mu_target.  Returns (m_s, sigma)."""
    m = run_mass(m_s_2gev, 2.0, mu_target, mc, mb, mt, as_mz)
    if not math.isfinite(m) or m_s_2gev <= 0:
        return (float("nan"), float("nan"))
    sigma = sigma_m_s_2gev * (m / m_s_2gev)
    return m, sigma


# --- Fallback 1-loop functions (used only if rundec unavailable) ---

def _beta_0(nf):
    return (33.0 - 2.0 * nf) / 3.0


def _alpha_s_1loop_fallback(mu, mc, mb, mt, as_mz):
    """One-loop alpha_s with threshold matching.  Fallback only."""
    if mu <= 0:
        return float("nan")
    nf = 5
    as_cur = as_mz
    mu_cur = M_Z
    if mu < mb and mu_cur >= mb:
        inv = 1.0/as_cur + (_beta_0(5)/(2*math.pi)) * math.log(mb/mu_cur)
        as_cur = 1.0/inv if inv > 0 else float("nan")
        mu_cur = mb
        nf = 4
    if mu < mc and mu_cur >= mc:
        inv = 1.0/as_cur + (_beta_0(4)/(2*math.pi)) * math.log(mc/mu_cur)
        as_cur = 1.0/inv if inv > 0 else float("nan")
        mu_cur = mc
        nf = 3
    if not math.isfinite(as_cur) or as_cur <= 0:
        return float("nan")
    inv = 1.0/as_cur + (_beta_0(nf)/(2*math.pi)) * math.log(mu/mu_cur)
    return 1.0/inv if inv > 0 else float("nan")


def _run_mass_1loop_fallback(m_ref, mu_ref, mu_target, mc, mb, mt, as_mz):
    """One-loop mass running with continuous threshold matching."""
    if m_ref <= 0 or mu_ref <= 0 or mu_target <= 0:
        return float("nan")
    if mu_ref == mu_target:
        return m_ref
    thresholds = sorted([mc, mb, mt])
    lo, hi = min(mu_ref, mu_target), max(mu_ref, mu_target)
    interior = [th for th in thresholds if lo < th < hi]
    path = [mu_ref] + (interior if mu_target > mu_ref else interior[::-1]) + [mu_target]
    m_cur = m_ref
    mu_cur = mu_ref
    for mu_nxt in path[1:]:
        as_cur = _alpha_s_1loop_fallback(mu_cur, mc, mb, mt, as_mz)
        as_nxt = _alpha_s_1loop_fallback(mu_nxt, mc, mb, mt, as_mz)
        if not math.isfinite(as_cur) or not math.isfinite(as_nxt):
            return float("nan")
        nf_mid = _nf_of(math.sqrt(mu_cur*mu_nxt), mc, mb, mt)
        exp = GAMMA_0 / (2.0 * _beta_0(nf_mid))
        m_cur = m_cur * (as_nxt/as_cur) ** exp
        mu_cur = mu_nxt
    return m_cur


if __name__ == "__main__":
    # Self-tests
    print(f"HAS_RUNDEC = {HAS_RUNDEC}")
    mu_star = 1.88310
    m, s = m_s_at(mu_star)
    print(f"m_s(2 GeV) = 93.44 +/- 0.68 MeV")
    print(f"m_s({mu_star:.4f} GeV) = {m*1000:.4f} +/- {s*1000:.4f} MeV")
    print(f"Paper:                   95.07 +/- 0.69 MeV")
    print(f"alpha_s(m_tau) = {alpha_s_at(1.77693):.5f}  (PDG 4-loop: ~0.33)")
    # Threshold crossing sanity check
    m_at_10 = run_mass(0.09344, 2.0, 10.0)
    print(f"m_s(10 GeV, crosses m_b) = {m_at_10*1000:.4f} MeV")
    m_at_200 = run_mass(0.09344, 2.0, 200.0)
    print(f"m_s(200 GeV, crosses m_b, m_t) = {m_at_200*1000:.4f} MeV")
    m_at_1 = run_mass(0.09344, 2.0, 1.0)
    print(f"m_s(1 GeV, crosses m_c) = {m_at_1*1000:.4f} MeV")
