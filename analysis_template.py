"""
3iExit — Heliospheric Amplification Test (starter)
Inputs (expected in repo root or /mnt/data):
  - area_angle_timeseries.csv  (from geometry pipeline)
  - solar_wind.csv   (time, Vsw_kms, B_tot_nT, Bz_nT, n_p_cc)
  - geomag.csv       (time, H_nT)  # or |B|
  - elf.csv          (time, band_7_8, band_14, band_20_40)
Outputs:
  - predictions_summary.csv
Notes:
  - Replace ELF proxy with raw PSD workflow when available.
"""
import os, numpy as np, pandas as pd
from scipy.signal import welch
from scipy.stats import linregress

DATA_DIR = "."  # adjust if needed

def load_csv(name):
    path = os.path.join(DATA_DIR, name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path, parse_dates=["time"])
    df["time"] = pd.to_datetime(df["time"], utc=True)
    return df.sort_values("time")

def high_k_slope(signal, fs_hz, band_hz):
    f, Pxx = welch(signal, fs=fs_hz, nperseg=min(len(signal)//2, 2048))
    m = (f >= band_hz[0]) & (f <= band_hz[1]) & np.isfinite(Pxx) & (Pxx>0) & (f>0)
    if m.sum() < 8: return np.nan, np.nan
    x = np.log10(f[m]); y = np.log10(Pxx[m])
    s, i, r, p, se = linregress(x, y)
    return s, r**2

def moving_beta(ts, col, fs_hz, band_hz, window_h=1, step_h=1):
    out=[]; t=ts["time"].to_numpy(); v=ts[col].to_numpy()
    if len(t)<2: return pd.DataFrame(columns=["time_center","beta","r2"])
    dt = np.median(np.diff(t).astype('timedelta64[s]').astype(float))
    sph = max(1, int(round(3600.0/dt)))
    win = max(sph*window_h, 8); step=max(sph*step_h,1)
    for i in range(0, len(v)-win+1, step):
        seg=v[i:i+win]; beta,r2=high_k_slope(seg, fs_hz, band_hz)
        out.append((t[i+win//2], beta, r2))
    return pd.DataFrame(out, columns=["time_center","beta","r2"])

def closure_residual(phi, rhoI):
    dphi = np.gradient(phi); term = rhoI*dphi; return np.gradient(term)

def main():
    geom = load_csv("area_angle_timeseries.csv")
    s = geom["A_scaled"].rolling(3, center=True, min_periods=1).mean()
    mins = [i for i in range(1,len(s)-1) if s.iloc[i]<=s.iloc[i-1] and s.iloc[i]<=s.iloc[i+1]]
    minima = geom.iloc[mins][["time","A_scaled","theta_deg"]].rename(columns={"time":"t_min"})

    sw = load_csv("solar_wind.csv")
    gm = load_csv("geomag.csv")
    elf = load_csv("elf.csv")

    def windowed(df, t0, span_h):
        t0 = pd.to_datetime(t0, utc=True)
        return df[(df["time"]>=t0-pd.Timedelta(hours=span_h)) & (df["time"]<=t0+pd.Timedelta(hours=span_h))]

    def forcing_mask(sw_win):
        sw1 = sw_win.copy().set_index("time").resample("1H").mean().reset_index()
        c1 = sw1["B_tot_nT"]>=8; c2 = sw1["Vsw_kms"]>=500; c3 = sw1["Bz_nT"]<=-5
        gate = (c1|c2|c3).rolling(3, min_periods=1).max().astype(bool)
        sw1["forcing_gate"]=gate.values; return sw1

    WINDOW_H=48; results=[]
    for _, row in minima.iterrows():
        tmin=row["t_min"]
        sw_win=windowed(sw,tmin,WINDOW_H)
        gm_win=windowed(gm,tmin,WINDOW_H)
        elf_win=windowed(elf,tmin,WINDOW_H)

        sw_gate=forcing_mask(sw_win)
        forced=sw_gate["forcing_gate"].any()

        # geomag slopes: assume 1 min cadence -> fs=1/60 Hz; band 0.03–0.18 Hz
        gm1 = gm_win.copy().set_index("time").resample("1T").interpolate().reset_index()
        gm_beta = moving_beta(gm1,"H_nT", fs_hz=1/60.0, band_hz=(0.03,0.18))

        # closure residual: Phi = low-freq geomag; rhoI = z(Vsw,|B|,n_p)
        gm_low = gm1["H_nT"].rolling(60, center=True, min_periods=1).mean()
        swz = sw_gate.copy().set_index("time").reindex(gm1["time"]).interpolate().reset_index()
        for c in ["Vsw_kms","B_tot_nT","n_p_cc"]:
            swz[c+"_z"] = (swz[c]-swz[c].mean())/swz[c].std(ddof=0)
        swz["rhoI"]=swz[["Vsw_kms_z","B_tot_nT_z","n_p_cc_z"]].mean(axis=1)
        R = closure_residual(gm_low.to_numpy(), swz["rhoI"].to_numpy())

        beta_hit = (gm_beta["beta"]<=-3.7).sum()>=3
        R_peak = np.nanpercentile(np.abs(R),90)
        R_hit = (np.nanmax(np.abs(R))>=R_peak)

        results.append({
            "t_min": str(tmin),
            "forcing_present": bool(forced),
            "beta_geomag_hit": bool(beta_hit),
            "R_residual_peakish": bool(R_hit),
            "note": "ELF PSD not computed in this starter; replace with raw ELF PSD when available."
        })

    pd.DataFrame(results).to_csv("predictions_summary.csv", index=False)
    print("Wrote predictions_summary.csv")

if __name__=="__main__":
    import pandas as pd
    main()
