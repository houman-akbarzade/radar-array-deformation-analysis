"""Public-safe, parameterized cylindrical-sector array simulation.

Synthetic parameters only. Not a validated reproduction of any operational radar.
Run: python src/python/run_demo.py
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal.windows import taylor

C0 = 299_792_458.0

@dataclass(frozen=True)
class Config:
    frequency_hz: float = 5.8e9
    frequency_min_hz: float = 5.2e9
    frequency_max_hz: float = 5.8e9
    phase_center_power: float = 1.0
    bending_shape: str = 'uniform'
    mechanical_model: str = 'base'
    wind_az_deg: float = 0.0
    height_gradient: float = 0.0
    total_column_count: int = 60
    panel_count: int = 20
    row_count: int = 16
    az_spacing_wavelength: float = 0.52
    vertical_spacing_wavelength: float = 0.48
    seed: int = 2026
    trials: int = 100
    panel_correlation: float = 0.35
    within_panel_correlation: float = 0.995
    displacement_std_wavelength: float = 0.02
    mean_displacement_wavelength: float = 0.0
    radial_coupling: float = 0.0001
    vertical_coupling: float = 0.02
    steer_az_deg: float = 0.0
    target_az_deg: float = 1.5
    beam_el_deg: float = 3.0
    analog_el_deg: float = 0.0
    calibration_half_span_deg: float = 5.0
    peak_half_span_deg: float = 6.0

def geometry(c):
    """Equal arc spacing, with radius inferred from wavelength and full circumference."""
    assert c.total_column_count >= c.panel_count >= 4 and c.panel_count % 2 == 0 and c.row_count >= 3
    lam = C0 / c.frequency_hz
    arc_spacing = c.az_spacing_wavelength * lam
    radius = c.total_column_count * arc_spacing / (2 * np.pi)
    dz = c.vertical_spacing_wavelength * lam
    phi_panels = (np.arange(c.panel_count) - (c.panel_count-1)/2) * (2*np.pi/c.total_column_count)
    panel = np.tile(np.arange(c.panel_count), c.row_count)
    row = np.repeat(np.arange(c.row_count), c.panel_count)
    phi = phi_panels[panel]
    z = (row - (c.row_count-1)/2) * dz
    xyz = np.column_stack((radius*np.cos(phi), radius*np.sin(phi), z))
    return dict(xyz=xyz, phi=phi, panel=panel, row=row, radius=radius, dz=dz, lam=lam)

def displacement(c, g, rng):
    """Open-sector exponential panel correlation and local residual."""
    n = c.panel_count
    K = c.panel_correlation ** np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    # Cholesky supports 0 <= rho < 1; use a spectral square root for rho=1.
    eigval, eigvec = np.linalg.eigh(K)
    shared = eigvec @ (np.sqrt(np.maximum(eigval, 0)) * rng.standard_normal(n))
    q = np.sqrt(c.within_panel_correlation)*shared[g['panel']] + np.sqrt(1-c.within_panel_correlation)*rng.standard_normal(len(g['panel']))
    xi = ((c.frequency_hz-c.frequency_min_hz)/(c.frequency_max_hz-c.frequency_min_hz))**c.phase_center_power
    if not 0 <= xi <= 1: raise ValueError('Frequency outside modeled band')
    if c.bending_shape == 'uniform':
        gain = xi**2*(6-4*xi+xi**2)/3
    elif c.bending_shape == 'tip_force':
        gain = xi**2*(3-xi)/2
    else: raise ValueError('Unknown bending shape')
    scale = np.ones(len(q))
    if c.mechanical_model == 'wind_parametric':
        scale = np.sin(np.deg2rad(c.wind_az_deg)-g['phi'])*(1+c.height_gradient*g['xyz'][:,2]/max(abs(g['xyz'][:,2])))
    elif c.mechanical_model != 'base': raise ValueError('Unknown mechanical model')
    # Tip statistics use the upper-band wavelength; effective motion follows bending gain.
    s = gain*scale*(C0/c.frequency_max_hz)*(c.mean_displacement_wavelength + c.displacement_std_wavelength*q)
    r = c.radial_coupling*s
    z = c.vertical_coupling*s
    phi = g['phi']
    delta = np.column_stack((-np.sin(phi)*s+np.cos(phi)*r,
                             np.cos(phi)*s+np.sin(phi)*r, z))
    return delta, s

def steering(c, xyz, az_deg, el_deg):
    az = np.atleast_1d(np.asarray(az_deg, float))
    el = np.broadcast_to(np.asarray(el_deg, float), az.shape)
    a, e = np.deg2rad(az), np.deg2rad(el)
    u = np.vstack((np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)))
    return np.exp(1j*(2*np.pi*c.frequency_hz/C0)*(xyz @ u))

def weights(c, g, architecture):
    """Row-wise left/right combining, then digital elevation weighting."""
    if architecture not in ('steered_analog', 'fixed_analog'):
        raise ValueError('Unknown architecture')
    n, m = c.panel_count, c.row_count
    azwin = np.hamming(n)
    elwin = taylor(m, nbar=min(4,m), sll=30, norm=False)
    elwin /= elwin.sum()
    analog_el = c.beam_el_deg if architecture == 'steered_analog' else c.analog_el_deg
    k = 2*np.pi/c.lam
    horizontal = g['xyz'][:,0]*np.cos(np.deg2rad(c.steer_az_deg))+g['xyz'][:,1]*np.sin(np.deg2rad(c.steer_az_deg))
    analog = azwin[g['panel']]*np.exp(-1j*k*np.cos(np.deg2rad(analog_el))*horizontal)
    zrow = (np.arange(m)-(m-1)/2)*g['dz']
    digital = elwin*np.exp(-1j*k*zrow*np.sin(np.deg2rad(c.beam_el_deg)))/azwin.sum()
    combined = analog*digital[g['row']]
    sign = np.where(g['panel'] < n//2, -1., 1.)
    return combined, combined*sign

def response(c, xyz, weights_pair, az_deg, el_deg):
    a = steering(c, xyz, az_deg, el_deg)
    return weights_pair[0] @ a, weights_pair[1] @ a

def calibration(c, nominal_xyz, w):
    az = c.steer_az_deg + np.linspace(-c.calibration_half_span_deg, c.calibration_half_span_deg, 1001)
    s, d = response(c, nominal_xyz, w, az, c.beam_el_deg)
    ratio = np.imag(d/s)
    valid = (np.abs(s) > .1*np.max(np.abs(s))) & np.isfinite(ratio)
    center = np.argmin(np.abs(az-c.steer_az_deg))
    slope = np.sign(ratio[center+1]-ratio[center])
    lo = hi = center
    while lo>0 and valid[lo-1] and slope*(ratio[lo]-ratio[lo-1])>1e-12: lo-=1
    while hi<len(az)-1 and valid[hi+1] and slope*(ratio[hi+1]-ratio[hi])>1e-12: hi+=1
    if hi-lo < 10: raise RuntimeError('No monotonic calibration region')
    order = np.argsort(ratio[lo:hi+1])
    return az[lo:hi+1][order], ratio[lo:hi+1][order], .1*np.max(np.abs(s))

def estimate(c, xyz, w, cal):
    s, d = response(c, xyz, w, c.target_az_deg, c.beam_el_deg)
    v = float(np.imag(d[0]/s[0]))
    az, ratio, min_sum = cal
    if abs(s[0]) <= min_sum or not np.isfinite(v) or v<ratio[0] or v>ratio[-1]:
        return np.nan
    return float(np.interp(v, ratio, az))

def peak(c, xyz, w):
    grid = c.steer_az_deg + np.arange(-c.peak_half_span_deg, c.peak_half_span_deg+.001, .25)
    s, _ = response(c, xyz, w, grid, c.beam_el_deg)
    idx = int(np.argmax(np.abs(s)**2))
    if idx in (0, len(grid)-1): return np.nan, np.nan
    def objective(a):
        ss, _ = response(c, xyz, w, a, c.beam_el_deg)
        return -float(abs(ss[0])**2)
    opt = minimize_scalar(objective, bounds=(grid[idx-1],grid[idx+1]), method='bounded',
                          options={'xatol':1e-8})
    return float(opt.x), -float(opt.fun)

def run(c=Config()):
    assert c.frequency_min_hz < c.frequency_max_hz and c.frequency_min_hz <= c.frequency_hz <= c.frequency_max_hz and c.trials>=1
    assert 0<=c.panel_correlation<=1 and 0<=c.within_panel_correlation<=1
    g = geometry(c)
    nominal = g['xyz']
    wa = weights(c,g,'steered_analog')
    wb = weights(c,g,'fixed_analog')
    cal = calibration(c,nominal,wb)
    reference_peak, reference_power = peak(c,nominal,wb)
    if not np.isfinite(reference_peak): raise RuntimeError('Reference peak at search boundary')
    rng = np.random.default_rng(c.seed)
    records = []
    for iteration in range(c.trials):
        delta, _ = displacement(c,g,rng)
        actual = nominal+delta
        for name, xyz, w in (
            ('nominal_steered',nominal,wa), ('nominal_fixed',nominal,wb),
            ('deformed_steered',actual,wa), ('deformed_fixed',actual,wb)):
            beam, power = peak(c,xyz,w)
            mono = estimate(c,xyz,w,cal)
            records.append(dict(trial=iteration,scenario=name,
                beam_error_deg=beam-reference_peak,
                monopulse_error_deg=mono-c.target_az_deg,
                relative_peak_db=10*np.log10(power/reference_power) if power>0 else np.nan))
    return c,g,records

def main():
    c,g,records = run()
    out=Path('results')
    out.mkdir(exist_ok=True)
    with (out/'synthetic_demo.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f'Synthetic demo: {len(records)} rows; radius={g["radius"]:.4f} m; output={out/"synthetic_demo.csv"}')

if __name__ == '__main__':
    main()
