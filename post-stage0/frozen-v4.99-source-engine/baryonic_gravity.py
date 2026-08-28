#!/usr/bin/env python3
"""Newtonian baryonic acceleration constructor for the WALLABY pretarget lane.

Primary model: axisymmetric razor-thin stellar + gas discs.  The stellar disc is
an exponential-equivalent profile fixed by total Mstar and Rd,*; the gas disc is
piecewise constant on the *baryonic* annuli produced by the disjoint annular
cross-fit.  The gas potential is integrated directly, then differentiated with a
convergence check.  This module never reads observed velocities.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk, iv, kv
from wallaby_source_geometry import RingProfile

G_KPC_KMS2_MSUN = 4.300917270036279e-6
KPC_TO_M = 3.085677581491367e19
KMS2_TO_M2S2 = 1.0e6
ACC_KMS2_PER_KPC_TO_MS2 = KMS2_TO_M2S2 / KPC_TO_M
HELIUM_FACTOR = 1.35  # frozen to Deg et al. 2024 convention


def arcsec_to_kpc(radius_arcsec, distance_mpc: float):
    if not np.isfinite(distance_mpc) or distance_mpc<=0: raise ValueError("distance_mpc must be positive")
    return np.asarray(radius_arcsec,float) * (distance_mpc*1e3) / 206264.80624709636


def exponential_disc_v2_thin(R_kpc, mass_msun: float, rd_kpc: float):
    """Freeman razor-thin exponential disc circular-speed squared [(km/s)^2]."""
    R=np.asarray(R_kpc,float)
    if mass_msun<=0 or rd_kpc<=0: raise ValueError("mass and Rd must be positive")
    y=R/(2.0*rd_kpc)
    out=np.zeros_like(y)
    pos=y>0
    yp=y[pos]
    B=iv(0,yp)*kv(0,yp)-iv(1,yp)*kv(1,yp)
    out[pos]=2.0*G_KPC_KMS2_MSUN*mass_msun/rd_kpc * yp**2 * B
    return out


def _potential_piecewise_thin(R: float, profile: RingProfile, sigma_scale: float=1.0) -> float:
    """Midplane potential [(km/s)^2] of a piecewise-constant axisymmetric disc.

    profile radii are in kpc and Sigma in Msun/pc^2. sigma_scale can include the
    frozen helium factor.  The logarithmic ring singularity is integrable.
    """
    if R<0: raise ValueError("R must be nonnegative")
    total=0.0
    for a,b,s0 in zip(profile.rin,profile.rout,profile.sigma):
        if s0<=0 or b<=a: continue
        Sigma=float(s0)*1e6*float(sigma_scale)  # Msun/kpc^2
        def f(rp):
            den=R+rp
            if den<=0: return 0.0
            m=0.0 if R==0 else 4.0*R*rp/(den*den)
            m=min(max(m,0.0),1.0-1e-14)
            return rp*float(ellipk(m))/den
        points=[R] if a<R<b else None
        val=quad(f,float(a),float(b),points=points,epsabs=1e-8,epsrel=2e-9,limit=160)[0]
        total += Sigma*val
    return -4.0*G_KPC_KMS2_MSUN*total


def gas_disc_g_thin(R_kpc, profile: RingProfile, helium_factor: float=HELIUM_FACTOR, convergence_tol: float=1e-2):
    """Return gas radial acceleration [m/s^2] and a convergence diagnostic."""
    R=np.atleast_1d(np.asarray(R_kpc,float))
    if np.any(R<=0): raise ValueError("evaluation radii must be positive")
    # Characteristic length for derivative step.
    widths=np.asarray(profile.rout)-np.asarray(profile.rin)
    scale=float(np.median(widths[widths>0])) if np.any(widths>0) else float(np.median(R))
    vals=[]; errs=[]
    for x in R:
        # Keep derivative points away from zero and use two nested step sizes.
        h=max(1e-5*max(x,scale), min(2e-3*scale, 5e-4*max(x,scale)))
        h=min(h,0.2*x)
        def deriv(hh):
            return (_potential_piecewise_thin(x+hh,profile,helium_factor)-_potential_piecewise_thin(x-hh,profile,helium_factor))/(2.0*hh)
        d1=deriv(h); d2=deriv(h/2.0)
        g=max(0.0,d2)*ACC_KMS2_PER_KPC_TO_MS2
        rel=abs(d2-d1)/max(abs(d2),1e-12)
        vals.append(g); errs.append(rel)
    arr=np.asarray(vals); er=np.asarray(errs)
    return arr, {"relative_step_change":er,"converged_mask":er<=convergence_tol,"max_relative_step_change":float(np.max(er)),"converged":bool(np.all(er<=convergence_tol)),"tolerance":float(convergence_tol)}


def baryonic_acceleration(
    R_arcsec,
    distance_mpc: float,
    mstar_msun: float,
    rd_star_arcsec: float,
    gas_profile_arcsec: RingProfile,
    helium_factor: float=HELIUM_FACTOR,
):
    """Construct gbar and component accelerations in m/s^2."""
    Rk=arcsec_to_kpc(R_arcsec,distance_mpc)
    rd= float(arcsec_to_kpc([rd_star_arcsec],distance_mpc)[0])
    rin=arcsec_to_kpc(gas_profile_arcsec.rin,distance_mpc)
    rout=arcsec_to_kpc(gas_profile_arcsec.rout,distance_mpc)
    gp=RingProfile(rin,rout,np.asarray(gas_profile_arcsec.sigma,float),gas_profile_arcsec.selected_source_indices)
    v2s=exponential_disc_v2_thin(Rk,float(mstar_msun),rd)
    gs=(v2s/Rk)*ACC_KMS2_PER_KPC_TO_MS2
    gg,diag=gas_disc_g_thin(Rk,gp,helium_factor=helium_factor)
    return {"R_kpc":Rk,"gstar":gs,"ggas":gg,"gbar":gs+gg,"gravity_diagnostic":diag}


def observed_acceleration(R_arcsec, Vrot_kms, distance_mpc: float):
    Rk=arcsec_to_kpc(R_arcsec,distance_mpc)
    V=np.asarray(Vrot_kms,float)
    if Rk.shape!=V.shape or np.any(Rk<=0) or np.any(V<=0): raise ValueError("positive matching R/V arrays required")
    return (V**2/Rk)*ACC_KMS2_PER_KPC_TO_MS2


def normalize_hi_profile_to_mass(profile_arcsec: RingProfile, distance_mpc: float, mhi_msun: float) -> tuple[RingProfile,float]:
    """Renormalize a cross-fit H I profile to the independently measured total H I mass.

    R_HI,90 is invariant under this scaling.  The operation is used only for the
    baryonic gas model, so the galaxy-level morphology predictor does not inherit
    the total-flux normalization.
    """
    if not np.isfinite(mhi_msun) or mhi_msun<=0: raise ValueError("mhi_msun must be positive")
    rin=arcsec_to_kpc(profile_arcsec.rin,distance_mpc); rout=arcsec_to_kpc(profile_arcsec.rout,distance_mpc)
    mass0=float(np.sum(np.asarray(profile_arcsec.sigma,float)*1e6*math.pi*(rout**2-rin**2)))
    if not np.isfinite(mass0) or mass0<=0: raise ValueError("cross-fit H I profile has no finite mass")
    scale=float(mhi_msun/mass0)
    out=RingProfile(np.asarray(profile_arcsec.rin,float),np.asarray(profile_arcsec.rout,float),np.asarray(profile_arcsec.sigma,float)*scale,profile_arcsec.selected_source_indices)
    return out,scale
