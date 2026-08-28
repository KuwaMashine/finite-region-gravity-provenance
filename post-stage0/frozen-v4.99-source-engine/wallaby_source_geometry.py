#!/usr/bin/env python3
"""Target-free source-geometry constructors for the WALLABY Gate-A/B holdout.

Nothing in this module reads a rotation curve or residual.  It operates only on
optical surface-photometry profiles and deprojected H I surface-density rings.
The definitions are frozen prospectively in PRETARGET_ADDENDUM_v2.md.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
import numpy as np

EXP_HALF_LIGHT_X = 1.6783469900166605
OPTICAL_MUERR_CUT = 0.22  # mag arcsec^-2, published WALLABY optical truncation


@dataclass(frozen=True)
class OpticalScaleResult:
    r50_arcsec: float
    rd_star_arcsec: float
    n_used: int
    rmax_arcsec: float


@dataclass(frozen=True)
class RingProfile:
    rin: np.ndarray
    rout: np.ndarray
    sigma: np.ndarray
    selected_source_indices: tuple[int, ...]

    @property
    def area(self) -> np.ndarray:
        return math.pi * (self.rout**2 - self.rin**2)

    @property
    def mass_weights(self) -> np.ndarray:
        # Common multiplicative units cancel in mass fractions.
        return self.sigma * self.area


def _ring_edges_from_centres(radius: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = np.asarray(radius, dtype=float)
    if r.ndim != 1 or len(r) < 2 or not np.all(np.isfinite(r)):
        raise ValueError("radius must be a finite 1-D array with >=2 points")
    if np.any(r <= 0) or np.any(np.diff(r) <= 0):
        raise ValueError("radius centres must be positive and strictly increasing")
    mid = 0.5 * (r[:-1] + r[1:])
    rin = np.empty_like(r); rout = np.empty_like(r)
    rin[1:] = mid; rout[:-1] = mid
    rin[0] = max(0.0, r[0] - 0.5*(r[1]-r[0]))
    rout[-1] = r[-1] + 0.5*(r[-1]-r[-2])
    return rin, rout


def optical_exponential_equivalent_scale(
    radius_arcsec,
    mu_r_mag_arcsec2,
    muerr_r_mag_arcsec2,
    ellipticity=None,
) -> OpticalScaleResult:
    """Return the target-free exponential-equivalent stellar scale.

    The published WALLABY optical workflow uses AutoProf r-band geometry and
    truncates the grz profiles when the SB uncertainty reaches 0.22 mag/arcsec^2.
    We retain the r-band profile strictly before the first such point, construct
    a curve of growth, and define Rd,* = R50,r / 1.6783469900166605.  Thus no
    outer-disc fit interval is selected from target behaviour.
    """
    r=np.asarray(radius_arcsec,float)
    mu=np.asarray(mu_r_mag_arcsec2,float)
    emu=np.asarray(muerr_r_mag_arcsec2,float)
    if not (r.shape==mu.shape==emu.shape) or r.ndim!=1:
        raise ValueError("optical arrays must be same-length 1-D")
    if ellipticity is None:
        ell=np.zeros_like(r)
    else:
        ell=np.asarray(ellipticity,float)
        if ell.shape!=r.shape: raise ValueError("ellipticity shape mismatch")
    ok=np.isfinite(r)&np.isfinite(mu)&np.isfinite(emu)&np.isfinite(ell)
    r,mu,emu,ell=r[ok],mu[ok],emu[ok],ell[ok]
    o=np.argsort(r); r,mu,emu,ell=r[o],mu[o],emu[o],ell[o]
    if len(r)<5 or np.any(np.diff(r)<=0) or r[0]<=0:
        raise ValueError("need >=5 strictly increasing positive optical radii")
    bad=np.flatnonzero(emu>=OPTICAL_MUERR_CUT)
    n=int(bad[0]) if len(bad) else len(r)
    r,mu,emu,ell=r[:n],mu[:n],emu[:n],ell[:n]
    if len(r)<5: raise ValueError("fewer than 5 optical points survive the 0.22-mag uncertainty cut")
    q=1.0-ell
    if np.any((q<=0)|(q>1.0)): raise ValueError("ellipticity must imply 0<q<=1")
    rin,rout=_ring_edges_from_centres(r)
    # Relative intensity is sufficient for a half-light radius.
    inten=10.0**(-0.4*(mu-np.nanmin(mu)))
    # Local elliptical-annulus approximation using AutoProf's local q(r).
    flux=inten * math.pi*q*(rout**2-rin**2)
    if not np.all(np.isfinite(flux)) or flux.sum()<=0:
        raise ValueError("invalid optical curve of growth")
    target=0.5*float(flux.sum()); cum=np.cumsum(flux)
    j=int(np.searchsorted(cum,target,side='left'))
    prev=float(cum[j-1]) if j else 0.0
    frac=(target-prev)/float(flux[j])
    # Constant surface brightness and q inside the annulus => flux ~ R^2.
    r50=math.sqrt(rin[j]**2 + frac*(rout[j]**2-rin[j]**2))
    return OpticalScaleResult(r50, r50/EXP_HALF_LIGHT_X, len(r), float(rout[-1]))


def _phase_bit(galaxy_id: str) -> int:
    # Target-independent and stable under row ordering.
    return hashlib.sha256(str(galaxy_id).encode('utf-8')).digest()[0] & 1


def paired_annular_crossfit(radius, sigma_hi, galaxy_id: str) -> tuple[RingProfile,RingProfile,dict]:
    """Construct disjoint predictor and baryonic H I ring profiles.

    Adjacent source annuli are paired.  In each pair exactly one observed annulus
    supplies the predictor-fold surface density and the other supplies the
    baryonic-fold surface density.  The choice alternates by pair and starts from
    a galaxy-ID hash bit, so neither fold is systematically the inner/outer member.
    Each selected density is promoted to the full paired annular area.  If the
    source profile has odd length, the innermost annulus is dropped *before* the
    split so both folds retain the outer radial support needed by R_HI,90.
    """
    r=np.asarray(radius,float); s=np.asarray(sigma_hi,float)
    ok=np.isfinite(r)&np.isfinite(s)&(r>0)&(s>=0)
    r,s=r[ok],s[ok]
    o=np.argsort(r); r,s=r[o],s[o]
    if len(r)<8: raise ValueError("annular cross-fit requires >=8 finite H I rings")
    if np.any(np.diff(r)<=0): raise ValueError("H I radii must be strictly increasing")
    original_indices=np.arange(len(r),dtype=int)
    dropped=None
    if len(r)%2:
        dropped=int(original_indices[0]); r=r[1:]; s=s[1:]; original_indices=original_indices[1:]
    rin0,rout0=_ring_edges_from_centres(r)
    phase=_phase_bit(galaxy_id)
    pr=[]; br=[]; pin=[]; pout=[]
    for j in range(0,len(r),2):
        pair=j//2
        p_local=j + ((phase+pair)&1)
        b_local=j + (1-((phase+pair)&1))
        pr.append(float(s[p_local])); br.append(float(s[b_local]))
        pin.append(float(rin0[j])); pout.append(float(rout0[j+1]))
    pin=np.asarray(pin); pout=np.asarray(pout)
    pred=RingProfile(pin,pout,np.asarray(pr),tuple(int(original_indices[j+((phase+j//2)&1)]) for j in range(0,len(r),2)))
    bary=RingProfile(pin,pout,np.asarray(br),tuple(int(original_indices[j+(1-((phase+j//2)&1))]) for j in range(0,len(r),2)))
    if set(pred.selected_source_indices)&set(bary.selected_source_indices):
        raise AssertionError("cross-fit source annuli overlap")
    meta={"phase_bit":phase,"dropped_source_index":dropped,"n_source_used":len(r),"n_paired_bins":len(pin)}
    return pred,bary,meta


def mass_fraction_radius(profile: RingProfile, fraction: float=0.9) -> float:
    if not (0<fraction<1): raise ValueError("fraction must be in (0,1)")
    w=profile.mass_weights
    total=float(w.sum())
    if not np.isfinite(total) or total<=0: raise ValueError("profile has no positive finite mass")
    target=fraction*total; cum=np.cumsum(w)
    j=int(np.searchsorted(cum,target,side='left'))
    prev=float(cum[j-1]) if j else 0.0
    sigma=float(profile.sigma[j])
    if sigma<=0: raise ValueError("target mass fraction falls in zero-density annulus")
    # w = pi Sigma (Rout^2-Rin^2)
    r2=profile.rin[j]**2 + (target-prev)/(math.pi*sigma)
    return math.sqrt(max(r2,profile.rin[j]**2))


def kappa_hi_from_source_profiles(radius, sigma_hi, galaxy_id: str, rd_star_arcsec: float):
    if not np.isfinite(rd_star_arcsec) or rd_star_arcsec<=0:
        raise ValueError("rd_star_arcsec must be positive")
    pred,bary,meta=paired_annular_crossfit(radius,sigma_hi,galaxy_id)
    r90=mass_fraction_radius(pred,0.90)
    return {
        "R_HI90_arcsec":r90,
        "R_dstar_arcsec":float(rd_star_arcsec),
        "kappa_HI":r90/float(rd_star_arcsec),
        "predictor_profile":pred,
        "baryonic_profile":bary,
        "crossfit_meta":meta,
    }
