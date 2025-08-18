"""
Calibration preparation utilities for the KOALA data‑reduction manager
=====================================================================

This module provides a small orchestration layer around *pykoala* to
prepare the set of calibration objects required to reduce a full night of
observations.

Overview
--------
A :class:`CalibrationSet` bundles the individual correction objects in the
order they should be applied to RSS files and datacubes. The typical
calibration products are:

- **Fibre throughput map** – corrects fibre‑to‑fibre sensitivity variations.
- **Spectrograph throughput / spectral response** – accounts for wavelength‑dependent
  efficiency of the full system.
- **Atmospheric telluric absorption model** – corrects molecular absorption bands.
- **Atmospheric extinction** – removes the smooth, airmass‑dependent extinction curve.
- **(Optional) Wavelength offset correction** – fixes small global shifts
  (e.g. when using twilight flats or standards).

Data sources
------------
To build those products, users typically rely on

- Fibre throughput: twilight‑flat or dome‑flat RSS files.
- Spectral response: spectrophotometric standard stars.

Quick start
-----------
You can build a :class:`CalibrationSet` from a YAML configuration
(see :meth:`CalibrationSet.from_config_yml`). Example configuration::

    workdir: "."
    CalibrationSet:
      ThroughputCorrection:
        rss_set: ["/path/twilight_RSS_1.fits", "/path/twilight_RSS_2.fits"]
        # optional wavelength‑offset correction per RSS
        WaveOffsetCorrect:
          plot: true
          nsigma: 5
      AtmosphericExtinctionCorrection:
        file: default  # or a path to an extinction curve text file
      StandardStarsCal:
        # Either a flat list (single star observed multiple times) ...
        # rss_set: ["/path/std1_rss1.fits", "/path/std1_rss2.fits"]
        # ... or a mapping from star name to list of RSS files
        rss_set:
          Feige110: ["/path/Feige110_rss1.fits", "/path/Feige110_rss2.fits"]
        SubstractBackground: true
        WaveOffsetCorrect:
          plot: true

The call

>>> calset = CalibrationSet.from_config_yml("config.yml")

returns a ready‑to‑use object containing the available corrections. Apply to a
list of RSS objects (in place order) using::

    rss_corr = calset.apply(rss_list)

Design notes
------------
- The class only stores correction **objects** from *pykoala* and applies
  them in a user‑defined order (``correct_order``).
- All constructors perform minimal validation and log non‑fatal issues via
  :func:`koala_drpipe.vprint`.

"""
from __future__ import annotations

import os
from typing import Any, Dict, Iterable, List, Optional, Sequence
import yaml

from astropy import units as u

from pykoala.corrections.flux_calibration import FluxCalibration
from pykoala.corrections.throughput import ThroughputCorrection
from pykoala.corrections.wavelength import (
    WavelengthCorrection,
    TelluricWavelengthCorrection,
)
from pykoala.corrections.sky import (
    TelluricCorrection,
    combine_telluric_corrections,
    SkyFromObject,
    SkySubsCorrection,
)
from pykoala.corrections.atmospheric_corrections import AtmosphericExtCorrection, get_adr
from pykoala.corrections.astrometry import AstrometryCorrection
from pykoala.cubing import CubeInterpolator, build_wcs_from_rss

from koala_drpipe import vprint, instrument_config

__all__ = [
    "CalibrationSet",
    "create_throughput",
    "create_stellar_cal_set",
    "create_telluric",
]

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def get_kwargs(config: Any) -> Dict[str, Any]:
    """Return ``config.get('kwargs', {})`` if *config* is a mapping.

    This is a convenience that lets YAML allow both::

        WaveOffsetCorrect: {plot: true, kwargs: {nsigma: 5}}

    and the simplified::

        WaveOffsetCorrect: {plot: true, nsigma: 5}
    """
    if isinstance(config, dict):
        # allow the user to set parameters at the top level or within "kwargs"
        out = dict(config)
        out.update(config.get("kwargs", {}))
        out.pop("kwargs", None)
        return out
    return {}


class CalibrationSet(object):
    """Bundle of calibration/correction objects for a KOALA night.

    Parameters
    ----------
    aaomega_config : instrument_config.AAOmegaConfig or ``None``
        Spectrograph configuration. Stored for traceability. It is **not**
        used directly by this module but may be useful for downstream code.
    **kwargs : Any
        Additional attributes to attach to the instance. Most importantly,
        you may pass previously built corrections under the keys used in
        :attr:`correct_order` (see below).

    Attributes
    ----------
    correct_order : list[str]
        The attribute names (on this instance) that will be applied in order
        by :meth:`apply`. Default order is::

            ["throughput_corr", "telluric_corr", "atm_ext_corr", "flux_cal_corr"]

    throughput_corr : :class:`ThroughputCorrection` or ``None``
    telluric_corr   : :class:`TelluricCorrection` or ``None``
    wave_corr       : :class:`WavelengthCorrection` or ``None``
    flux_cal_corr   : :class:`FluxCalibration` or ``None``
    atm_ext_corr    : :class:`AtmosphericExtCorrection` or ``None``
    """

    # --- aaomega configuration ------------------------------------------------
    @property
    def aaomega_config(self) -> Optional[instrument_config.AAOmegaConfig]:
        """AAOmega configuration for this calibration set."""
        return getattr(self, "_aaomega_config", None)

    @aaomega_config.setter
    def aaomega_config(self, value: Optional[instrument_config.AAOmegaConfig]):
        self._aaomega_config = value

    # --- correction objects ---------------------------------------------------
    @property
    def throughput_corr(self) -> Optional[ThroughputCorrection]:
        return getattr(self, "_throughput_corr", None)

    @throughput_corr.setter
    def throughput_corr(self, value: Any) -> None:
        if isinstance(value, ThroughputCorrection) or value is None:
            self._throughput_corr = value
        else:
            vprint("ThroughputCorrection not set: wrong type")

    @property
    def telluric_corr(self) -> Optional[TelluricCorrection]:
        return getattr(self, "_telluric_corr", None)

    @telluric_corr.setter
    def telluric_corr(self, value: Any) -> None:
        if isinstance(value, TelluricCorrection) or value is None:
            self._telluric_corr = value
        else:
            vprint("TelluricCorrection not set: wrong type")

    @property
    def wave_corr(self) -> Optional[WavelengthCorrection]:
        return getattr(self, "_wave_corr", None)

    @wave_corr.setter
    def wave_corr(self, value: Any) -> None:
        if isinstance(value, WavelengthCorrection) or value is None:
            self._wave_corr = value
        else:
            vprint("WavelengthCorrection not set: wrong type")

    @property
    def flux_cal_corr(self) -> Optional[FluxCalibration]:
        return getattr(self, "_flux_cal_corr", None)

    @flux_cal_corr.setter
    def flux_cal_corr(self, value: Any) -> None:
        if isinstance(value, FluxCalibration) or value is None:
            self._flux_cal_corr = value
        else:
            vprint("FluxCalibration not set: wrong type")

    @property
    def atm_ext_corr(self) -> Optional[AtmosphericExtCorrection]:
        return getattr(self, "_atm_ext_corr", None)

    @atm_ext_corr.setter
    def atm_ext_corr(self, value: Any) -> None:
        if isinstance(value, AtmosphericExtCorrection) or value is None:
            self._atm_ext_corr = value
        else:
            vprint("AtmosphericExtCorrection not set: wrong type")

    # --- ctor -----------------------------------------------------------------
    def __init__(self, aaomega_config: Optional[instrument_config.AAOmegaConfig], **kwargs: Any) -> None:
        self.aaomega_config = aaomega_config
        for key, val in kwargs.items():
            setattr(self, key, val)
        # default correction order (can be overridden via kwargs)
        if not hasattr(self, "correct_order"):
            self.correct_order = [
                "throughput_corr",
                "telluric_corr",
                "atm_ext_corr",
                "flux_cal_corr",
            ]

    def apply(self, rss: Sequence[Any], correct_order: Optional[List[str]] = None) -> List[Any]:
        """Apply the configured corrections to each RSS in *rss*.

        The method loops over the input sequence and returns a new list with
        corrected objects. Corrections that are ``None`` are silently skipped.

        Parameters
        ----------
        rss : sequence
            A sequence of RSS‑like objects compatible with the ``apply``
            methods of the underlying *pykoala* correction classes.
        correct_order : list[str], optional
            If given, overrides the instance :attr:`correct_order` for this call.

        Returns
        -------
        list
            The list of corrected RSS objects.
        """
        order = list(correct_order) if correct_order is not None else list(self.correct_order)
        out: List[Any] = []
        for item in rss:
            corrected = item
            for corr_name in order:
                corr = getattr(self, corr_name, None)
                if corr is None:
                    vprint(f"Correction {corr_name} not set – skipping")
                    continue
                vprint(f"Applying correction {corr_name}")
                corrected = corr.apply(corrected)
            out.append(corrected)
        return out

    @classmethod
    def from_config_yml(cls, yaml_file):
        """Create a CalibrationSet from an input configuration yaml file."""
        with open(yaml_file, "r") as file:
            config = yaml.safe_load(file)
        return cls.from_config_dict(config)

    @classmethod
    def from_config_dict(cls, config):
        vprint("Initialising Corrections from config file")
        # Initialise throughput
        if "workdir" in config:
            workdir = config["workdir"]
        else:
             workdir="."

        # if "from_rss" in config["AAOmegaConfig"]:
        #     aaomega_config = instrument_config.AAOMegaConfig.from_fits(
        #         config["AAOmegaConfig"]["from_rss"])
        # else:
        #     raise ValueError("User must include a configuration of AAOmega")

        # Initialise corrections
        if "ThroughputCorrection" in config["CalibrationSet"]:
            throughput_corr = create_throughput(
                config["CalibrationSet"]["ThroughputCorrection"],
                workdir=workdir)
        else:
            throughput_corr = None
        
        if "AtmosphericExtinctionCorrection" in config["CalibrationSet"]:
            if config["CalibrationSet"]["AtmosphericExtinctionCorrection"
                                        ]["file"].lower() == "default":
                atm_ext_corr = AtmosphericExtCorrection.from_text_file(
                    AtmosphericExtCorrection.default_extinction)
            else:
                atm_ext_corr = AtmosphericExtCorrection.from_text_file(
                    config["CalibrationSet"]["AtmosphericExtinctionCorrection"
                                        ]["file"])
        else:
            atm_ext_corr = None

        if "StandardStarsCal" in config["CalibrationSet"]:
            print("Preparing StandardStarsCal set")
            telluric_corr, flux_cal_corr = create_stellar_cal_set(
                config["CalibrationSet"]["StandardStarsCal"],
                throughput_corr,
                atm_ext_corr,
                workdir=workdir)

        else:
            telluric_corr = None
            flux_cal_corr = None

        return cls(None,
                   throughput_corr=throughput_corr,
                   atm_ext_corr=atm_ext_corr,
                   telluric_corr=telluric_corr,
                   flux_cal_corr=flux_cal_corr)


def create_throughput(config, workdir="."):
    """Create a thoughput correction from an input configuration.
    
    Parameters
    ----------
    config : dict
        Configuration for creating the :class:`ThroughputCorrection`.

    Returns
    -------
    throuput_corr : :class:`ThroughputCorrection`
    """
    if "file" in config:
        ThroughputCorrection.from_file(config["file"])
    elif "rss_set" in config:
        rss_set = [
            instrument_config.koala_ifu.koala_rss(fl) for fl in config["rss_set"]]
        # Correct wavelength shifts when using Twilight exposures
        if "WaveOffsetCorrect" in config:
            kwargs = get_kwargs(config["WaveOffsetCorrect"])
            for ith, rss in enumerate(rss_set):
                wave_corr, figures = TelluricWavelengthCorrection.from_rss(
                    rss, **kwargs)
                rss = wave_corr.apply(rss)
                if figures is not None:
                    figures[0].savefig(
                        os.path.join(
                            workdir,
                            f"throughput_wavecorr_{ith}_rss_{rss.info['name']}_wave_offset.png"),
                            dpi=200, bbox_inches="tight")
                    figures[1].savefig(
                        os.path.join(
                            workdir,
                            f"throughput_wavecorr_{ith}_rss_{rss.info['name']}_offset_fibre_map.png"),
                            dpi=200, bbox_inches="tight")
        kwargs = get_kwargs(config)
        throughput_corr = ThroughputCorrection.from_rss(rss_set, **kwargs)
    return throughput_corr

def create_stellar_cal_set(config, throughput_corr=None, atm_ext_corr=None, workdir="."):
    """Build the calibration set from standard stars."""
    vprint("Preparing calibrations from standard stars")
    if isinstance(config["rss_set"], list):
        rss_set = [instrument_config.koala_ifu.koala_rss(fl) for fl in config["rss_set"]]
        stars = {rss_set[0].info["name"]: rss_set}
    else:
        stars = {name: [instrument_config.koala_ifu.koala_rss(
            fl) for fl in files] for name, files in config["rss_set"].items()}

    # rss_config = [instrument_config.ObservationConfig.from_fits(
    #     fl) for fl in config["rss_set"]]
    vprint(f"Number of input std. stars: {len(stars)}")
    star_cubes = []
    for std_name, star_rss in stars.items():
        # If requested, first correct for wavelength offset

        for ith, rss in enumerate(star_rss):
            if "WaveOffsetCorrect" in config:
                wave_corr, figures = TelluricWavelengthCorrection.from_rss(
                    rss, **get_kwargs(config["WaveOffsetCorrect"]))
                rss = wave_corr.apply(rss)
                if figures is not None:
                    figures[0].savefig(
                        os.path.join(
                            workdir,
                            f"{std_name}_std_wavecorr_{ith}_rss_wave_offset.png"),
                            dpi=200, bbox_inches="tight")
                    figures[1].savefig(
                        os.path.join(
                            workdir,
                            f"{std_name}_std_wavecorr_{ith}_rss_offset_fibre_map.png"),
                            dpi=200, bbox_inches="tight")
            if throughput_corr is not None:
                vprint("Applying ThroughputCorrection to RSS")
                rss = throughput_corr.apply(rss)
            if atm_ext_corr is not None:
                vprint("Applying AtmosphericExtinctionCorrection to RSS")
                rss = atm_ext_corr.apply(rss)

            if config.get("SubstractBackground", True):
                
                skymodel = SkyFromObject(rss, bckgr_estimator='mad',
                                         source_mask_nsigma=3, remove_cont=False)
                skycorrection = SkySubsCorrection(skymodel)
                rss, fig = skycorrection.apply(rss, plot=True)
                fig.savefig(
                    os.path.join(
                        workdir,
                        f"{std_name}_std_{ith}_rss_sky_substract.png"),
                        dpi=200, bbox_inches="tight")

        # Register the RSS
        astrom_corr = AstrometryCorrection()
        adr_corr_set = []

        offsets, fig = astrom_corr.register_centroids(star_rss,
                                                      object_name=std_name,
                                                qc_plot=True, centroider='gauss')
        fig.savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{std_name}_register.png"),
                        dpi=200, bbox_inches="tight")

        for ith, (rss, offset) in enumerate(zip(star_rss, offsets)):
            astrom_corr.apply(rss, offset=offset)
            adr_pol_ra, adr_pol_dec, fig = get_adr(rss, max_adr=0.5, pol_deg=2,
                                                   plot=True)
            adr_corr_set.append([adr_pol_ra, adr_pol_dec])
            fig.savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{std_name}_adr.png"),
                        dpi=200, bbox_inches="tight")

        wcs = build_wcs_from_rss(star_rss, spatial_pix_size= 0.5 * u.arcsec,
                                spectra_pix_size=rss.wavelength[1] - rss.wavelength[0])
        interpolator = CubeInterpolator(rss_set=star_rss, wcs=wcs, kernel_scale=1.0,
                                        qc_plots=True)
        cube = interpolator.build_cube(cube_info=dict(name=std_name))
        interpolator.cube_plots["stack_cube"].savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{std_name}_cube_qc.png"),
                        dpi=200, bbox_inches="tight")
        interpolator.cube_plots["weights"].savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{std_name}_cube_weights.png"),
                        dpi=200, bbox_inches="tight")
        # Telluric correction
        telluric_corr, fig = TelluricCorrection.from_model(cube, plot=True)
        fig.savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{cube.info['name']}_telluric_correction.png"),
                        dpi=200, bbox_inches="tight")
        cube = telluric_corr.apply(cube)
        star_cubes.append(cube)

    # Flux calibration
    extract_args = dict(wave_range=None, wave_window=None, plot=True)
    response_params = dict(pol_deg=7, spline=False, median_filter_n=10,
                           plot=True)

    flux_cal_results, _, master_flux_corr = FluxCalibration.auto(
        data=star_cubes,
        calib_stars=list(stars.keys()),
        fnames=None,
        extract_args=extract_args,
        response_params=response_params,
        combine=True)

    for name in stars.keys():
        flux_cal_results[name]['extraction']['figure'].savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{name}_extraction.png"),
                        dpi=200, bbox_inches="tight")
        flux_cal_results[name]['response_fig'].savefig(
                    os.path.join(
                        workdir,
                        f"stdcalset_{name}_spectral_response.png"),
                        dpi=200, bbox_inches="tight")
    return telluric_corr, master_flux_corr


def create_telluric(config, workdir="."):
    if "file" in config:
        TelluricCorrection.from_text_file(config["file"])
    elif "rss_set" in config:
        rss_set = [
            instrument_config.koala_ifu.koala_rss(fl) for fl in config["rss_set"]]
        
        # Correct wavelength shifts when using Twilight exposures
        if "WaveOffsetCorrect" in config and config["WaveOffsetCorrect"]:
            for ith, rss in enumerate(rss_set):
                wave_corr, figures = TelluricWavelengthCorrection.from_rss(
                    rss, plot=True)
                rss = wave_corr.apply(rss)
                figures[0].savefig(
                    os.path.join(
                        workdir,
                        f"throughput_wavecorr_{ith}_rss_{rss.info['name']}_wave_offset.png"),
                        dpi=200, bbox_inches="tight")
                figures[1].savefig(
                    os.path.join(
                        workdir,
                        f"throughput_wavecorr_{ith}_rss_{rss.info['name']}_offset_fibre_map.png"),
                        dpi=200, bbox_inches="tight")
        if "kwargs" in config:
            kwargs = config["kwargs"]
        else:
            kwargs = {}
        throughput_corr = ThroughputCorrection.from_rss(rss_set, **kwargs)
    return throughput_corr
