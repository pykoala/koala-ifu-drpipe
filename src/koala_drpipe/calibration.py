"""
This module contains the tools for preparing the set of calibration objects required for
reducing a full night of observations.

The data that conform a `CalibrationSet` are:
- Fibre throughpupt map: used to correct the inhomogeneities of fibre efficiency.
- Spectrograph throughput: used to account for wavelength-dependent variations of the
spectrograph sensitivity.
- Atmospheric telluric absorption model: a model that corrects the absoprtion of
light by molecules in the atmosphere.


To build such calibration products, users may use the following data:
- Fibre throughput:
    - Twilight flat observations
    - Dome flat observations
- Spectrograph throughput:
    - Observations of standard stars

"""
import os
import yaml

from astropy import units as u

from pykoala.corrections.flux_calibration import FluxCalibration
from pykoala.corrections.throughput import Throughput, ThroughputCorrection
from pykoala.corrections.wavelength import WavelengthCorrection, TelluricWavelengthCorrection
from pykoala.corrections.sky import TelluricCorrection, combine_telluric_corrections
from pykoala.corrections.sky import SkyFromObject, SkySubsCorrection
from pykoala.corrections.atmospheric_corrections import AtmosphericExtCorrection
from pykoala.corrections.astrometry import AstrometryCorrection
from pykoala.corrections.atmospheric_corrections import get_adr
from pykoala.cubing import CubeInterpolator, build_wcs_from_rss

from koala_drpipe import vprint
from koala_drpipe import instrument_config


class CalibrationSet(object):
    """Calibration data.
    
    This class represents the set of calibrations required for reducing
    scientific data. A particular NightCalibration is build for a certain configuration
    of the spectrograph.


    """

    @property
    def aaomega_config(self) -> instrument_config.AAOMegaConfig:
        """Night calibration AAOmega configuration."""
        return self._aaomega_config
    
    @aaomega_config.setter
    def aaomega_config(self, value):
        self._aaomega_config = value

    @property
    def throughput_set(self) -> list:
        """List of :class:`ThroughputCorrection`."""
        return getattr(self, "_throughput_set", None)

    @throughput_set.setter
    def throughput_set(self, value):
        if isinstance(value, ThroughputCorrection):
            self._throughput_set = [value]
        elif isinstance(value, list):
            self._throughput_set = value
        elif value is None:
            return
        else:
            print(f"Unrecognized type {value.__class__}")

    @property
    def telluric_corr_set(self) -> list:
        """List of :class:`TelluricCorrection`."""
        return self._telluric_corr_set
        
    @telluric_corr_set.setter
    def telluric_corr_set(self, value):
        if isinstance(value, TelluricCorrection):
            self._telluric_corr_set = [value]
        elif isinstance(value, list):
            self._telluric_corr_set = value
        elif value is None:
            return
        else:
            print(f"Unrecognized type {value.__class__}")

    @property
    def wave_corr_set(self) -> list:
        """List of :class:`WavelengthCorrection`."""
        return self._wave_corr_set
        
    @wave_corr_set.setter
    def wave_corr_set(self, value):
        if isinstance(value, WavelengthCorrection):
            self._wave_corr_set = [value]
        elif isinstance(value, list):
            self._wave_corr_set = value
        elif value is None:
            return
        else:
            print(f"Unrecognized type {value.__class__}")

    @property
    def flux_cal_corr_set(self) -> list:
        """List of :class:`WavelengthCorrection`."""
        return self._flux_cal_corr_set
        
    @flux_cal_corr_set.setter
    def flux_cal_corr_set(self, value):
        if isinstance(value, FluxCalibration):
            self._flux_cal_corr_set = [value]
        elif isinstance(value, list):
            self._flux_cal_corr_set = value
        elif value is None:
            return
        else:
            print(f"Unrecognized type {value.__class__}")

    @property
    def atm_ext_corr(self) -> list:
        """List of :class:`WavelengthCorrection`."""
        return self._atm_ext_corr
        
    @atm_ext_corr.setter
    def atm_ext_corr(self, value):
        if isinstance(value, AtmosphericExtCorrection):
            self._atm_ext_corr = value
        elif value is None:
            return
        else:
            print(f"Unrecognized type {value.__class__}")


    def __init__(self, aaomega_config, **kwargs):
        self.aaomega_config = aaomega_config

        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

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

        if "from_rss" in config["AAOmegaConfig"]:
            aaomega_config = instrument_config.AAOMegaConfig.from_fits(
                config["AAOmegaConfig"]["from_rss"])
        else:
            raise ValueError("User must include a configuration of AAOmega")

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
            telluric_corr, flux_cal_corr = create_stellar_cal_set(
                config["CalibrationSet"]["StandardStarsCal"],
                throughput_corr,
                atm_ext_corr,
                workdir=workdir)

        else:
            telluric_corr = None

        return cls(aaomega_config,
                   throughput_set=throughput_corr,
                   telluric_corr_set=telluric_corr,
                   flux_cal_corr_set=flux_cal_corr)


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

def create_stellar_cal_set(config, throughput_corr=None, atm_ext_corr=None, workdir="."):
    """Build the calibration set from standard stars."""
    vprint("Preparing calibrations from standard stars")
    rss_set = [instrument_config.koala_ifu.koala_rss(
        fl) for fl in config["rss_set"]]
    rss_config = [instrument_config.ObservationConfig.from_fits(
        fl) for fl in config["rss_set"]]
    vprint(f"Number of input std. stars: {len(rss_set)}")

    if "WaveOffsetCorrect" in config and config["WaveOffsetCorrect"]:
        for ith, rss in enumerate(rss_set):
            wave_corr, figures = TelluricWavelengthCorrection.from_rss(
                rss, plot=True, median_smooth=10, pol_fit_deg=3)
            rss = wave_corr.apply(rss)
            figures[0].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_wavecorr_{ith}_rss_{rss.info['name']}_wave_offset.png"),
                    dpi=200, bbox_inches="tight")
            figures[1].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_wavecorr_{ith}_rss_{rss.info['name']}_offset_fibre_map.png"),
                    dpi=200, bbox_inches="tight")

    if throughput_corr is not None:
        vprint("Applying ThroughputCorrection to RSS")
        for i in range(len(rss_set)):
            rss_set[i] = throughput_corr.apply(rss_set[i])

    if atm_ext_corr is not None:
        vprint("Applying AtmosphericExtinctionCorrection to RSS")
        for i in range(len(rss_set)):
            rss_set[i] = atm_ext_corr.apply(rss_set[i])

    if config.get("SubstractBackground", True):
        for ith in range(len(rss_set)):
            skymodel = SkyFromObject(rss_set[ith], bckgr_estimator='mad',
                                         source_mask_nsigma=3, remove_cont=False)
            skycorrection = SkySubsCorrection(skymodel)
            rss_set[ith], fig = skycorrection.apply(rss_set[ith], plot=True)
            fig.savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{ith}_rss_{rss_set[ith].info['name']}_sky_substract.png"),
                    dpi=200, bbox_inches="tight")

    # Register the RSS
    astrom_corr = AstrometryCorrection()
    star_name = rss_set[0].info['name'].split()[0]
    adr_corr_set = []

    offsets, fig = astrom_corr.register_centroids(rss_set, object_name=star_name,
                                            qc_plot=True, centroider='gauss')
    fig.savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{star_name}_register.png"),
                    dpi=200, bbox_inches="tight")
    for ith, (rss, offset) in enumerate(zip(rss_set, offsets)):
        astrom_corr.apply(rss, offset=offset)
        adr_pol_ra, adr_pol_dec, fig = get_adr(rss, max_adr=0.5, pol_deg=2,
                                               plot=True)
        adr_corr_set.append([adr_pol_ra, adr_pol_dec])
        fig.savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{star_name}_adr.png"),
                    dpi=200, bbox_inches="tight")

    wcs = build_wcs_from_rss(rss_set, spatial_pix_size= 0.5 * u.arcsec,
                             spectra_pix_size=rss.wavelength[1] - rss.wavelength[0])
    interpolator = CubeInterpolator(rss_set=rss_set, wcs=wcs, kernel_scale=1.0,
                                    qc_plots=True)
    cube = interpolator.build_cube(cube_info=dict(name=star_name))
    interpolator.cube_plots["stack_cube"].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{star_name}_cube_qc.png"),
                    dpi=200, bbox_inches="tight")
    interpolator.cube_plots["weights"].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{star_name}_cube_weights.png"),
                    dpi=200, bbox_inches="tight")
    # Telluric correction
    telluric_corr, fig = TelluricCorrection.from_model(cube, plot=True)
    fig.savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{cube.info['name']}_telluric_correction.png"),
                    dpi=200, bbox_inches="tight")
    cube = telluric_corr.apply(cube)
    # Flux calibration
    extract_args = dict(wave_range=None, wave_window=5, plot=True)
    response_params = dict(pol_deg=7, spline=False, median_filter_n=10,
                           plot=True)

    flux_cal_results, _, master_flux_corr = FluxCalibration.auto(
        data=[cube],
        calib_stars=[cube.info['name']],
        fnames=None,
        extract_args=extract_args,
        response_params=response_params,
        combine=True)

    flux_cal_results[cube.info['name']]['extraction']['figure'].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{cube.info['name']}_extraction.png"),
                    dpi=200, bbox_inches="tight")
    flux_cal_results[cube.info['name']]['response_fig'].savefig(
                os.path.join(
                    workdir,
                    f"stdcalset_{cube.info['name']}_spectral_response.png"),
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
