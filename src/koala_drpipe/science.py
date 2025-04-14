
import os
import yaml

from astropy import units as u

from pykoala.data_container import RSS, Cube
from pykoala.instruments.koala_ifu import koala_rss
from pykoala.corrections.atmospheric_corrections import get_adr
from pykoala.cubing import CubeInterpolator, build_wcs_from_rss

from koala_drpipe import vprint
from koala_drpipe import instrument_config


class ScienceSet:
    """"""    
    def __init__(self, files, is_pykoala=False, workdir="."):
        self.files = files
        self.is_pykoala = {key : is_pykoala for key in self.files.keys()}
        self.files_history = {key : [] for key in self.files.keys()}
        self.workdir = workdir

    def calibrate_object(self, name, calibration_set,
                       calibration_corr_order=None, save=False, overwrite=False):
        vprint(f"Processing object: {name}")
        if name not in self.files:
            raise NameError("Input object name not found")
        
        if self.is_pykoala[name]:
            print("READING PYKOALA")
            rss_set = [RSS.from_fits(file) for file in self.files[name]]
        else:
            rss_set = [koala_rss(file) for file in self.files[name]]

        rss_set = calibration_set.apply(rss_set, calibration_corr_order)
        if save:
            self.files_history[name].append(self.files[name].copy())
            self.is_pykoala[name] = True
            for ith, (file, rss) in enumerate(zip(self.files[name], rss_set)):
                if not overwrite:
                    new_file = ".".join(file.split(".")[:-1]) + "_cal.fits"
                    self.files[name][ith] = new_file
                else:
                    new_file = file
                rss.to_fits(new_file, overwrite=True)
        else:
            return rss_set

    def register_object(self, rss):
        pass

    def cube_object(self, rss):
        target_wcs = build_wcs_from_rss(rss)
        interpolator = CubeInterpolator(rss_set=rss, wcs=target_wcs)
        return interpolator.build_cube()

    @classmethod
    def from_config_yml(cls, yaml_file):
        """Create a ScienceSet from an input configuration yaml file."""
        with open(yaml_file, "r") as file:
            config = yaml.safe_load(file)
        if "ScienceSet" not in config:
            raise KeyError("ScienceSet key not found in yml file")
        return cls.from_config_dict(config["ScienceSet"])

    @classmethod
    def from_config_dict(cls, config):
        vprint("Initialising ScienceSet from config file")
        return cls(files=config["files"])