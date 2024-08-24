
from pykoala import vprint
from astropy.io import fits

def print_aaomega_gratings():
    print("Blue arm gratings: ", "\n".join(["580V" , "1500V" ,"1700B" , "3200B" , "2500V"]))
    print("Red arm gratings: ", "\n".join(["385R","1000R","2000R", "1000I", "1700D","1700I"]))

AAOMEGA_GRATINGS = {
    "blue_arm": ["580V" , "1500V" ,"1700B" , "3200B" , "2500V"],
    "red_arm": ["385R","1000R","2000R", "1000I", "1700D","1700I"]}


class AAOMegaConfig(object):
    def __init__(self, arm, grating, dichroic, exptime, *args, **kwargs):
        self.arm = arm
        self.grating = grating
        self.dichroic = dichroic
        self.exptime = exptime

    @classmethod
    def from_header(cls, header):
        """Build a configuration using the information provided in the header."""
        arm = header["SPECTID"]
        grating = header["GRATID"]
        dichroic = header["DICHROIC"]
        exptime = header["EXPOSED"]

        return cls(arm, grating, dichroic, exptime)

    @classmethod
    def from_fits(cls, path, extension=0):
        header = fits.getheader(path, extension)
        return cls.from_header(header)

    def show_config(self):
        vprint(f"AAOMega configuration:\n - Arm: {self.arm}"
               + f"\n - Grating: {self.grating}\n - Dichroic: {self.dichroic}"
               + f"\n - Exp. Time: {self.exptime}")


class KOALAConfig(object):
    def __init__(self, mode, rotator_pa=0.0, *args, **kwargs):
        self.mode = mode
        self.rotator_pa = rotator_pa

    @classmethod
    def from_header(cls, header):
        """Build a configuration using the information provided in the header."""
        rotator_pa = header["TEL_PA"]
        mode = header.get("FOV", "N/A")
        return cls(mode, rotator_pa)

    @classmethod
    def from_fits(cls, path, extension=0):
        header = fits.getheader(path, extension)
        return cls.from_header(header)

    def show_config(self):
        vprint(f"KOALA configuration:\n - Mode: {self.mode}"
               + f"\n - Rotator PA: {self.rotator_pa}")


class ObservationConfig(object):
    def __init__(self, aaomega_config, koala_config, mean_ra, mean_dec, utdate, **kwargs):
        self.aaomega_config = aaomega_config
        self.koala_config = koala_config
        self.mean_ra, self.mean_dec = mean_ra, mean_dec
        self.utdate = utdate
        # Convert to python datetime
        self.utstart = kwargs.get("utstart")
        self.utend = kwargs.get("utend")
        self.utmjd = kwargs.get("utmjd")
        self.zenital_distance = kwargs.get("zdstart")
        self.ver_2dfdr = kwargs.get("ver_2dfdr")
        self.is_raw = kwargs.get("is_raw")
    
    @classmethod
    def from_header(cls, header):
        aaomega_config = AAOMegaConfig.from_header(header)
        koala_config = KOALAConfig.from_header(header)

        obs_args = {"mean_ra": header["MEANRA"], "mean_dec": header["MEANDEC"],
                    "utdate": header["UTDATE"], "utmjd": header["UTMJD"], 
                    "utstart": header["UTSTART"], "utend": header["UTEND"],
                    "is_raw": False}
        obs_args["ver_2dfdr"] = header.get('2DFDRVER', None)
        if obs_args["ver_2dfdr"] is None:
            obs_args["is_raw"] = True

        return cls(aaomega_config, koala_config, **obs_args)

    @classmethod
    def from_fits(cls, path, extension=0):
        header = fits.getheader(path, extension)
        return cls.from_header(header)
    
    def show_config(self):
        vprint(f"Observation configuration:\n - Mean Ra/Dec: {self.mean_ra}/{self.mean_dec}"
               + f"\n - UT date: {self.utdate}"
               + f"\n - UT start: {self.utstart}"
               + f"\n - Is raw file: {self.is_raw}"
               + f"\n - 2dfdr version: {self.ver_2dfdr}")
        self.aaomega_config.show_config()
        self.koala_config.show_config()


def group_observations(list_of_obs_config, criteria):
    pass

