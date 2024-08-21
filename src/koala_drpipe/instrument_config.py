

def print_aaomega_gratings():
    print("Blue arm gratings: ", "\n".join(["580V" , "1500V" ,"1700B" , "3200B" , "2500V"]))
    print("Red arm gratings: ", "\n".join(["385R","1000R","2000R", "1000I", "1700D","1700I"]))

AAOMEGA_GRATINGS = {
    "blue_arm": ["580V" , "1500V" ,"1700B" , "3200B" , "2500V"],
    "red_arm": ["385R","1000R","2000R", "1000I", "1700D","1700I"]}

class AAOMegaConfig(object):
    def __init__(self, blue_arm_grating, red_arm_grating, *args, **kwargs):
        self.blue_arm_grating

    @classmethod
    def from_header(self):
        """Build a configuration using the information provided in the header."""
        pass


class KOALAConfig(object):
    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def from_header(self):
        """Build a configuration using the information provided in the header."""
        pass


class ObservationConfig(object):
    def __init__(aaomega_config, koala_config, *args, **kwargs):
        pass
