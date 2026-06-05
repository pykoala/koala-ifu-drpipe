import logging
from pykoala import __version__, pykoala_logger

KOALA_LOGGER = pykoala_logger.getChild("koala-dr-pipe")

def vprint(msg, *args, **kwargs):
    """
    Convenience function for using with the pykoala generic logger.
    """
    print_method = getattr(KOALA_LOGGER, kwargs.get('level', 'info').lower())
    print_method(msg, *args, **kwargs)
