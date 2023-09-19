import os
from pathlib import Path
from importlib.resources import files
from shutil import copyfile
import platformdirs
import yaml

CONFIG_VERSION = 1

CONFIG_ORIGIN = files('varnaapi').joinpath('config-settings-v{}.yml'.format(CONFIG_VERSION))
CONFIG_DIR = platformdirs.user_config_path('varnaapi', ensure_exists=True)
CONFIG_USER = Path(CONFIG_DIR, 'config-settings-v{}.yml'.format(CONFIG_VERSION))

CONFIG = {}

def check_settings_exists():
    if not CONFIG_USER.exists():
        copyfile(CONFIG_ORIGIN, CONFIG_USER)


def load_settings():
    global CONFIG
    CONFIG = yaml.load(CONFIG_USER.open(), Loader=yaml.Loader)


def dump_settings():
    yaml.dump(CONFIG, CONFIG_USER.open('w'))


def enable_hack(enable=True, write=False):
    """Enable hack mode.
    When enable, most of the validity check is turned off, except for color.
    """
    global CONFIG
    CONFIG['hackmode'] = enable
    if write:
        dump_settings()

def set_VARNA(path, write=True):
    """Set VARNA location
    """
    global CONFIG
    CONFIG['varnapath'] = path
    if write:
        dump_settings()
