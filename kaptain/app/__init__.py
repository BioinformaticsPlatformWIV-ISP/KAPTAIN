from importlib.metadata import version, PackageNotFoundError

_pkg_name = "kaptain"

try:
    __version__ = version(_pkg_name)
except PackageNotFoundError:
    __version__ = "0.1.0"
