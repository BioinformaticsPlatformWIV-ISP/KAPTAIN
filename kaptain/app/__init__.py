from importlib.metadata import version, PackageNotFoundError

_pkg_name = "kaptain"

try:
    __version__ = version(_pkg_name)
except PackageNotFoundError:
    raise RuntimeError(
        "Package metadata not found. "
        "Install the package before importing it."
    )
