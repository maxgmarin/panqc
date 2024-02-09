try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    # for python v3.7 or earlier
    from importlib_metadata import PackageNotFoundError, version  # type: ignore

try:
    __version__ = version("pgqc")
except PackageNotFoundError:
    __version__ = "uninstalled"
