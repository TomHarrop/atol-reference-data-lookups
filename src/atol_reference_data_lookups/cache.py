import hashlib
import shelve
from pathlib import Path

from atol_reference_data_lookups import logger


def compute_sha256(file_path):
    logger.debug(f"Computing sha256 checksum for {file_path}.")
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256.update(block)
    hex_digest = sha256.hexdigest()
    logger.debug(f"Checksum: {hex_digest}")
    return hex_digest


def open_cache(cache_dir, name):
    """Open a shelve cache file, ensuring the directory exists."""
    cache_file = Path(cache_dir, name)
    Path.mkdir(cache_file.parent, exist_ok=True, parents=True)
    return cache_file