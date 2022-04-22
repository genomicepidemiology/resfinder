import io
import os.path

def read_version():
    with io.open(
        os.path.join(os.path.dirname(__file__), "VERSION"), encoding="utf-8"
    ) as f:
        return f.read()

__version__ = read_version()
