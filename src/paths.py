import os
import sys

# Compute project-relative paths used across modules
utilitiesPath = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "./utilities")
)

# Check for SQANTI3_DIR environment variable first
if "SQANTI3_DIR" in os.environ:
    sqantiqcPath = os.path.abspath(os.environ["SQANTI3_DIR"])
else:
    sqantiqcPath = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../SQANTI3")
    )

# Ensure SQANTI3 python sources are importable
if sqantiqcPath not in sys.path:
    sys.path.insert(0, sqantiqcPath)

# Also ensure utilities path is importable
if utilitiesPath not in sys.path:
    sys.path.insert(0, utilitiesPath)

__all__ = [
    'utilitiesPath',
    'sqantiqcPath',
]


