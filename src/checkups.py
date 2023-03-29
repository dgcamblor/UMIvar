# checkups.py | Functions for checking dependencies and directories

import sys, os


def check_dependencies():
    """Check if all dependencies are installed and in PATH."""
    if os.system("samtools --version > /dev/null") != 0:
        print("samtools is not installed or not in PATH")
        sys.exit(1)


def check_dirs():
    """Create the required directories if they don't exist."""
    os.makedirs("tmp", exist_ok=True)