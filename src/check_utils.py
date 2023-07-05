# check_utils.py | Functions for checking dependencies and directories

import sys, os


def check_dependencies(extra_depends=[]):
    """Check if all dependencies are installed and in PATH.
    
    Args:
        extra_depends (list): Extra dependencies to check.
    """
    if os.system("samtools --version > /dev/null") != 0:
        print("samtools is not installed or not in PATH")
        sys.exit(1)
    if "umi_tools" in extra_depends and os.system("umi_tools --version > /dev/null") != 0:
        print("umi_tools is not installed or not in PATH")
        sys.exit(1)


def check_dirs():
    """Create the required directories if they don't exist."""
    os.makedirs("tmp", exist_ok=True)