import re

# Colored output
BLUE="\033[94m"
GREEN="\033[92m"
YELLOW="\033[93m"
RED="\033[91m"
PURPLE="\033[95m"
END="\033[0m"

NUCLEOTIDES = ["A", "C", "G", "T"]
QUALITIES = [*range(35, 42)]  # Insertion qualities
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHPX=])")