import re
import numpy as np

# Bash scripts
TO_FASTQ = "./src/to_fastq.sh"

# Colored output
BLUE="\033[94m"
GREEN="\033[92m"
YELLOW="\033[93m"
RED="\033[91m"
PURPLE="\033[95m"
END="\033[0m"

# Random mode
ALLELIC_FREQS = np.linspace(0.05, 1, 10000)
POS_DIFF = 170  # Minimum distance between positions
EXTRA_POS = 20  # Extra positions to be determined (for random choice of positions)
VAR_TYPES = ["SNV", "INS", "DEL"]
INDEL_LENGTHS = [*range(1, 4)]  # Insertion and deletion lengths

# Base functioning
CHR_ORDER = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
NUCLEOTIDES = ["A", "C", "G", "T"]
INS_QUALITIES = [*range(35, 42)]  # Insertion qualities
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHPX=])")