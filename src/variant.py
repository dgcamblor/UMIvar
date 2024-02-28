# variant.py | Variant class with basic information about a genetic variant

from src.constants import CHR_ORDER


class Variant:
    def __init__(self, chr, pos1, ref, alt, af=None):
        """This class represents basic information about a variant.

        Args:
            chr (str): chromosome name (e.g. "chr1")
            pos1 (int): 1-based position of the variant
            ref (str): ref nucleotide/s
            alt (str): alt nucleotide/s
            af (float, optional): allelic frequency of the variant. Defaults to None.
        """
        self.chr = chr
        self.pos1 = pos1  # 1-based position (usual notation in VCF)
        self.pos0 = pos1 - 1  # 0-based position (python notation)

        self.ref = ref
        self.alt = alt

        self.af = af

        self.type = self.get_type()

        if self.type == "SNV":
            self.length = 1
        elif self.type == "MNV":
            self.length = len(self.ref)
        elif self.type == "INS":
            self.length = len(self.alt) - len(self.ref)
            self.indel_seq = self.alt[1:]
        elif self.type == "DEL":
            self.length = len(self.ref) - len(self.alt)
            self.indel_seq = self.ref[1:]

    def __str__(self):
        if self.af:
            return f"{self.chr}:{self.pos1}{self.ref}>{self.alt} (AF={self.af:.2f})"
        else:
            return f"{self.chr}:{self.pos1}{self.ref}>{self.alt}"
        
    def __eq__(self, other):
        return self.chr == other.chr and \
               self.pos1 == other.pos1 and \
               self.ref == other.ref and \
               self.alt == other.alt
    
    def __lt__(self, other):
        if CHR_ORDER.index(self.chr) < CHR_ORDER.index(other.chr):
            return True
        elif CHR_ORDER.index(self.chr) == CHR_ORDER.index(other.chr):
            return self.pos1 < other.pos1
        else:
            return False
    
    def __hash__(self):
        return hash((self.chr, self.pos1, self.ref, self.alt))
    
    def get_type(self):
        """Returns the type of the variant."""
        if len(self.ref) == len(self.alt) == 1:
            return "SNV"
        elif len(self.ref) == len(self.alt) != 1:
            return "MNV"
        elif len(self.ref) > len(self.alt) and self.ref[0] == self.alt[0]:
            return "DEL"
        elif len(self.ref) < len(self.alt) and self.ref[0] == self.alt[0]:
            return "INS"
        else:
            raise ValueError(f"Invalid variant type for variant: {self}. UMIvar only supports SNVs, MNVs, insertions and deletions.")