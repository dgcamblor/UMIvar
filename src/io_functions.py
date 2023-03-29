# io_functions.py | Functions for reading and writing in umivar

from subprocess import run
from src.variant import Variant
from src.adding_variants import get_read_id
from multiprocessing import cpu_count


def create_ontarget(bedfile, bamfile):
    """Create a BAM file with only the reads that overlap the target regions.
    Shortens computation time and memory usage.
    
    Args:
        bedfile (file): path to BED file with the target regions.
        bamfile (file): path to BAM file with the reads.
    """
    input_name = bamfile.split("/")[-1].split(".")[0]

    run(f"samtools view -b -L {bedfile} {bamfile} > tmp/{input_name}_ontarget.bam", shell=True)
    run(f"samtools sort tmp/{input_name}_ontarget.bam -o tmp/{input_name}_ontarget.sorted.bam", shell=True)
    run(f"samtools index tmp/{input_name}_ontarget.sorted.bam", shell=True)


def extract_variants(varfile):
    """Extract variants from a CSV file.
    
    Args:
        varfile (file): File with the variants. Each line must have the following format:
            chr,pos1,ref,alt,af
    Returns:
        list: List of Variant objects
    """
    variants = []
    for line in varfile:
        fields = line.rstrip().split(",")
        variants.append(Variant(chr = fields[0], 
                                pos1 = int(fields[1]),
                                ref = fields[2], 
                                alt = fields[3], 
                                af = float(fields[4])))
    return variants


def write_reads(in_bamfile, out_bamfile, record, stats):
    """Write the reads to the output file. If a read has a variant, write the mutated read instead.

    Args:
        in_bamfile (pysam.AlignmentFile): Input BAM file.
        out_bamfile (pysam.AlignmentFile): Output BAM file.
        record (dict): Dictionary with the reads that have a variant.
    Returns:
        dict: Dictionary with the number of reads and variants written.
    """
    cpus = cpu_count()  # Shortens computation time of counting reads
    total = int(run(f"samtools view -@ {cpus} -c {in_bamfile.filename.decode('utf-8')}", shell=True, capture_output=True).stdout.decode("utf-8").rstrip())
    written_reads, variant_reads = 0, 0
    
    for read in in_bamfile:
        if read.reference_name in record and read.reference_start in record[read.reference_name]:
            read_id = get_read_id(read)
            if read_id not in record[read.reference_name][read.reference_start]:
                out_bamfile.write(read)
                written_reads += 1
            else:
                mut_read = record[read.reference_name][read.reference_start][read_id]
                out_bamfile.write(mut_read)
                written_reads += 1
                variant_reads += 1
        else:
            out_bamfile.write(read)
            written_reads += 1

        print(f"Progress: {written_reads}/{total} ({written_reads/total*100:.2f}%)", end="\r")

    stats["written_reads"] = written_reads
    stats["variant_reads"] = variant_reads