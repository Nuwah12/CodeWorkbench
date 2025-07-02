import pysam
import argparse
from collections import defaultdict

def fetch_read_ids(bamfile, chrom, start, window=5000):
    """Fetch unique read IDs mapping to a region ±window bp"""
    read_ids = set()
    for read in bamfile.fetch(chrom, start - window, start + window):
        if not read.is_unmapped:
            read_ids.add(read.query_name)
    return read_ids

def extract_alignments(bamfile, read_ids):
    """Get all alignments for given read IDs"""
    read_alignments = defaultdict(list)
    for read in bamfile.fetch(until_eof=True):
        if read.query_name in read_ids:
            read_alignments[read.query_name].append({
                'chrom': read.reference_name,
                'pos': read.reference_start + 1,  # 1-based
                'cigar': read.cigarstring,
                'flag': read.flag,
                'mapq': read.mapping_quality
            })
    return read_alignments

def main():
    parser = argparse.ArgumentParser(description="Find interchromosomal split reads")
    parser.add_argument("bam", help="Input BAM file (indexed)")
    parser.add_argument("region1", help="First region (e.g., chr3:84873000)")
    parser.add_argument("region2", help="Second region (e.g., chr18:25249000)")
    parser.add_argument("-w", "--window", type=int, default=5000, help="Window size around breakpoints [default: 5000]")
    parser.add_argument("-o", "--output", default="interchr_reads.tsv", help="Output TSV filename")
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    
    bad_name = True
    if bad_name:
        chrom1, y, pos1 = args.region1.split(":")
        chrom2, x, pos2 = args.region2.split(":")
        chrom1 = chrom1+":"+y
        chrom2 = chrom2+":"+x
    else:
        chrom1, pos1 = args.region1.split(":")
        chrom2, pos2 = args.region2.split(":")
    pos1 = int(pos1)
    pos2 = int(pos2)

    print(f"Fetching reads near {args.region1} and {args.region2}...")

    ids1 = fetch_read_ids(bam, chrom1, pos1, args.window)
    ids2 = fetch_read_ids(bam, chrom2, pos2, args.window)

    shared_ids = ids1 & ids2

    print(f"Found {len(shared_ids)} supporting reads.")

    aln_data = extract_alignments(bam, shared_ids)

    with open(args.output, "w") as out:
        out.write("read_id\tchrom\tpos\tcigar\tflag\tmapq\n")
        for rid in sorted(aln_data):
            for aln in aln_data[rid]:
                out.write(f"{rid}\t{aln['chrom']}\t{aln['pos']}\t{aln['cigar']}\t{aln['flag']}\t{aln['mapq']}\n")

    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()

