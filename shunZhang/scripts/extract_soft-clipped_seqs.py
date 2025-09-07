import pysam
import sys

bam = pysam.AlignmentFile(sys.argv[1], "rb")

for read in bam:
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 4:  # 4 = soft clip
                print(f"{read.query_name}\t{read.cigarstring}\t{read.query_sequence}")
                break

