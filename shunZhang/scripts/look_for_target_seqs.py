#!/usr/bin/python

tag = "GTTTAATTGAGTTGTCATATGTTAATAACGGTAT"
target_seq = "GGAGTCACATGGGAGTCACAGGG"

with open("CEP290_1.r2.rc.fastq") as f:
    for i, line in enumerate(f):
        if i % 4 == 1:  # sequence line in FASTQ
            line = line.strip().upper()
            idx = line.find(tag)

            if idx != -1:
                # Extract 30 bp upstream and downstream
                left = line[max(0, idx - 30):idx]
                right = line[idx + len(tag):idx + len(tag) + 30]

                joined = left + right

                #print("joined:", joined)
                if target_seq in joined:
                    print("✅ MATCH to target:", target_seq)
                #else:
                #    print("❌ no match")
