#!/usr/bin/env python3
"""
CRISPR_indel_detect.py - Fixed Version with Error Handling
Detects insertions and deletions at CRISPR cut sites from merged FASTQ reads
"""

import re
import gzip
import os
import sys
from collections import defaultdict, Counter
import argparse

def parse_reference_sequences(ref_file):
    """Parse reference sequences without BioPython dependency"""
    references = {}
    
    print(f"Reading reference file: {ref_file}")
    
    if not os.path.exists(ref_file):
        print(f"ERROR: Reference file not found: {ref_file}")
        sys.exit(1)
    
    try:
        with open(ref_file, 'r') as f:
            content = f.read().strip()
        
        if not content.startswith('>'):
            print(f"ERROR: Reference file should be in FASTA format (start with '>')")
            sys.exit(1)
        
        # Parse FASTA manually
        sequences = {}
        current_seq = ""
        current_id = ""
        
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    sequences[current_id] = current_seq.upper()
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = ""
            elif line:
                current_seq += line
        
        # Don't forget the last sequence
        if current_id and current_seq:
            sequences[current_id] = current_seq.upper()
        
        print(f"Found {len(sequences)} reference sequences")
        
        # For each sequence, find the cut site (assume middle for now)
        for seq_id, sequence in sequences.items():
            # Based on your original data, cut sites should be around position 100
            # Let's try to find the actual cut site by looking for guide sequences
            cut_pos = len(sequence) // 2  # Default to middle
            
            # Your guide sequences (without NGG PAM)
            guide_seqs = {
                'AAVS1': 'GTCCCTAGTGGCCCCACTGT',
                'CEP290': 'CCCTGTGACTCCCATGTGACTCC',  # This is the reverse complement
                'TRAC': 'CCCTGACCCTGCCGTGTACCAGC'   # This is the reverse complement
            }
            
            # Try to find guide sequence in reference to locate cut site
            for guide_name, guide_seq in guide_seqs.items():
                if guide_name.upper() in seq_id.upper():
                    pos = sequence.find(guide_seq)
                    if pos != -1:
                        # Cut is 3bp upstream of PAM, which is after the 20bp guide
                        cut_pos = pos + 17  # Cut between positions 17 and 18 of guide
                        print(f"Found {guide_name} guide in {seq_id} at position {pos}, cut at {cut_pos}")
                        break
                    else:
                        # Try reverse complement
                        rc_guide = reverse_complement(guide_seq)
                        pos = sequence.find(rc_guide)
                        if pos != -1:
                            cut_pos = pos + 17
                            print(f"Found {guide_name} guide (RC) in {seq_id} at position {pos}, cut at {cut_pos}")
                            break
            
            references[seq_id] = {
                'sequence': sequence,
                'cut_position': cut_pos,
                'left_flank': sequence[:cut_pos],
                'right_flank': sequence[cut_pos:]
            }
            
            print(f"  {seq_id}: {len(sequence)}bp, cut at position {cut_pos}")
            
    except Exception as e:
        print(f"ERROR: Cannot read reference file: {e}")
        sys.exit(1)
        
    return references

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in seq[::-1])

def find_best_alignment(read_seq, ref_left, ref_right, max_indel_size=100):
    """
    Find best alignment of read to reference, allowing for indels at cut site
    """
    read_seq = read_seq.upper()
    ref_left = ref_left.upper()
    ref_right = ref_right.upper()
    
    best_alignment = None
    best_score = -1
    
    min_flank_match = 15  # Minimum bases that must match on each flank
    
    # Try different left flank lengths
    for left_len in range(min_flank_match, min(len(ref_left), len(read_seq) - min_flank_match, 50) + 1):
        left_ref = ref_left[-left_len:]  # Take last left_len bases from left flank
        
        # Find where left flank matches in read
        left_match_pos = read_seq.find(left_ref)
        if left_match_pos == -1:
            continue
            
        # Position after left match
        left_end = left_match_pos + left_len
        
        # Try different right flank lengths
        for right_len in range(min_flank_match, min(len(ref_right), len(read_seq) - left_end, 50) + 1):
            right_ref = ref_right[:right_len]  # Take first right_len bases from right flank
            
            # Find where right flank should start
            right_start_pos = read_seq.find(right_ref, left_end)
            if right_start_pos == -1:
                continue
                
            # Calculate indel
            expected_right_start = left_end  # For perfect match
            actual_right_start = right_start_pos
            indel_size = actual_right_start - expected_right_start
            
            if abs(indel_size) > max_indel_size:
                continue
                
            # Calculate alignment score
            score = left_len + right_len
            
            if score > best_score:
                if indel_size > 0:
                    # Insertion
                    indel_seq = read_seq[left_end:right_start_pos]
                    alignment_type = "insertion"
                elif indel_size < 0:
                    # Deletion
                    indel_seq = ""
                    alignment_type = "deletion"
                else:
                    # Perfect match
                    indel_seq = ""
                    alignment_type = "perfect"
                    
                best_alignment = {
                    'type': alignment_type,
                    'size': abs(indel_size),
                    'sequence': indel_seq,
                    'score': score
                }
                best_score = score
    
    return best_alignment

def analyze_reads(fastq_file, references, min_quality=20, max_indel_size=100):
    """Analyze FASTQ reads for indels at cut sites"""
    
    if not os.path.exists(fastq_file):
        print(f"ERROR: FASTQ file not found: {fastq_file}")
        return {}
    
    print(f"File size: {os.path.getsize(fastq_file):,} bytes")
    
    results = {}
    
    for ref_name, ref_data in references.items():
        print(f"\nAnalyzing {ref_name}...")
        
        indel_counts = Counter()
        total_reads = 0
        analyzed_reads = 0
        
        try:
            # Open FASTQ file
            if fastq_file.endswith('.gz'):
                file_handle = gzip.open(fastq_file, 'rt')
            else:
                file_handle = open(fastq_file, 'r')
            
            # Read FASTQ manually (no BioPython dependency)
            line_count = 0
            for line in file_handle:
                line_count += 1
                
                # Every 4th line starting from line 2 is a sequence
                if line_count % 4 == 2:  
                    total_reads += 1
                    read_seq = line.strip().upper()
                    
                    # Quick pre-filter: check if read contains parts of target region
                    left_flank = ref_data['left_flank']
                    right_flank = ref_data['right_flank']
                    
                    # Look for partial matches
                    left_partial = left_flank[-25:] if len(left_flank) >= 25 else left_flank
                    right_partial = right_flank[:25] if len(right_flank) >= 25 else right_flank
                    
                    if left_partial not in read_seq and right_partial not in read_seq:
                        continue
                    
                    analyzed_reads += 1
                    
                    # Try forward orientation
                    alignment = find_best_alignment(read_seq, left_flank, right_flank, max_indel_size)
                    
                    # Try reverse complement
                    if alignment is None or alignment['score'] < 25:
                        rc_read = reverse_complement(read_seq)
                        rc_alignment = find_best_alignment(rc_read, left_flank, right_flank, max_indel_size)
                        if rc_alignment and (alignment is None or rc_alignment['score'] > alignment['score']):
                            alignment = rc_alignment
                    
                    if alignment and alignment['score'] >= 25:  # Minimum score threshold
                        if alignment['type'] == 'perfect':
                            indel_counts['perfect_match'] += 1
                        elif alignment['type'] == 'insertion':
                            key = f"ins_{alignment['size']}bp"
                            if alignment['sequence'] and len(alignment['sequence']) <= 20:
                                key += f"_{alignment['sequence']}"
                            indel_counts[key] += 1
                        elif alignment['type'] == 'deletion':
                            key = f"del_{alignment['size']}bp"
                            indel_counts[key] += 1
                    
                    # Progress update
                    if total_reads % 50000 == 0:
                        print(f"  Processed {total_reads:,} reads, analyzed {analyzed_reads:,}")
                        
            file_handle.close()
            
        except Exception as e:
            print(f"ERROR: Problem reading FASTQ file: {e}")
            return {}
        
        results[ref_name] = {
            'total_reads': total_reads,
            'analyzed_reads': analyzed_reads,
            'indel_counts': indel_counts
        }
        
        print(f"  Total reads: {total_reads:,}")
        print(f"  Analyzed reads: {analyzed_reads:,}")
        print(f"  Reads with events: {sum(indel_counts.values()):,}")
    
    return results

def print_results(results):
    """Print analysis results"""
    
    for ref_name, data in results.items():
        print(f"\n{'='*60}")
        print(f"Results for {ref_name}")
        print(f"{'='*60}")
        
        total_reads = data['total_reads']
        analyzed_reads = data['analyzed_reads']
        indel_counts = data['indel_counts']
        
        print(f"Total reads processed: {total_reads:,}")
        print(f"Reads analyzed (containing target region): {analyzed_reads:,}")
        
        if analyzed_reads > 0:
            print(f"Analysis coverage: {analyzed_reads/total_reads*100:.2f}%")
            
            total_events = sum(indel_counts.values())
            print(f"Total editing events detected: {total_events:,}")
            print(f"Editing efficiency: {total_events/analyzed_reads*100:.2f}%")
            
            if total_events > 0:
                print("\nEvent breakdown:")
                for event_type, count in indel_counts.most_common():
                    percentage = count/analyzed_reads*100
                    print(f"  {event_type}: {count:,} ({percentage:.2f}%)")
            else:
                print("\nNo editing events detected")
        else:
            print("No reads found containing target regions")

def main():
    parser = argparse.ArgumentParser(description='Detect CRISPR indels from FASTQ reads')
    parser.add_argument('fastq', help='Input FASTQ file (can be gzipped)')
    parser.add_argument('reference', help='Reference sequences FASTA file')
    parser.add_argument('--min-quality', type=int, default=20, help='Minimum average read quality')
    parser.add_argument('--max-indel', type=int, default=100, help='Maximum indel size to detect')
    parser.add_argument('--output', help='Output file for results (optional)')
    
    args = parser.parse_args()
    
    # Get sample name from output file for cleaner logging
    sample_name = args.output.replace('_indels.txt', '') if args.output else 'Unknown'
    
    print(f"\n{'='*80}")
    print(f"CRISPR Indel Detection Analysis - Sample: {sample_name}")
    print(f"{'='*80}")
    
    # Parse reference sequences
    print("Loading reference sequences...")
    references = parse_reference_sequences(args.reference)
    
    # Analyze reads
    print(f"\nAnalyzing reads from: {args.fastq}")
    results = analyze_reads(args.fastq, references, args.min_quality, args.max_indel)
    
    if not results:
        print("ERROR: No results generated")
        sys.exit(1)
    
    # Print results
    print_results(results)
    
    # Save results
    if args.output:
        try:
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            
            with open(args.output, 'w') as f:
                f.write(f"CRISPR Indel Detection Results - Sample: {sample_name}\n")
                f.write("=" * 80 + "\n")
                f.write(f"Input file: {args.fastq}\n")
                f.write(f"Reference file: {args.reference}\n")
                f.write(f"Min quality: {args.min_quality}\n")
                f.write(f"Max indel size: {args.max_indel}\n\n")
                
                for ref_name, data in results.items():
                    f.write(f"\n{ref_name}\n")
                    f.write("-" * len(ref_name) + "\n")
                    f.write(f"Total reads: {data['total_reads']:,}\n")
                    f.write(f"Analyzed reads: {data['analyzed_reads']:,}\n")
                    total_events = sum(data['indel_counts'].values())
                    if data['analyzed_reads'] > 0:
                        f.write(f"Editing efficiency: {total_events/data['analyzed_reads']*100:.2f}%\n")
                        f.write(f"Coverage: {data['analyzed_reads']/data['total_reads']*100:.2f}%\n")
                    f.write("\nEvents:\n")
                    if total_events > 0:
                        for event_type, count in data['indel_counts'].most_common():
                            if data['analyzed_reads'] > 0:
                                percentage = count/data['analyzed_reads']*100
                                f.write(f"  {event_type}: {count:,} ({percentage:.2f}%)\n")
                    else:
                        f.write("  No editing events detected\n")
            print(f"\nResults saved to: {args.output}")
            
        except Exception as e:
            print(f"ERROR: Cannot write output file: {e}")
    
    print(f"\nCompleted analysis for {sample_name}")
    print("-" * 80)

if __name__ == "__main__":
    main()
