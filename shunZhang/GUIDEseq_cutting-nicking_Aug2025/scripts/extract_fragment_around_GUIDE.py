#!/usr/bin/env python3
"""
CRISPR Guide Sequence Finder
Finds CRISPR guide sequences in FASTA files and extracts them with flanking regions.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(sequence).reverse_complement())

def find_guide_matches(fasta_file, guide_sequence, flanking_bp=100, search_reverse=True):
    """
    Find guide sequence matches in FASTA file and extract with flanking regions.
    
    Args:
        fasta_file (str): Path to FASTA file
        guide_sequence (str): CRISPR guide sequence to search for
        flanking_bp (int): Number of base pairs to extract on each side
        search_reverse (bool): Whether to also search reverse complement
    
    Returns:
        list: List of dictionaries containing match information
    """
    matches = []
    guide_sequence = guide_sequence.upper()
    guide_rc = reverse_complement(guide_sequence) if search_reverse else None
    
    try:
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence = str(record.seq).upper()
                seq_id = record.id
                seq_description = record.description
                
                # Search forward strand
                start_pos = 0
                while True:
                    pos = sequence.find(guide_sequence, start_pos)
                    if pos == -1:
                        break
                    
                    # Extract flanking sequence
                    extract_start = max(0, pos - flanking_bp)
                    extract_end = min(len(sequence), pos + len(guide_sequence) + flanking_bp)
                    extracted_seq = sequence[extract_start:extract_end]
                    
                    # Calculate relative position of guide within extracted sequence
                    guide_start_in_extract = pos - extract_start
                    guide_end_in_extract = guide_start_in_extract + len(guide_sequence)
                    
                    matches.append({
                        'sequence_id': seq_id,
                        'sequence_description': seq_description,
                        'strand': '+',
                        'position_in_genome': pos + 1,  # 1-based position
                        'guide_sequence': guide_sequence,
                        'extracted_sequence': extracted_seq,
                        'extract_start': extract_start + 1,  # 1-based position
                        'extract_end': extract_end,  # 1-based position
                        'guide_start_in_extract': guide_start_in_extract + 1,  # 1-based position
                        'guide_end_in_extract': guide_end_in_extract,  # 1-based position
                        'extracted_length': len(extracted_seq)
                    })
                    
                    start_pos = pos + 1
                
                # Search reverse complement if requested
                if search_reverse and guide_rc:
                    start_pos = 0
                    while True:
                        pos = sequence.find(guide_rc, start_pos)
                        if pos == -1:
                            break
                        
                        # Extract flanking sequence
                        extract_start = max(0, pos - flanking_bp)
                        extract_end = min(len(sequence), pos + len(guide_rc) + flanking_bp)
                        extracted_seq = sequence[extract_start:extract_end]
                        
                        # Calculate relative position of guide within extracted sequence
                        guide_start_in_extract = pos - extract_start
                        guide_end_in_extract = guide_start_in_extract + len(guide_rc)
                        
                        matches.append({
                            'sequence_id': seq_id,
                            'sequence_description': seq_description,
                            'strand': '-',
                            'position_in_genome': pos + 1,  # 1-based position
                            'guide_sequence': guide_rc,
                            'extracted_sequence': extracted_seq,
                            'extract_start': extract_start + 1,  # 1-based position
                            'extract_end': extract_end,  # 1-based position
                            'guide_start_in_extract': guide_start_in_extract + 1,  # 1-based position
                            'guide_end_in_extract': guide_end_in_extract,  # 1-based position
                            'extracted_length': len(extracted_seq)
                        })
                        
                        start_pos = pos + 1
                        
    except FileNotFoundError:
        print(f"Error: FASTA file '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    return matches

def print_results(matches, output_file=None):
    """Print or save the results."""
    output = []
    
    if not matches:
        output.append("No matches found.")
    else:
        output.append(f"Found {len(matches)} match(es):\n")
        
        for i, match in enumerate(matches, 1):
            output.append(f"=== Match {i} ===")
            output.append(f"Sequence ID: {match['sequence_id']}")
            output.append(f"Description: {match['sequence_description']}")
            output.append(f"Strand: {match['strand']}")
            output.append(f"Position in genome: {match['position_in_genome']}")
            output.append(f"Guide sequence: {match['guide_sequence']}")
            output.append(f"Extracted region: {match['extract_start']}-{match['extract_end']} ({match['extracted_length']} bp)")
            output.append(f"Guide position in extract: {match['guide_start_in_extract']}-{match['guide_end_in_extract']}")
            output.append("")
            output.append("Extracted sequence:")
            
            # Format sequence with guide highlighted
            seq = match['extracted_sequence']
            guide_start = match['guide_start_in_extract'] - 1  # Convert to 0-based
            guide_end = match['guide_end_in_extract'] - 1
            
            # Print sequence in chunks of 80 characters with line numbers
            chunk_size = 80
            for j in range(0, len(seq), chunk_size):
                chunk = seq[j:j+chunk_size]
                line_start = j + match['extract_start']
                output.append(f"{line_start:>8}: {chunk}")
                
                # Add guide highlighting if it's in this chunk
                if guide_start < j + chunk_size and guide_end >= j:
                    highlight_start = max(0, guide_start - j)
                    highlight_end = min(len(chunk), guide_end - j + 1)
                    highlight = " " * 10 + " " * highlight_start + "^" * (highlight_end - highlight_start)
                    if highlight.strip():
                        output.append(highlight)
            
            output.append("\n" + "="*50 + "\n")
    
    result_text = "\n".join(output)
    
    if output_file:
        try:
            with open(output_file, 'w') as f:
                f.write(result_text)
            print(f"Results saved to {output_file}")
        except Exception as e:
            print(f"Error writing to output file: {e}")
            print(result_text)
    else:
        print(result_text)

def save_fasta_results(matches, output_prefix):
    """Save extracted sequences as FASTA files."""
    if not matches:
        print("No matches to save.")
        return
    
    fasta_file = f"{output_prefix}_extracted.fasta"
    
    try:
        with open(fasta_file, 'w') as f:
            for i, match in enumerate(matches, 1):
                header = (f">{match['sequence_id']}_match{i} "
                         f"strand:{match['strand']} "
                         f"pos:{match['position_in_genome']} "
                         f"region:{match['extract_start']}-{match['extract_end']}")
                f.write(f"{header}\n")
                
                # Write sequence in 80-character lines
                seq = match['extracted_sequence']
                for j in range(0, len(seq), 80):
                    f.write(f"{seq[j:j+80]}\n")
        
        print(f"Extracted sequences saved to {fasta_file}")
        
    except Exception as e:
        print(f"Error writing FASTA file: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Find CRISPR guide sequences in FASTA files and extract with flanking regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python crispr_finder.py input.fasta GTCCCTAGTGGCCCCACTGTGGG
  python crispr_finder.py input.fasta GTCCCTAGTGGCCCCACTGTGGG --flanking 200
  python crispr_finder.py input.fasta GTCCCTAGTGGCCCCACTGTGGG --output results.txt --fasta-output extracted
  python crispr_finder.py input.fasta GTCCCTAGTGGCCCCACTGTGGG --no-reverse
        """
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('guide_sequence', help='CRISPR guide sequence to search for')
    parser.add_argument('--flanking', '-f', type=int, default=100, 
                       help='Number of flanking base pairs to extract on each side (default: 100)')
    parser.add_argument('--output', '-o', help='Output file for results (default: print to stdout)')
    parser.add_argument('--fasta-output', help='Prefix for output FASTA file with extracted sequences')
    parser.add_argument('--no-reverse', action='store_true', 
                       help='Do not search for reverse complement')
    
    args = parser.parse_args()
    
    # Validate guide sequence
    valid_bases = set('ATCG')
    if not all(base.upper() in valid_bases for base in args.guide_sequence):
        print("Error: Guide sequence contains invalid characters. Only A, T, C, G are allowed.")
        sys.exit(1)
    
    print(f"Searching for guide sequence: {args.guide_sequence}")
    print(f"Flanking region: {args.flanking} bp on each side")
    print(f"Search reverse complement: {not args.no_reverse}")
    print(f"Input file: {args.fasta_file}\n")
    
    # Find matches
    matches = find_guide_matches(
        args.fasta_file, 
        args.guide_sequence, 
        args.flanking,
        search_reverse=not args.no_reverse
    )
    
    # Print results
    print_results(matches, args.output)
    
    # Save FASTA if requested
    if args.fasta_output:
        save_fasta_results(matches, args.fasta_output)

if __name__ == "__main__":
    main()
