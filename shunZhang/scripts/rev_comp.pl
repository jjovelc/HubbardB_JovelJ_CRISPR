#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# Command line options
my $input_file = '';
my $output_file = '';
my $help = 0;

GetOptions(
    'input|i=s'  => \$input_file,
    'output|o=s' => \$output_file,
    'help|h'     => \$help
);

# Help message
if ($help || !$input_file) {
    print_usage();
    exit(0);
}

# Set output to STDOUT if no output file specified
my $out_fh;
if ($output_file) {
    open($out_fh, '>', $output_file) or die "Cannot open output file $output_file: $!\n";
} else {
    $out_fh = \*STDOUT;
}

# Open input file (handle compressed files)
my $in_fh;
if ($input_file =~ /\.gz$/) {
    open($in_fh, "gunzip -c $input_file |") or die "Cannot open compressed file $input_file: $!\n";
} else {
    open($in_fh, '<', $input_file) or die "Cannot open input file $input_file: $!\n";
}

# Process FASTQ file
my $line_count = 0;
my ($header, $sequence, $plus, $quality);

while (my $line = <$in_fh>) {
    chomp $line;
    
    my $position = $line_count % 4;
    
    if ($position == 0) {
        # Header line (starts with @)
        if ($line !~ /^@/) {
            die "Error: Expected header line starting with '@' at line " . ($line_count + 1) . "\n";
        }
        $header = $line;
    }
    elsif ($position == 1) {
        # Sequence line
        $sequence = reverse_complement($line);
    }
    elsif ($position == 2) {
        # Plus line (starts with +)
        if ($line !~ /^\+/) {
            die "Error: Expected '+' line at line " . ($line_count + 1) . "\n";
        }
        $plus = $line;
    }
    elsif ($position == 3) {
        # Quality line
        $quality = reverse($line);  # Reverse quality scores to match reversed sequence
        
        # Print the complete reverse complemented read
        print $out_fh "$header\n";
        print $out_fh "$sequence\n";
        print $out_fh "$plus\n";
        print $out_fh "$quality\n";
    }
    
    $line_count++;
}

# Close file handles
close($in_fh);
if ($output_file) {
    close($out_fh);
}

print STDERR "Processed " . ($line_count / 4) . " reads\n";

# Subroutine to reverse complement DNA sequence
sub reverse_complement {
    my $seq = shift;
    
    # Convert to uppercase for consistency
    $seq = uc($seq);
    
    # Define complement mapping
    my %complement = (
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'N' => 'N',
        'R' => 'Y',  # A or G -> T or C
        'Y' => 'R',  # C or T -> G or A
        'S' => 'S',  # G or C -> C or G
        'W' => 'W',  # A or T -> T or A
        'K' => 'M',  # G or T -> C or A
        'M' => 'K',  # A or C -> T or G
        'B' => 'V',  # C or G or T -> G or C or A
        'D' => 'H',  # A or G or T -> T or C or A
        'H' => 'D',  # A or C or T -> T or G or A
        'V' => 'B'   # A or C or G -> T or G or C
    );
    
    # Reverse and complement the sequence
    my $rev_comp = '';
    for my $base (reverse split //, $seq) {
        if (exists $complement{$base}) {
            $rev_comp .= $complement{$base};
        } else {
            # Keep unknown characters as is
            $rev_comp .= $base;
        }
    }
    
    return $rev_comp;
}

# Print usage information
sub print_usage {
    print << "EOF";
Usage: $0 -i input.fastq [-o output.fastq]

Options:
    -i, --input     Input FASTQ file (required)
    -o, --output    Output FASTQ file (optional, defaults to STDOUT)
    -h, --help      Show this help message

Description:
    This script reverse complements DNA sequences in a FASTQ file.
    Both the sequence and quality scores are reversed to maintain
    proper base-quality correspondence.

Examples:
    # Basic usage
    perl $0 -i input.fastq -o output_revcomp.fastq
    
    # Output to STDOUT
    perl $0 -i input.fastq > output_revcomp.fastq
    
    # Works with compressed files
    perl $0 -i input.fastq.gz -o output_revcomp.fastq

EOF
}
