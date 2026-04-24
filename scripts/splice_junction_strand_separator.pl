#!/usr/bin/env perl
#splice_junction_strand_separator.pl

# This program takes a STAR "SJ.out.tab" junction file and reformats it to a bed format for visualization on an IGV browser. It keeps all of the information from the STAR output and puts it in the junction name.

# Notes on interpretation of output files:
# Motif key:
# 0 = non-canonical
# 1 = GT/AG
# 2 = CT/AC
# 3 = GC/AG
# 4 = CT/GC
# 5 = AT/AC
# 6 = GT/AT

# Column 3 key:
# Junction_name
# Uniq:x is Number of uniquely mapping reads crossing the junction
# multi:x is Number of multi-mapping reads crossing the junction
# MaxOH:x is Maximum spliced alignment overhang
# Motif:x is the splice donor/acceptor motif for this junction (see "Motif key" above)
# Annot:x indicates whether the splice junction is already annotated or not. 0:unannotated, 1:annotated

# Column 4 is number of reads spanning the junction

# USAGE:
# perl /PATH/splice_junction_strand_separator.pl /PATH/inputfiledirectory/*


foreach my $inputfile(@ARGV) {
    open (INF, "<$inputfile") or die "couldn't open input file";
	open (OUT1, ">$inputfile.positive.bed") or die "couldn't open file";
	open (OUT2, ">$inputfile.negative.bed") or die "couldn't open file";
    print OUT1 "#track name=junctions positive_strand\n";
    print OUT2 "#track name=junctions negative_strand\n";
	while (my $line = <INF>) {
		chomp($line);

        my @split_line = split("\t", $line);
        if ($split_line[3] == 1) {
            print OUT1 $split_line[0], "\t", $split_line[1]-1, "\t", $split_line[2], "\tjnct\t", $split_line[6], "\t+\n";
        }
        elsif ($split_line[3] == 2) {
            print OUT2 $split_line[0], "\t", $split_line[1]-1, "\t", $split_line[2], "\tjnct\t", $split_line[6], "\t-\n";
        }
    }
    close(INF);
    close(OUT1);
    close(OUT2);
}
