#! /usr/bin/perl

#length+GC.pl
#Gene W. Tyson

### Input files ###

$n_args = @ARGV;

if ($n_args != 2) {
    print "\nThis program takes a fasta file, extracts length\n"; 
    print "and %GC information\n";
    print "Usage: ./length+GC.pl\n\n";
    print "e.g ./";
    exit;
}

open (LIST, $ARGV[0]) || die "Couldn't open $ARGV[1]\n";

while (my $a = <LIST>) {
    chomp $a;
    my ($name) = $a =~ /^(\S+)/;
#    print "$name\n";
    push (@list, $name);
}
close LIST;

open (CONTIGS, $ARGV[1]) || die "Couldn't open $ARGV[1]\n";

$/= ">";

while (my $b = <CONTIGS>) {
    chomp $b;
    next unless $b;
    my ($name, @sequence) = split (/\n/, $b);
    my ($new_name) = $name =~ /(\S+)/;
    my $seq = join ("", @sequence);
    $sequences{$new_name} = $seq;
}
close CONTIGS;

foreach $value (@list) {
#    print "$value\n";
    if (exists ($sequences{$value})) {
	$seq = $sequences{$value};
	@sequence = (split //, $seq);
	my $size = @sequence;
	print "$value\t$size\n";
    }
    else {
	print "$value\tunknown\n";
    }
}
exit;
