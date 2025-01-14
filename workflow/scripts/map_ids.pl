#!/bin/env perl

# Written by James Fitch and revised for flexibility in input handling.

$quant_file = shift;
$gentrome = shift;

print "Input file: $quant_file\n";

foreach $line (`cat $quant_file`) {
    chomp $line;
    if ($line !~ /^Name/) {
	if ($line =~ /.*\|.*/) {
	    print "Already converted to old ids, aborting\n";
	    exit;
	} else {
	    last;
	}
    }
}

print "Renaming $quant_file -> ${quant_file}.orig\n";
`mv $quant_file ${quant_file}.orig`;

print "Reading $gentrome\n";
foreach $line (`zcat $gentrome`) {
    chomp $line;
    if ($line =~ m/^>(.*?)\|.*/) {
	$old_id = $line;
	$old_id = substr($line,1);
	$new_id = $1;
	
	$new_to_old{$new_id} = $old_id;
    }
}

open OUT, ">$quant_file";

foreach $line (`cat ${quant_file}.orig`) {
    chomp $line;

    if ($line =~ /^Name/) {
	print OUT "$line\n";
    } else {
	@cols = split (/\t/, $line);

	$old_id = $new_to_old{$cols[0]};

	if (! defined($new_to_old{$cols[0]}) ) {
	    $old_id = $cols[0];
	}

	print OUT "$old_id\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\n";
    }
}

close OUT;
