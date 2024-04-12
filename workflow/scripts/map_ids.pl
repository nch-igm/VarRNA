#!/bin/env perl

#./map_to_old_ids.pl /igm/projects/CH_results/JLE-032/JLE_032_M21-6286_S-2021-032149-R/salmon/quant.sf

$quant_file = shift;

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

print "Reading /igm/apps/salmon/gencode_v28_transcripts_salmon-1.9.0/gentrome.fa.gz\n";
foreach $line (`zcat /igm/apps/salmon/gencode_v28_transcripts_salmon-1.9.0/gentrome.fa.gz`) {
    chomp $line;
    if ($line =~ m/^>(.*?)\|.*/) {
	$old_id = $line;
	$old_id = substr($line,1);
	$new_id = $1;
	
	$new_to_old{$new_id} = $old_id;
	#print "$old_id -> $new_id\n";
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
	    #print "Can't find old id for: $cols[0]\n";
	    $old_id = $cols[0];
	}

	print OUT "$old_id\t$cols[1]\t$cols[2]\t$cols[3]\t$cols[4]\n";
    }
}

close OUT;
