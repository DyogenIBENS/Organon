#!/usr/bin/env perl

use warnings;
#use strict;

open(EXO_FILE, $ARGV[0]) || die("Could not open file.") ;

#print "${ARGV[0]}";
#die();

while (<EXO_FILE>) {

	if (/^\s+Target:/) {
		chomp ;
		$target_name = $_ ;
		$target_name =~ s/^\s+Target: //;
	} elsif (/^\s+Query range:/) {
		chomp ;
		$query_range = $_ ;
		#$query_range =~ s/^\s+Query range: // ;
		$query_range =~ s/^\s+// ;
	} elsif (/Target range:/) {
		chomp ;
		$target_range = $_ ;
		#$query_range =~ s/^\s+Query range: // ;
		$target_range =~ s/^\s+// ;
		<EXO_FILE>;<EXO_FILE>;
		while ($_ !~ /^vulgar:/) {
		  <EXO_FILE>;
		  chomp($target_aa .= <EXO_FILE>);
		  chomp($target_na .= <EXO_FILE>);
		  <EXO_FILE>;
		  chomp($_=<EXO_FILE>);
		}
	  }
}

$target_na =~ s/[a-z\s\:\d\-\.]//g;
$target_na =~ s/{.+?}//g;
$target_aa =~ s/[\s\-\+]//g;
$target_aa =~ s/{.+?}//g;

print ">$target_name $query_range $target_range\n$target_na\n";
$target_aa = $target_na = "";
