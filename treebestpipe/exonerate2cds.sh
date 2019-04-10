#!/bin/bash
# From https://github.com/sujaikumar/assemblage/blob/master/README-meloidogyne.md

set -euo pipefail

exonerate_out=$1

cat $exonerate_out | perl -ne '
  if (/Target range:/) {
    <>;<>;
    while ($_ !~ /^vulgar:/) {
      <>;
      chomp($target_aa .= <>);
      chomp($target_na .= <>);
      <>;
      chomp($_=<>);
    }
    $target_na =~ s/[a-z\s\:\d\-\.]//g;
    $target_na =~ s/{.+?}//g;
    $target_aa =~ s/[\s\-\+]//g;
    $target_aa =~ s/{.+?}//g;
    #print "$_\t$target_aa\t$target_na\n";
    print "$target_na\n";
    $target_aa = $target_na = "";
  }
'
