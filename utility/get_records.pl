#!/usr/bin/perl
$i=0;

open HA, ">usa_HA.fasta";
open NA, ">usa_NA.fasta";

while(($title = (<>))) {
	$seq=<>;
	chop($title);
	chop($seq);

	last unless ($title=~/^>/ && length($seq)>0);
	@a=split (/ /,$title);
	$data=$a[0];
	$accnum=$a[1];
	$strain=$a[2];
	$prot=$a[3];
	$country=$a[4];
	$proteins{$prot}++;

	$strains{$strain}{$prot}=$seq;
	$strains{$strain}{"data"}=$data;
	$i++;
}
@k = keys %strains;
print "got $i records\n";
foreach $kk (@k){
	@pp = keys %{$strains{$kk}};

	$data = $strains{$kk}{"data"};
	$seq_na = $strains{$kk}{"NA"};
	$seq_ha = $strains{$kk}{"HA"};
	print "Missing NA! $kk\n" if(length($seq_na)<10);
	print "Missing HA! $kk\n" if(length($seq_ha)<10);

	print NA "$data|$kk|NA\n$seq_na\n\n";
	print HA "$data|$kk|HA\n$seq_ha\n\n";
}
