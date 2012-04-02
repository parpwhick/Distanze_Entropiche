#!/usr/bin/perl
$i=0;

open HA, ">usa_HA.fasta";
open NA, ">usa_NA.fasta";
open DATE, ">sequenze_date.txt";

while(($title = (<>))) {
	$seq=<>;
	chop($title);
	chop($seq);
	last unless ($title=~/^>/ && length($seq)>0);
	#eliminazione > dal titolo
	$title = substr($title, 1);
	@a=split (/\|/,$title);
	$data=$a[0];
	$accnum=$a[1];
	$strain=$a[2];
	$prot=$a[3];
	$country=$a[4];

	if($old=$strains{$strain}{$prot}){
		$duplicati++;
		print "Duplicato: $strain $data $prot: FAIL\n" unless $old==$seq;
	}
	$strains{$strain}{$prot}=$seq;
	$strains{$strain}{"data"}=$data;
	$i++;
}
print "read $i records, $duplicati repetitions\n";
$i=0;

#usa file di indicativi
if(0){
	open ID, "<identificativi";
	while($id=<ID>){
		chop($id);
		$ids[$i++]=$id;
	}
	print "read $i identifiers\n";
	$i=0;
	@k=@ids;
}
#stampa tutti quelli trovati
else{
	@k = keys %strains;
}


foreach $kk (@k){
	$data = $strains{$kk}{"data"};
	if( $data=~m#(\d+)/(\d+)/(\d{4})#){
		$mm=sprintf("%02d",$1);
		$dd=sprintf("%02d",$2);
		$yy=$3;
		$data="$yy/$mm/$dd";
	}
	if( !($data=~m#\d{4}/\d{2}/\d{2}#)){
		print "Data incompleta, $kk: $data\n";
		next;
	}
	$seq_na = $strains{$kk}{"NA"};
	$seq_ha = $strains{$kk}{"HA"};
	if(length($seq_na)<10){
		print "Missing NA! $kk, $data\n" ;
		next;
	}
	if(length($seq_ha)<10){
		print "Missing HA! $kk, $data\n" ;
		next;
	}
	$i++;
	print NA ">$data|$kk|NA\n$seq_na\n";
	print HA ">$data|$kk|HA\n$seq_ha\n";
	print DATE "$data\n";
}

print "Found and printed $i records\n";
