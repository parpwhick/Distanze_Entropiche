#!/usr/bin/perl
$i=0;

open HA, ">usa_HA.fasta";
open NA, ">usa_NA.fasta";
open DATE, ">sequenze_date.txt";
open VAC, ">vaccine_data.txt";

$nomefile=$ARGV[0];
$duplicati=0;
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
	$other=$a[5];

	if($old=$strains{$strain}{$prot}){
		$duplicati++;
		print "Duplicato: $strain $data $prot: FAIL\n" unless $old==$seq;
	}
	$strains{$strain}{$prot}=$seq;
	$strains{$strain}{"data"}=$data;
	$strains{$strain}{"other"}=$a[5] if $a[5];
	$i++;
}
print "Read $i records, $duplicati harmless repetitions\n";
if($i < 10) {
	print<<EOF;
Pochi records: affinche il programma funzioni, eseguire prima:
  perl -i -pe '\$i++; chop unless /^>/; s/>/\\n>/ unless (\$i==1);' $nomefile
per formattare correttamente le sequenze!
EOF
	exit(1);
}
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

@k = map  { $_->[0] }
     sort { $a->[1] cmp $b->[1] }
	 #map  { [$_, foo($_)] }
	 map  { [$_, $strains{$_}{"data"}] }
     @k;

$skipped=0;

foreach $kk (@k){
	$data = $strains{$kk}{"data"};
	if( $data=~m#(\d+)/(\d+)/(\d{4})#){
		$mm=sprintf("%02d",$1);
		$dd=sprintf("%02d",$2);
		$yy=$3;
		$data="$yy/$mm/$dd";
	}
	$seq_na = $strains{$kk}{"NA"};
	$seq_ha = $strains{$kk}{"HA"};
	if($other=$strains{$kk}{"other"}){
		$other="|$other";
		#$kk=~m#./(.*)/(\d+)/(\d+)#;
		#       ^------------------ A
		#          ^--------------- Location
		#                ^--------- Seq number
		#                      ^--- Year
		@specials=split("/",$kk);
		#$specialname="$specials[1]$specials[3]";
		$specialname=substr($specials[1],0,6);
	}else{
		$other="";
	}
	if( !($data=~m#\d{4}/\d{2}/\d{2}#) && !$other){
		print "Data incompleta, $kk: $data\n";
		#next;
	}
	if(length($seq_na)<10){
		print "Missing NA! $kk, $data\n" ;
		next;
	}
	if(length($seq_ha)<10){
		print "Missing HA! $kk, $data\n" ;
		next;
	}
	$n_na=$seq_count{$seq_na}++;
	$n_ha=$seq_count{$seq_ha}++;
	if($n_na > 1 && $n_ha > 1 && !$other){
		$skipped++;
			next;
	}
	$i++;
	print NA ">$data|$kk|NA$other\n$seq_na\n";
	print HA ">$data|$kk|HA$other\n$seq_ha\n";
	print DATE "$data\n";
	print "$i is a vaccine of $specialname\n" if $other;
	print VAC "$i $specialname\n" if $other;
}

print "Printed $i records, skipped $skipped (" . sprintf("%.2f%%)\n", $skipped / $#k * 100);
