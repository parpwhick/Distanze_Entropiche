perl -i -pe 'chop unless /^>/; s/>/\n>/;' FASTA_1749_chopped.fa 
grep -v '>' FASTA_1746.fa | perl -ne 'print length($_)."\n";'
grep -v '>' FASTA_1746.fa | perl -ne 'print length($_)."\n";' | sort | uniq -c
cat FASTA_1746.fa | grep '>' | grep   -v -n -P '\d\d\d\d/\d\d/\d\d' 
cat FASTA_1746.fa | grep '>' > sequenze_date.txt 
perl -i -pe 's#(\d\d\d\d)// #$1/01/01 #' sequenze_date.txt 
perl -i -pe 's#(\d\d\d\d/\d\d)/ #$1/01 #' sequenze_date.txt 
 grep -o  -P '\d\d\d\d/\d\d/\d\d' sequenze_date.txt > sequenze_date2.txt
mv sequenze_date2.txt sequenze_date.txt 
