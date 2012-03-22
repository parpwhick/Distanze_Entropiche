  571  perl -i -pe 'chop unless /^>/; s/>/\n>/;' FASTA_1749_chopped.fa 
  586  grep -v '>' FASTA_1746.fa | perl -ne 'print length($_)."\n";'
  587  grep -v '>' FASTA_1746.fa | perl -ne 'print length($_)."\n";' | sort | uniq -c
  623  cat FASTA_1746.fa | grep '>' | grep   -v -n -P '\d\d\d\d/\d\d/\d\d' 
  624  cat FASTA_1746.fa | grep '>' > sequenze_date.txt 
  628  perl -i -pe 's#(\d\d\d\d)// #$1/01/01 #' sequenze_date.txt 
  632  perl -i -pe 's#(\d\d\d\d/\d\d)/ #$1/01 #' sequenze_date.txt 
  639   grep -o  -P '\d\d\d\d/\d\d/\d\d' sequenze_date.txt > sequenze_date2.txt
  640  mv sequenze_date2.txt sequenze_date.txt 
