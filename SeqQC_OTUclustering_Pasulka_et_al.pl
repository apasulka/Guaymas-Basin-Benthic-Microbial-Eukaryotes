#! /usr/bin/perl -w
open INFILE, "sampleprefix.txt"; #custom sample prefix text file
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;
#Additional scripts required: join_otu_repset.pl, seqlength_cutoff.pl
#QIIME v1.9.1 (09/24/2016)
for $i(@prefix)
{
#Using QIIME to merge
print "join_paired_ends.py -f /galadriel/sarah/deepsea_Pasulka/raw_data/",$i,"_L001_R1_001.fastq -r /galadriel/sarah/deepsea_Pasulka/raw_data/",$i,"_L001_R2_001.fastq -o joined_",$i," -j 25\n"; 
print "mv joined_",$i,"/fastqjoin.join.fastq ",$i,"_merged.fastq\n";

#Split library in QIIME, remove N
print "split_libraries_fastq.py -i ",$i,"_merged.fastq -m map_files/",$i,"_PE_map.txt --barcode_type 'not-barcoded' --sample_ids ",$i," -q 29 -n 0 -o split_",$i,"\n";
print "mv split_",$i,"/seqs.fna ",$i,".merged.Q30.fasta\n";


#Clip primers - allows seqs to have either forward or reverse primers
#search for forward primer, discard seqs without those primers
print "cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --discard-untrimmed -m 10 -o ",$i,".assembled.clipped.regF.fasta ",$i,".merged.Q30.fasta >> ",$i,".filter.log\n";

#perform again, but this time write a new file with the discarded seqs
print "cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --untrimmed-output discard_regF.fasta -m 10 -o tmp.fasta ",$i,".merged.Q30.fasta >> ",$i,".filter.log\n";
print "rm tmp.fasta\n"; #remove tmp (which is a dup)

#Check previously discarded seqs for reverse primer, save them if these is a reverse primer
print "cutadapt -a TYRATCAAGAACGAAAGT -O 3 --discard-untrimmed -m 10 -o ",$i,".assembled.clipped.regR.fasta discard_regF.fasta>> ",$i,".filter.log\n";

print "cat ",$i,".assembled.clipped.regF.fasta ",$i,".assembled.clipped.regR.fasta >> ",$i,".assembled.clipped.fasta\n"; #put them all together!!
print "mv ",$i,"*reg* split_",$i,"/\n"; #move excess files from trimming to split file
print "mv ",$i,".filter.log split_",$i,"/\n";

#Length filter: seqlength_cutoff.pl [input.fasta] [min] [max] [output.fasta] 
print "/beleriand/python_fun/seqlength_cutoff.pl ",$i,".assembled.clipped.fasta 150 500 ",$i,".assembled.clipped.len.fasta\n";

#chimera check with vsearch (uchime)
print "vsearch --uchime_ref ",$i,".assembled.clipped.len.fasta --db /galadriel/sarah/PR2/pr2.qiime.fasta --uchimeout ",$i,".uchimeinfo_ref --chimeras ",$i,".chimeras_ref.fasta --strand plus --nonchimeras ",$i,".assembled.clipped.len.nc.final.fasta \n";

print "mv ",$i,".chimeras_ref.fasta split_",$i,"/\n"; #move excess chimera files to split dir
print "mv ",$i,".uchimeinfo_ref split_",$i,"/\n";

print "cat ",$i,".assembled.clipped.len.nc.final.fasta >> allseqs_deepsea.fasta\n";

} 

#######Written as: clean_getstats.pl
#Cleaning up results:
print "abyss-fac -s 0 *_merged.fastq >>stats_merged.txt\n";
print "abyss-fac -s 0 *.merged.Q30.fasta >> stats_Q30.txt\n";
print "abyss-fac -s 0 *.assembled.clipped.fasta >> stats_primerclipped.txt\n";
print "abyss-fac -s 0 *.assembled.clipped.len.fasta >> stats_length_cutoff.txt\n";
print "abyss-fac -s 0 *.assembled.clipped.len.nc.final.fasta >> stats_chimeras.txt\n";
#import above text files to create seq stats doc via R
# See 'Seq_results.R' script

#start another loop - cleans up directory
open INFILE, "sampleprefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;

for $i(@prefix)
{
#move excess files to split_dir for each sample
print "mv ",$i,"_merged.fastq split_",$i,"/\n";
print "mv ",$i,".merged.Q30.fasta split_",$i,"/\n";
print "mv ",$i,".assembled.clipped.fasta split_",$i,"/\n";
print "mv ",$i,".assembled.clipped.len.fasta split_",$i,"/\n";
}
#########OTU_cluster.pl
#clustering steps
print "pick_de_novo_otus.py -i allseqs_deepsea.fasta -o pick_denovo_uclust.97_deepsea -a -O 4\n";
#Approx 4 million sequences <2 hrs
#Uclust automatically makes BIOM file and assigns taxonomy. This is based on a lame 16S database. So just ignore those!
print "cd pick_denovo_uclust.97_deepsea\n";

#Assign Taxonomy based on 18S PR2!
print "assign_taxonomy.py -i rep_set/allseqs_deepsea_rep_set.fasta -t /galadriel/sarah/PR2/ids.names.2.txt -r /galadriel/sarah/PR2/pr2.qiime.fasta -m uclust --similarity 0.90 --uclust_max_accepts 3 --min_consensus_fraction 0.51 -o uclust_tax_0.90_3_0.51\n";
   # 3 = look up the 3 closest hits in the database
   # 0.90 = in order for it to be a match, it has to be 90 percent similar
   # 0.51 = assign taxonomy which is consistent with 51 percent of the hits (consistent with 2 out of 3)

#OTU table
print "make_otu_table.py -i uclust_picked_otus/allseqs_deepsea_otus.txt -o deepsea_uclust_wPR2.biom -t /galadriel/sarah/deepsea_Pasulka/pick_denovo_uclust.97_deepsea/uclust_tax_0.90_3_0.51/allseqs_deepsea_rep_set_tax_assignments.txt\n";

print "/beleriand/python_fun/join_otu_repset.pl deepsea_uclust_wPR2.biom rep_set/allseqs_deepsea_rep_set.fasta > OTU_table_wRepSeqs.txt\n";