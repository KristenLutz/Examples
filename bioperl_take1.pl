#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use BIO::DB::Genbank;
use Bio::SeqFeatureI;
use Bio::Seq;
#check for cli exit loop if they are not found
if (scalar @ARGV == 2){
	my $file = $ARGV[0];
	my $gene = $ARGV[1];
	my $count = 0;
	my (@genbankID, @translated, $aa_seq, $name, $start, $stop, $strand, $short_seq, $revcom);


	open (my $fh_list, "<", $file) or die "cannot open the file";
	my $out_dna = new Bio::SeqIO(-file => ">dna_kristen_$gene.fa", -format=>'fasta');
	my $out_aa = new Bio::SeqIO(-file => ">aa_kristen_$gene.fa", -format=>'fasta');

	while (<$fh_list>){
		#generates an id list with values in cli 1 file list of numbers
		push(@genbankID, $_);
	}
	#open each id in the small list
	foreach my $i (@genbankID){
		my $db_o = Bio::DB::GenBank -> new;
		my $seq_o = $db_o->get_Seq_by_acc("$i");
		#continue if a genbank db object was successfully found
		if($seq_o){
			#get the species name 
			$name = $seq_o -> species -> binomial('FULL');
			$name =~ s/ /_/g;
			my @seq_object = $seq_o->get_SeqFeatures;

#iterate through sequence features for each CDS
			foreach my $feat_object (@seq_object) {  
#only following actions if it is a CDS object 
				if($feat_object->primary_tag eq 'CDS'){
#pull gene 
					my ($gene_name) = $feat_object->get_tag_values('gene');
#check gene name print add translated to list IF gene in this sequencing object is a gene
					if(($gene_name) eq $gene){
				#pull gene features
						$start = $feat_object->start;
						$stop = $feat_object->end;
						$strand = $feat_object->strand;
						@translated = $feat_object->get_tag_values("translation");
						$aa_seq = "@translated";
						$count++;
						$short_seq = $seq_o->subseq($start,$stop);
						#check for reverse compliment
						if ($strand == -1){
							$short_seq = $seq_o->trunc($start, $stop)->revcom->seq;
						}else{
							$short_seq = $seq_o->subseq($start,$stop);
						}
					}
				}
			}
			#when a gene is found count increments when 0 no gene entries
			#contained gene of interest count will be 0
			if($count==0){
				print "$gene could not be found in genbank entry $i \n";
			}	
		}
	#generate sequence feature and write to file
		$out_aa->write_seq(Bio::Seq->new(-seq=> $aa_seq, -id=> $name));
		$out_dna->write_seq(Bio::Seq->new(-seq=> $short_seq, -id=> $name));
	}
	

#if there are not two cli quit the program
}else{
	print "you need two cli arguments";
}
