#!/usr/bin/env perl
# developed by David Jebb
use warnings;
use strict ;

## Prior to runnign this script I ran this oneliner to do more work
## grep -P 'HLmyoMyo4_0{5}169' hq_transcripts.cds.gff3 \
## | gt gff3 -setsource isoseq -tidy - \
## | perl -ne '$l=$_;if ($l=~/\t(exon|CDS)\t/){$l=~s/Parent=gene/Parent=mRNA/;print $l}elsif($l=~/\tmRNA\t/){chomp $l;$l=~s/Parent=gene/Parent=mRNA/;\
## ($id)=$l=~/Parent=(\w+)/;$id=";ID=$id\n";print $l.$id}else{print $l}' > tmp

my $count ;
my $fh ;
if (scalar(@ARGV)>0){
	open $fh, "<$ARGV[0]";
}
else{
	$fh = *STDIN ;
}
while (<$fh>){
	my $l = $_ ;
	if ($l =~ /^##( |\w)/){
		print $l
	}
	elsif ($l =~ /^###/){
		$count = 0 ;
		print $l 
	}
	elsif ($l=~/\texon\t/){
		$count++ ;
		chomp $l ;
		my ($id)=$l=~/Parent=(\w+)/ ;
		$id =";ID=$id.$count" ;
		print $l.$id."\n" ;
	}
	elsif ($l=~/\tCDS\t/){
		chomp $l ;
		if ($l=~/ID/){
			print $l."\n"
		}
		else{
			my ($num)=$l=~/Parent=mRNA(\d+)/ ;
			my $id =";ID=CDS$num\n" ;
			print $l.$id ;
		}

	}
	else{
		print $l ;
	}
}	
