##organize the output and attach annotations
##the script appends the annotations to the matrix of significant motifs
##and produces the list of genes with a particular motif

use Bio::SeqIO;

sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

##creates the reverse complement of the given input
sub reverseComplement($_)
{
	$string = shift;
	$revString = reverse($string);
	$revCmp = "";
	for($i = 0; $i <= length($revString); $i++)
	{
		$char = substr($revString,$i,1);
		if($char eq "A"){$revCmp = $revCmp."T"}
		if($char eq "C"){$revCmp = $revCmp."G"}
		if($char eq "G"){$revCmp = $revCmp."C"}
		if($char eq "T" || $char eq "U"){$revCmp = $revCmp."A"}
	}
	return $revCmp;
}

$usage = "perl attachAnnotation.pl anno wmin wmax subset";
$anno  = shift or die $usage;
$wmin  = shift or die $usage;
$wmax  = shift or die $usage;
$subset= shift or die $usage;

##read the annotations
##expected as fasta input
$seqs = Bio::SeqIO->new(-format => 'fasta', -file => "$anno");
while (my $seq = $seqs->next_seq ) 
{
	my $sequence   = trim($seq->seq());
	my $identifier = trim($seq->display_id());
	push @annotations, [($identifier, $sequence)];
	push @annotationsR, [($identifier, reverseComplement($sequence))];
	##the annotation on the positive strand indicates that a motif stems from the transcript
	##the annotation on the negative strand indicates that a motif indicates the targetting location
}
print "full annotation set read with $#annotations entries; both positive and negative strand annotations were created\n";

$seqs = Bio::SeqIO->new(-format => 'fasta', -file => "$subset");
while (my $seq = $seqs->next_seq ) 
{
	my $sequence   = trim($seq->seq());
	my $identifier = trim($seq->display_id());
	push @ss, [($identifier,$sequence)];
}
print "full subset read: $#ss\n";

for($w = $wmin; $w <= $wmax; $w++)
{
	print "processing files for windows of length $w\n";
	
	open inp, "motifs".$w."_significant.csv"
		or die "cannot open input for windows of length $w";
	open out, ">motifs".$w."_significant_withAnno.csv"
		or die "cannot open output for windows with length $w";
	open outM, ">motifs".$w."_withGenes.csv"
		or die "cannot open output (with genes) for windows with length $w";
	print outM "motif; geneSetP; geneSetM\n";
	
	$header = 1;
	while(<inp>)
	{
		chomp;
		@d    = split(/,/);
		if($header == 0)
		{
			$d[1] = trim($d[1]);
			$motif= substr($d[1],1,length($d[1])-2);
			
			$currentAnno = ""; ##the annotation on the positive strand
			$currentAnnoR= ""; ##the annotation on the negative strand
			for($i = 0; $i <= $#annotations; $i++)
			{
				##check whether the annotation transcript (e.g. miRNA, TF) is longer than the motif
				if(length($annotations[$i][1]) > length($motif))
				{
					if(index($annotations[$i][1],$motif) != -1 && length($annotations[$i][0]) > 5) {$currentAnno = $currentAnno.",".$annotations[$i][0]}
					if(index($annotationsR[$i][1],$motif) != -1 && length($annotations[$i][0]) > 5) {$currentAnnoR = $currentAnnoR.",".$annotationsR[$i][0]}
				}
				else
				{
					if(index($motif,$annotations[$i][1]) != -1 && length($annotations[$i][0]) > 5) {$currentAnno = $currentAnno.",".$annotations[$i][0]}
					if(index($motif,$annotationsR[$i][1]) != -1 && length($annotations[$i][0]) > 5) {$currentAnnoR = $currentAnnoR.",".$annotationsR[$i][0]}
				}
			}
			for($j = 0; $j <= $#d; $j++)
			{
				$add[$j] = $d[$j];
				print out "$d[$j],";
			}
			$currentAnno = substr($currentAnno,1);
			$add[$#d + 1] = $currentAnno;
			print out "$currentAnno, $currentAnnoR\n";
			
#			print "$_ $currentAnno, $currentAnnoR\n";
#			sleep(1);
			
			push @inpAnno, [@add];
		}
		else
		{
			print out "$_ annotationP, annotationN\n";
			$header = 0;
		}
	
	}##end while
		
	##for each motif, determine the list of genes associated with it.
	##one list for the positive strand and one list for the negative strand
	for($ii = 0; $ii <= $#inpAnno; $ii++)
	{
		$motif= substr($inpAnno[$ii][1],1,length($inpAnno[$ii][1])-2);
		#print "$ii $motif \n";sleep(1);
		
		$geneSetP = ""; ##list of genes with the motif on positive strand
		$geneSetN = ""; ##list of genes with the motif on negative strand
		for($k = 0; $k <= $#ss; $k++)
		{
			if(index($ss[$k][1],$motif) != -1){$geneSetP = $geneSetP.",".$ss[$k][0]}
			if(index($ss[$k][1],reverseComplement($motif)) != -1){$geneSetN = $geneSetN.",".$ss[$k][0]}
		}
		print outM "$motif;$geneSetP; $geneSetM\n";
	}
	#print "\n\n";
	
	close(inp);
	close(out);
	close(outM);
}

exit;