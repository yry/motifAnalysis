##input: whole and subset of transcripts
##potantial inputs are 5'/3' UTRs, genes
##the subsets can be genes belonging to a particular pathway
##or known to be co expressed from other experiments.

use Bio::SeqIO;
use Algorithm::Loops 'NestedLoops';

sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

##function that returns 1 if the parameters has low sequence complexity
##return 0 otherwise
##the complexity threshold is set to 0.5 i.e. if > 50% of the sequence is low complexity then the sequence is excluded
sub lowSequenceComplexity($)
{
	$seqIn = shift;
	$complexityThr = 0.5;
	
	##create the single and di-nucleotides motifs
	#single nucleotide motifs
	$len = 1;
	@chars  = ('A','C','T','G');
	$get_combo = NestedLoops([ (\@chars) x $len ]);

	#create "motif sequences"
	@combo;
	while ( @combo = $get_combo->() ) {
		$motif = "";
		for($j = 0; $j <= $#combo; $j++)
		{
			$motif = $motif.trim($combo[$j]);
		}
		$motif = trim($motif);
		push @motifs, [($motif)];
	}

	#di nucleotide motifs
	$len = 2;
	$get_combo = NestedLoops([ (\@chars) x $len ]);

	#create "motif sequences"
	@combo;
	while ( @combo = $get_combo->() ) {
		$motif = "";
		for($j = 0; $j <= $#combo; $j++)
		{
			$motif = $motif.trim($combo[$j]);
		}
		$motif = trim($motif);
		push @motifs, [($motif)];
	}
	
	$accepted = 1; # the original accept, this remains universal to the next loops, and governs the filtering process
	$sequence = $seqIn;
	
    for($i = 0; $i <= 3; $i++) # filter 1, based on the first 4 lines of the motif array - single nucleotide
    {
		$count1 = () = $sequence =~ /$motifs[$i][0]/g;
        $prop1  = $count1 / length($sequence); # length of full sequence as single nucleotide
        
		if($prop1 > $complexityThr)
        {
            $accepted = 0;
        }
    }
    
	$sequence = $seqIn;
    for($i = 4; $ i <= 19; $i++) # second filter level based on 4-19 of $motifs - di nucleotide
    {
		$count2 = () = $sequence =~ /$motifs[$i][0]/g;
        $prop2  = $count2 / (length($sequence)/2); # seq length/2 as di nucleotide
		
        if($prop2 > $complexityThr)
        {
            $accepted = 0;
        }
    }
	
	##if accepted == 1 then the sequence is not low complexity
	##if accepted == 0 then the sequence is low complexity
	return $accepted;
}

##########################################################################################
#### MAIN SCRIPT
##########################################################################################

$usage = "perl windowCounting.pl full.set subset win.min win.max";
$full = shift or die $usage;
$sub  = shift or die $usage;
$wmin = shift or die $usage;
$wmax = shift or die $usage;

$seqs = Bio::SeqIO->new(-format => 'fasta', -file => "$full");
while (my $seq = $seqs->next_seq ) 
{
	my $sequence = $seq->seq();
	my $identifier = $seq->display_id();
	push @fullSeq, [(trim($identifier), trim($sequence))];
}
print "full transcript set read\n";

$seqs = Bio::SeqIO->new(-format => 'fasta', -file => "$sub");
while (my $seq = $seqs->next_seq ) 
{
	my $sequence = $seq->seq();
	my $identifier = $seq->display_id();
	push @subSeq, [(trim($identifier), trim($sequence))];
}
print "transcript subset read\n";

##for each entry in the transcript set
##create the motifs with the defined window length 

for($w = $wmin; $w <= $wmax; $w++)
{
	print "working with a window of length $w\n";
	
	##reinitialize the set of motifs
	$#motifsT = -1;
	
	for($k=0; $k<=$#subSeq; $k++)
	{
		$lenT = length($subSeq[$k][1]);
#		print "processing transcript: $subSeq[$k][0] with length $lenT\n";
		for($l = 0; $l < $lenT - $w; $l++)
		{
			$motifT = substr($subSeq[$k][1],$l, $w); ##motif from the transcript
			
			$acc = lowSequenceComplexity($motifT);
			if($acc == 1) ## the accept flag in the subroutine is 1
			{
				push @motifsT,[($motifT)];
			}
		}
	}
	print "motifs after the processing of the subset of transcripts: $#motifsT\n";
	
	##create the set of unique motifs
	@ms = sort{uc($a->[0]) cmp uc($b->[0])}@motifsT;
	$m = $ms[0][0];
	$uniqueMotifs = -1; ##reinitialize for every set of windows
	for($k = 0; $k <= $#ms; $k++)
	{
		if($m ne $ms[$k][0])
		{
			push @uniqueMotifs,[($m)];
			$m = $ms[$k][0];
		}
	}
	push @uniqueMotifs,[($m)];
	
	$allPossible = 4 ** $w;
	print "number of unique motifs: $#uniqueMotifs out of $allPossible possible ones\n";
	
	##count the presence/absence on all genes and the subset genes
	##the significance of these will be assessed using the Fisher exact test
	
	$outputFile = "motifs".$w.".csv";
	open out, ">".$outputFile
		or die "cannot open $outputFile";
	print out "motif, #fullSet, #motifsFullSet, #subSet, #motifsSubSet\n";
	for($k = 0; $k <= $#uniqueMotifs; $k++)
	{
		$motifCountSubset = 0;
		$motifCountFullset= 0;
		
		for($l1 = 0; $l1 <= $#subSeq; $l1++)
		{
			if(index($subSeq[$l1][1],$uniqueMotifs[$k][0]) != -1){$motifCountSubset++}
		}
		for($l2 = 0; $l2 <= $#fullSeq; $l2++)
		{
			if(index($fullSeq[$l2][1],$uniqueMotifs[$k][0]) != -1){$motifCountFullset++}
		}
		print out "$uniqueMotifs[$k][0], $#fullSeq, $motifCountFullset, $#subSeq, $motifCountSubset\n";
	}
	close(out);
	
	##count the motif frequency per gene
	##the significance of these will be assess using the Z score against the number of occurrences on all genes.
}

exit;