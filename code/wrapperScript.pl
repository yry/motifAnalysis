##wrapper script for the window analysis

#use Statistics::R;

$usage = "window analysis v 1.0;  author: Irina Mohorianu\n
perl wrapperScript.pl fullSet subset win.min win.max pval annotation
fullset = the full set of transcripts e.g. all genes, all 3'/5' UTRs
subset  = the subset on which the window analysis is being done
win.min = start length/ minimum length of the window
win.max = stop length/ maximum length of the window
p.val   = p value threshold for the Fisher significance test
anno    = annotations to be appended after the significance analysis\n
this analysis makes use of Perl (Bio:SeqIO and NestedLoops packages) and R\n\n";

$fullset = shift or die $usage;
$subset  = shift or die $usage;
$wmin    = shift or die $usage;
$wmax    = shift or die $usage;
$pval    = shift or die $usage;
$anno    = shift or die $usage;

##motif based analysis in Perl
##system call on the external Perl script
print "perl windowCounting.pl $fullset $subset $wmin $wmax\n";
#system("perl windowCounting.pl $fullset $subset $wmin $wmax");

##enrichment analysis using a Fisher exact test
##the output consits of entries with p.val < 0.05
#print "executing the R script with the enrichment component\n";
print "Rscript enrichmentAnalysis.R $wmin $wmax $pval\n";
#system("Rscript enrichmentAnalysis.R $wmin $wmax $pval");

##attach annotations and organize outputs
print "perl attachAnnotation.pl $anno $wmin $wmax $subset\n";
system("perl attachAnnotation.pl $anno $wmin $wmax $subset");

exit;