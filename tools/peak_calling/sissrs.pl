#!/usr/bin/perl 

#====================================================================
#
# SISSRs: Site Identification from Short Sequence Reads
#
# Author: Raja Jothi
#         National Institutes of Health
#         jothi@mail.nih.gov
#
# SISSRs is a tool for identifying genome-wide transcription factor
# binding sites from ChIP-Seq data generated. In ChIP-Seq, which
# combines chromatin immunoprecipitation (ChIP) with massively
# parallel sequencing, the DNA fragments obtained from ChIP are
# directly sequenced using next-generation genome sequencers.
#
# A binding site is a region on the DNA to which specific proteins
# including, but not limited to, transcription factors bind in vivo.
# A typical binding site be anywhere from ~5-20 nucleotides in length.
#
# SISSRs takes as input a file containing genomic coordinates of the
# mapped short sequence reads in BED file format, and outputs a list
# of binding sites. Each inferred binding site is reported as a
# genomic coordinate representing the center of the binding motif.
#
# If you use this program in your research, please cite:
#
# Raja Jothi, Suresh Cuddapah, Artem Barski, Kairong Cui, Keji Zhao
# Genome-wide identification of in vivo protein-DNA binding sites
# from ChIP-Seq Data, Nucleic Acids Research, 36(16):5221-31 (2008)
#
#====================================================================

my $old_fh = select(STDOUT);
$| = 1;
select($old_fh); 
$numArgs = $#ARGV + 1;
use Cwd;
use Getopt::Std;
my %Options;
$ok = getopts('i:o:s:F:D:b:p:m:w:E:L:e:q:tacrux', \%Options);
if($numArgs < 6 || !ok || !exists($Options{i}) || !exists($Options{o}) || !exists($Options{s}) ) 
{ 
  printUsage();
  exit;
}


### SETTING DEFAULT PARAMETERS, IF NECESSARY

### Genome length
if ( $Options{s} < 1)	
  { print "\n\nERROR: -s option should be set to a positive integer value (>0)\n\n"; exit; }
else
  { $genomeLength = $Options{s}; }

### number of tags necessary on either side of transition point t (see Fig 1 in the paper)
if ( !exists($Options{E}) )
  { $numReadsOnEachSide = 2; }
elsif ( $Options{E} < 1)
  { print "\n\nERROR: -E option should be set to a positive integer value (>0)\n\n"; exit }
else
  { $numReadsOnEachSide = $Options{E}; }


### option to keep only one read if multiple reads aling the same genomic coordinate
if ( !exists($Options{a}) )
  { $keepOneRead = 0; }
else
  { $keepOneRead = 1; }

### option to report binding sites with tags mapped to only strand (see Fig 6B in the paper)
if ( !exists($Options{u}) )
  { $oneside = 0; }
else
  { $oneside = 1; }

### option to NOT print progress report during execution 
if ( !exists($Options{x}) )
  { $progress = 1; }
else
  { $progress = 0; }

### option to report identified binding sites as singular genomic coordinates
if ( !exists($Options{t}) )
  { $reportSingularCoord = 0; }
else
  { $reportSingularCoord = 1; }

### option to report identified binding sites as regions instead of a singular genomic coordinate
if ( !exists($Options{r}) )
  { $reportRegions = 0; }
else
  { $reportRegions= 1; }

### option to report identified binding sites as clusters of binding sites (currently set as default)
if ( !exists($Options{c}) )	### regardless of whether user sets -c or not, this is the default setting
  { $reportClusters = 0; }
else
  { $reportClusters = 1; }

### scanning window size, w (see Fig 1 in the paper)
if ( !exists($Options{w}) )
  { $interval = 20; }
elsif ( $Options{w}%2 != 0 || $Options{w} < 2)
  { print "\n\nERROR: -w option should be set to an *even* number greater than 1\n\n"; exit;}
else
  { $interval = $Options{w}; } 		### size of window to scan the reads; default 20

### upper-bound of the DNA fragment length
if ( !exists($Options{L}) )
  { $upperBound = 500; }
elsif ( $Options{L} < 1)
  { print "\n\nERROR: -L option should be set to a positive integer\n\n"; exit;}
else
  { $upperBound = $Options{L}; } 	

### fraction of genome that is mappable by sequenced reads
if ( !exists($Options{m}) )
  { $mappableFraction = 0.8; }
elsif ( $Options{m} < 0.0 || $Options{m} > 1.0)
  { print "\n\nERROR: -m option should be set to a *real* number between 0.0 and 1.0, inclusive\n\n"; exit; }
else
  { $mappableFraction = $Options{m}; } 	

$bedFileName = $Options{i};		### input (BED) file containing sequenced reads

### HASHING GENOMIC REGIONS (e.g., SATELLITE REPEAT COORDINATES) WHICH NEED TO BE IGNORED 
if ( exists($Options{q}) )
{
  $repeatsFileName = $Options{q};
  open(SAT, "$repeatsFileName") || die "can't open file $repeatsFileName\n";
  %repeatStarts = ();
  %repeatStops = ();
  $ignoredRegionLength = 0;
  while(<SAT>)
  {
    chomp;
    @words = split /\s+/, $_;	### 0-chr 1-start 2-end 
    $chr = $words[0]; 
    $start = $words[1]; 
    $end = $words[2];
    if ( !exists($repeatStarts{$chr}) )
    {
      $repeatStarts{$chr} = $start;
      $repeatStops{$chr} = $end;
    }
    else
    {
      $repeatStarts{$chr} .= " ".$start;
      $repeatStops{$chr} .= " ".$end;
    }
    $ignoredRegionLength += ($end-$start+1);
  }
  close SAT;
  $ignoreFraction = $ignoredRegionLength/$genomeLength;
}

print "\nProcessing...\n\n" if ($progress);

### HASHING BED COORDINATES (FROM DATA FILE)
print "Reading data file..." if ($progress);
$ctr = 0;
$tagsSelected = 0;
$numPositive = 0;
$numNegative = 0;
open(BED, "$bedFileName") || die "can't open file $bedFileName\n";
%bedStartsPos = ();
%bedStartsNeg = ();
%allReads = ();
while(<BED>)
{
  $ctr++;
  chomp;
  @words = split /\s+/, $_;	### 0-chrom 1-bedStart 2-bedEnd 5-strand

  $tagIsRepeat = 0;
  if ( exists($Options{q}) )
  {
    @myRepeatStarts = split /\s/, $repeatStarts{$words[0]};
    @myRepeatStarts = sort { $a <=> $b} @myRepeatStarts;
    @myRepeatStops= split /\s/, $repeatStops{$words[0]};
    @myRepeatStops= sort { $a <=> $b} @myRepeatStops;
    $startIndex = binarySearch(0, $words[1], \@myRepeatStarts);	# left off tag start
    $stopIndex = binarySearch(1, $words[2], \@myRepeatStarts);	# right off tag end 
    for($index = $startIndex; $index <= $stopIndex; $index++)
    {
      $satStart = $myRepeatStarts[$index];
      $satEnd = $myRepeatStops[$index];
      if ($words[1] < $satEnd && $words[2] > $satStart)
      {
        $tagIsRepeat = 1;
        last;
      }
    }
  }    

  if ( $tagIsRepeat==0 && $words[5] eq "+" )
  {
    $tmpStr = $words[0]."_".$words[1].$words[5];
    if ( !exists($bedStartsPos{$words[0]}) )
    {
      $bedStartsPos{$words[0]} = $words[1];
      $tagsSelected++;
      $numPositive++;
      $allReads{$tmpStr} = 1 if ($keepOneRead == 1);
    }
    else
    {
      if ($keepOneRead == 1)
      {
	if ( !exists($allReads{$tmpStr}) )
	{
          $bedStartsPos{$words[0]} .= " ".$words[1];
          $tagsSelected++;
          $numPositive++;
          $allReads{$tmpStr} = 1;
	}
      }
      else
      {
        $bedStartsPos{$words[0]} .= " ".$words[1];
        $tagsSelected++;
        $numPositive++;
      }
    }
  }
  elsif ( $tagIsRepeat==0 && $words[5] eq "-" )
  {
    $tmpStr = $words[0]."_".$words[2].$words[5];
    if ( !exists($bedStartsNeg{$words[0]}) )
    {
      $bedStartsNeg{$words[0]} = $words[2];
      $tagsSelected++;
      $numNegative++;
      $allReads{$tmpStr} = 1 if ($keepOneRead == 1);
    }
    else
    {
      if ($keepOneRead == 1)
      {
	if ( !exists($allReads{$tmpStr}) )
	{
          $bedStartsNeg{$words[0]} .= " ".$words[2];
          $allReads{$tmpStr} = 1;
          $tagsSelected++;
          $numNegative++;
	}
      }
      else
      {
        $bedStartsNeg{$words[0]} .= " ".$words[2];
        $tagsSelected++;
        $numNegative++;
      }
    }
  }
}
close BED;
undef %allReads;
print "done!\n" if ($progress);
$fracSelected = $tagsSelected/$ctr*100;
$posFrac = $numPositive/$tagsSelected*100;
$negFrac = $numNegative/$tagsSelected*100;

### HASHING CONTROL COORDINATES (BACKGROUND)
if ( exists($Options{b}) )
{
  print "Reading control file..." if ($progress);
  $controlTagsSelected = 0;
  open(CON, "$Options{b}") || die "can't open file $Options{b}\n";
  %bedStartsPosControl = ();
  %bedStartsNegControl = ();
  %allReadsControl = ();
  while(<CON>)
  {
    chomp;
    @words = split /\s+/, $_;	### 0-chrom 1-bedStart 2-bedEnd 5-strand
  
    $tagIsRepeat = 0;
    if ( exists($Options{q}) )
    {
      @myRepeatStarts = split /\s/, $repeatStarts{$words[0]};
      @myRepeatStarts = sort { $a <=> $b} @myRepeatStarts;
      @myRepeatStops= split /\s/, $repeatStops{$words[0]};
      @myRepeatStops= sort { $a <=> $b} @myRepeatStops;
      $startIndex = binarySearch(0, $words[1], \@myRepeatStarts);	# left off tag start
      $stopIndex = binarySearch(1, $words[2], \@myRepeatStarts);	# right off tag end 
      for($index = $startIndex; $index <= $stopIndex; $index++)
      {
        $satStart = $myRepeatStarts[$index];
        $satEnd = $myRepeatStops[$index];
        if ($words[1] < $satEnd && $words[2] > $satStart)
        {
          $tagIsRepeat = 1;
          last;
        }
      }
    }    
  
    if ( $tagIsRepeat==0 && $words[5] eq "+" )
    {
      $tmpStr = $words[0]."_".$words[1].$words[5];
      if ( !exists($bedStartsPosControl{$words[0]}) )
      {
        $bedStartsPosControl{$words[0]} = $words[1];
        $controlTagsSelected++;
        $allReadsControl{$tmpStr} = 1 if ($keepOneRead == 1);
      }
      else
      {
        if ($keepOneRead == 1)
        {
          if ( !exists($allReadsControl{$tmpStr}) )
          {
            $bedStartsPosControl{$words[0]} .= " ".$words[1];
            $controlTagsSelected++;
            $allReadsControl{$tmpStr} = 1;
          }
        }
        else
        { 
          $bedStartsPosControl{$words[0]} .= " ".$words[1];
          $controlTagsSelected++;
        }
      }
    }
    elsif ( $tagIsRepeat==0 && $words[5] eq "-" )
    {
      $tmpStr = $words[0]."_".$words[2].$words[5];
      if ( !exists($bedStartsNegControl{$words[0]}) )
      {
        $bedStartsNegControl{$words[0]} = $words[2];
        $controlTagsSelected++;
        $allReadsControl{$tmpStr} = 1 if ($keepOneRead == 1);
      }
      else
      {
        if ($keepOneRead == 1) 
        {
          if ( !exists($allReadsControl{$tmpStr}) )
          {
            $bedStartsNegControl{$words[0]} .= " ".$words[2];
            $controlTagsSelected++;
            $allReadsControl{$tmpStr} = 1;
          } 
        }
        else
        {
          $bedStartsNegControl{$words[0]} .= " ".$words[2];
          $controlTagsSelected++;
        }
      }
    }
  }
  close CON;
  undef %allReadsControl;

  ### p-value, which determines the minimal fold enrichment (real to background) necessary
  ### to call a candidate binding site as the true binding site
  if ( !exists($Options{p}))
    { $fold_pValue = 0.001; }
  elsif ( $Options{p} < 0 || $Options{p} > 1)
    { print "\n\nERROR: -p option should be set to a *real* number between 0 and 1, inclusive\n\n"; exit;}
  else
    { $fold_pValue = $Options{p}; }

  ### e-value, which determines the number of tags necessary to identify candidate binding sites
  if ( !exists($Options{e}))
    { $eValue = 10; }
  elsif ( $Options{e} < 0 )
    { print "\n\nERROR: -e option should be set to a *real* number >=0\n\n"; exit;}
  else
    { $eValue = $Options{e}; }
  print "done!\n" if ($progress);
}
else
{
  ### false discovery rate, if random background model based on Poisson probabilities need to be used
  if ( !exists($Options{D}))
    { $fdr = 0.001; }
  elsif ( $Options{D} < 0 || $Options{D} > 1)
    { print "\n\nERROR: -D option should be set to a *real* number between 0 and 1, inclusive\n\n"; exit;}
  else
    { $fdr = $Options{D}; }
}

### length of the sequenced DNA fragment
if ( !exists($Options{F}) )
{
  $fragmentLength = estimateFragmentLength();
}
elsif ( $Options{F} > 0)
{
  $fragmentLength = $Options{F};
}
else
{
  die "\n\nERROR:\t-F option should be set to an integer value >0\n\t\t\t\tor\n\t\tleft unset in which case the fragment length is estimated from reads\n"
}

$windowLookout = $fragmentLength;

$outputFile = $Options{o};	### file into which SISSRs output is stored

$mappableFraction -= $ignoreFraction if ( exists($Options{q}) ); ### subtracting ignored regions
$genomeLength = int($genomeLength*$mappableFraction);	### effective genome length used for the random model 
$numSlots = $genomeLength/(2*$fragmentLength);
$mean = ($tagsSelected/$genomeLength)*(2*$fragmentLength);

### initial estimate for R (see Fig 1 in the paper) if random background model is used
$evidNecessary = 2*$numReadsOnEachSide;	

### calculating the value for R based on e-value, if a background control data is used
if ( exists($Options{b}) )
{
  $i = $evidNecessary;
  $pValue = 1.0;
  for ($k = 0; $k <= $i-1; $k++)
  {
    $pValue -= Poisson($mean, $k);
  }
  $byChance = $numSlots*$pValue;
  $k = $i - 1;
  while($byChance > $eValue)
  {
    $k++;
    $i++;
    $tmp = Poisson($mean, $k);
    $pValue -= Poisson($mean, $k);
    $byChance = $numSlots*$pValue;
  }
  $evidNecessary = $i;
}

printSummary();

print "\nFinding candidate binding sites..." if ($progress);
findBindingSites();
print "done!\n" if ($progress);

system("sort -k 1.4,1.5 -k 2n $outputFile.tmp > $outputFile.tmp.tmp");
system("mv $outputFile.tmp.tmp $outputFile.tmp");

print "Printing the binding sites onto the output file..." if ($progress); 

if ( !exists($Options{b}) )
{
  ### Estimating the number of tags necessary to define binding sites
  ### based on random background model (Poisson probabilities)

  print "Estimating the number of tags necessary\n  to define binding sites based on random\n  background model (Poisson probabilities)..." if ($progress);

  open(FL, "$outputFile.tmp") || die "can't find file $outputFile.tmp\n";
  %support = ();
  $totalSites = 0;
  while(<FL>)
  {
    chomp;
    @words = split /\t/, $_;	### 0-chr 1-start 2-end 3-evid
    if (exists($support{$words[3]}))
    {
      $support{$words[3]}++;
    }
    else
    {
      $support{$words[3]} = 1;
    }
    $totalSites++;
  }
  close FL;
  
  $trueEvidNecessary = $tagsSelected;
  $i = $evidNecessary;
  while($totalSites > 0)
  {
    $pValue = 1.0;
    for ($k = 0; $k <= $i-1; $k++)
    {
      $pValue -= Poisson($mean, $k);
    }
    $byChance = $numSlots*$pValue;
    $trueFDR = $byChance/$totalSites;
    if ($trueFDR <= $fdr)
    {
      $trueEvidNecessary = $i;
      last;
    }
    else
    {
      $totalSites -= $support{$i} if ( exists($support{$i}) );
    }
    $i++;
  }
  if ($totalSites <= 0)
  {
    print "\nThe chosen FDR=$fdr is so stringent that no binding sites\n\tcould be identified. Try increasing the value for D\n\n";
    exit;
  }
  print "done!\n\n" if ($progress);
}
else
{
  $trueEvidNecessary = $evidNecessary;
}

$outputFile = $Options{o};
open(OUT, ">>$outputFile") || die "can't open file $outputFile\n";
print OUT "======================================================================\n";
print OUT "BINDING SITES\n";
print OUT "======================================================================\n";

if ( !exists($Options{b}) )
{
  printf OUT "Tags necessary to identify binding sites: $trueEvidNecessary (FDR < %.2e)\n", $fdr;
  print OUT "(at least E=$numReadsOnEachSide tags on each side)\n\n";
}
else
{
  $numTrials = 1000000;		### to estimate fold distribution between real and background data 
  $numTrials = int(1000/$fold_pValue) if ( 1000/$fold_pValue > $numTrials );
  @foldArray = ();
  pValueDistribution();	### to store the fold distribution values
  @foldArray = sort {$b <=> $a} @foldArray;	### non-decreasing order
  $foldCutoff = $foldArray[$fold_pValue*($#foldArray+1)-1];
  printf OUT "Tags necessary to identify binding sites: $evidNecessary (Fold >= %.2f)\n", $foldCutoff;
  printf OUT "(at least E=$numReadsOnEachSide tags on each side)\n\n"; 
}

if ($reportRegions == 0 && $reportSingularCoord == 1 && $reportClusters == 0)
{
  if ( !exists($Options{b}) )
  {
    print OUT "Chr\tBsite\tNumTags\n";
    print OUT "---\t-----\t-------\n";
  }
  else
  {
    print OUT "Chr\tBsite\tNumTags\tFold\tp-value\n";
    print OUT "---\t-----\t-------\t----\t-------\n";
  }
  open(FL, "$outputFile.tmp") || die "can't find file $outputFile.tmp\n";
  while(<FL>)
  {
    chomp;
    @words = split /\t/, $_;	### 0-chr 1-start 2-end 3-evid [4-fold]
    $midpoint = int(($words[1]+$words[2])/2);
    if ($words[3] >= $trueEvidNecessary && !exists($Options{b}) )
    {
      print OUT "$words[0]\t$midpoint\t$words[3]\n";
    }
    elsif ( exists($Options{b}) && $words[4] >= $foldCutoff)
    {
      $site_pValue = get_pValue($words[4]);
      printf OUT "$words[0]\t$midpoint\t$words[3]\t%.2f\t%.1e\n", $words[4], $site_pValue;
    }
  }
  close FL;
}
elsif ( $reportRegions == 1 && $reportClusters == 0)
{
  if ( !exists($Options{b}) )
  {
    print OUT "Chr\trStart\trEnd\tNumTags\n";
    print OUT "---\t------\t----\t-------\n";
  }
  else
  {
    print OUT "Chr\trStart\trEnd\tNumTags\tFold\tp-value\n";
    print OUT "---\t------\t----\t-------\t----\t-------\n";
  }
  open(FL, "$outputFile.tmp") || die "can't find file $outputFile.tmp\n";
  while(<FL>)
  {
    chomp;
    @words = split /\t/, $_;      ### 0-chrom 1-start 2-end 3-numEvid [4-fold]
    if ( !exists($Options{b}) && $words[3] >= $trueEvidNecessary )
    {
      print OUT "$_\n";
    }
    elsif ( exists($Options{b}) && $words[4] >= $foldCutoff)
    {
      $site_pValue = get_pValue($words[4]);
      printf OUT "$words[0]\t$words[1]\t$words[2]\t$words[3]\t%.2f\t%.1e\n", $words[4], $site_pValue;
    }
  }
  close FL;
}
else ### $reportClusters == 1 
{
  if ( !exists($Options{b}) )
  {
    print OUT "Chr\tcStart\tcEnd\tNumTags\n";
    print OUT "---\t------\t----\t-------\n";
  }
  else
  {
    print OUT "Chr\tcStart\tcEnd\tNumTags\tFold\tp-value\n";
    print OUT "---\t------\t----\t-------\t----\t-------\n";
  }
  open(FL, "$outputFile.tmp") || die "can't find file $outputFile.tmp\n";
  $line = <FL>;
  chomp($line);
  @words = split /\t/, $line;     ### 0-chr 1-start 2-end 3-numEvid [4-fold]
  while($words[3] < $trueEvidNecessary)
  {
    if ($line = <FL>)
    {
      chomp($line);
      @words = split /\t/, $line;     
    }
    else
    {
      if ( !exists($Options{$b}) )
      {
        print "\nFDR=$fdr is so stringent that no binding sites\n\tcould be identified. Try increasing the value for D\n\n";
      }
      else
      {
        print "\nE-value=$eValue is so stringent that no binding sites\n\tcould be identified. Try increasing the value for e\n\n";
      }
      exit;
    }
  }
  $prevChr = $words[0];
  $prevMidpoint = int( ($words[1]+$words[2])/2 );
  $evid = $words[3];
  $fold = $words[4] if ( exists($Options{b}) );
  $bsiteStart = $words[1];
  $bsiteEnd = $words[2];
  while(<FL>)
  {
    chomp;
    @words = split /\t/, $_;      ### 0-chrom 1-start 2-end 3-numEvid [4-fold]
    if ($words[3] >= $trueEvidNecessary)
    {
      $currChr = $words[0];
      $currMidpoint = int( ($words[1]+$words[2])/2 );
      $currStart = $words[1];
      $currEnd = $words[2];
      $currEvid = $words[3];
      $currFold = $words[4] if ( exists($Options{b}) );
      if ($currChr eq $prevChr && $currMidpoint - $prevMidpoint <= $fragmentLength) ### overlap
      {
        $bsiteEnd = $currEnd if ($bsiteEnd < $currEnd);
        $evid = $currEvid if ($currEvid > $evid);
        $fold = $currFold if ($currFold > $fold && exists($Options{b}) );
        $prevMidpoint = $currMidpoint;
      }
      else
      {
        if ( !exists($Options{b}) )
	{
          print OUT "$prevChr\t$bsiteStart\t$bsiteEnd\t$evid\n";
	}
	elsif ($fold >= $foldCutoff)
	{
          $site_pValue = get_pValue($fold);
          printf OUT "$prevChr\t$bsiteStart\t$bsiteEnd\t$evid\t%.2f\t%.1e\n", $fold, $site_pValue;
	}
        $prevChr = $currChr;
        $prevMidpoint = $currMidpoint;
        $evid = $currEvid;
        $fold = $currFold if ( exists($Options{b}) );
        $bsiteStart = $currStart;
        $bsiteEnd = $currEnd;
      }
    }
  }
  if ( !exists($Options{b}) )
  {
    print OUT "$prevChr\t$bsiteStart\t$bsiteEnd\t$evid\n";
  }
  elsif ($fold >= $foldCutoff)
  {
    $site_pValue = get_pValue($fold);
    printf OUT "$prevChr\t$bsiteStart\t$bsiteEnd\t$evid\t%.2f\t%.1e\n", $fold, $site_pValue;
  }
  close FL;
}
print OUT "======================================================================\n";
close OUT;
print "done!\n" if ($progress);

system("rm $outputFile.tmp");

#==================================== END MAIN ====================================


sub Poisson
{
  my($lambda, $k) = @_;
  my($exp) = exp 1;
  return $exp**(-1.0*$lambda) * $lambda**$k / Factorial($k);
}

sub Factorial
{ 
  my($num) = @_;
  if ($num == 0)
  { 
    return 1;
  }
  else
  { 
    return $num * Factorial($num-1);
  }
}


sub estimateFragmentLength	### estimates the average DNA fragment length
{
  my($avgLength) = 0; 
  my($cnt) = 0;
  my($chrom, @posTags, @negTags, $i, $retIndex, $fragmentLength, $sub);
  foreach $chrom (sort keys%bedStartsPos)
  {
    @posTags = split /\s/, $bedStartsPos{$chrom};
    @posTags = sort {$a <=> $b} @posTags;
    @negTags = split /\s/, $bedStartsNeg{$chrom};
    @negTags = sort {$a <=> $b} @negTags;
    $sub = 0;
    for ($i = 0; $i <= $#posTags; $i++)
    {
      $retIndex = binarySearch(1, $posTags[$i], \@negTags);
      if ($retIndex >= 0 && $retIndex <= $#negTags)
      {
        $fragLength = $negTags[$retIndex] - $posTags[$i] + 1;
        if ($fragLength <= $upperBound)	
        {
          $avgLength += $fragLength; 
          $cnt++;
          if ($i < $#posTags)
          {
            if ($posTags[$i+1] < $negTags[$retIndex] )
            {
              $sub++;
            }
            else
            {
              $avgLength -= ($sub+1)*($fragLength/2);
              $sub = 0;
            }
	  }
        }
      }
    }
  }
  undef @posTags;
  undef @negTags;
  $avgLength /= $cnt;
  return int(2*$avgLength+1.0);
}


sub binarySearch
{
  my ($begin, $query, $starts); 
  ($begin, $query, $starts) = @_;

  return 0 if ($query < $starts->[0] && $begin == 1);
  return -1 if ($query < $starts->[0] && $begin == 0);
  return $#$starts if ($query > $starts->[$#$starts] && $begin == 0);
  return $#$starts+1 if ($query > $starts->[$#$starts] && $begin == 1);

  my ($left) = 0;
  my ($right) = $#$starts;
  my ($center);
  my ($prevCenter) = -1;
  while(1)
  {
    $center = int( ($left + $right)/2 );
    if ($query == $starts->[$center])
    {
      if ($begin == 1)
      {
        while($center > 0 && $starts->[$center] == $starts->[$center-1]) 
          { $center =  $center - 1; }
      }
      else
      {
        while($center < $#$starts && $starts->[$center] == $starts->[$center+1]) 
          { $center =  $center + 1; }
      }
      return $center;
    }
    if ($center == $prevCenter)
    {
      return $right if ($begin == 1);
      return $right-1 if ($begin == 0);
    }
    $right = $center if ($query < $starts->[$center]);
    $left = $center if ($query > $starts->[$center]);
    $prevCenter = $center;
  }
}


sub binarySearchReverseSorted 
{
  my ($begin, $query, $starts); 
  ($begin, $query, $starts) = @_;

  ### Handling cases when query is out of bounds
  return 0 if ($query > $starts->[0] && $begin == 1);
  return -1 if ($query > $starts->[0] && $begin == 0);
  return $#$starts if ($query < $starts->[$#$starts] && $begin == 0);
  return $#$starts+1 if ($query < $starts->[$#$starts] && $begin == 1);

  my ($left) = 0;
  my ($right) = $#$starts;
  my ($center);
  my ($prevCenter) = -1;
  while(1)
  {
    $center = int( ($left + $right)/2 );
    if ($query == $starts->[$center])
    {
      if ($begin == 1)
      {
        while($center > 0 && $starts->[$center] == $starts->[$center-1])
          { $center =  $center - 1; }
      }
      else
      {
        while($center < $#$starts && $starts->[$center] == $starts->[$center+1])
          { $center =  $center + 1; }
      }
      return $center;
    }
    if ($center == $prevCenter)
    {
      return $right if ($begin == 1);
      return $right-1 if ($begin == 0);
    }
    $right = $center if ($query > $starts->[$center]);
    $left = $center if ($query < $starts->[$center]);
    $prevCenter = $center;
  }
}

sub findBindingSites
{
  my ($chrom, @starts, @myPos, @myNeg, @myPosControl, @myNegControl, %startsHash, %posStartsHash, %negStartsHash, $i, $tagCnt, $posTagCnt, $negTagCnt, $carryOver, $posCarryOver, $negCarryOver, $anchor, $prevAnchor, $nextAnchor, $startAnchor, $stopAncher, $prevBar, %windowCount, %posWindowCount, %negWindowCount, $key, $totalTagCnt, $posTotalTagCnt, $negTotalNegCnt, $numPosCount, $numNegCount, $leftEdgeCoord, $rightEdgeCoord, $leftIndex, $rightIndex, $evid, $evidControl, $firstTagIndex, $evidCnt, @negTags, @posTags, $bsitesMissedPos, $bsitesMissedNeg, $foldEnrichment, $numPosCountControl, $numNegCountControl, $leftEdgeCoordControl, $rightEdgeCoordControl, $leftIndexControl, $rightIndexControl);

  open(BED, ">$outputFile.tmp") || die "5 can't open file $outputFile.tmp\n";
  foreach $chrom (sort keys%bedStartsPos)
  {
    ### for every starting position, counting the number of tags in pos strand 

    @starts = split /\s/, $bedStartsPos{$chrom};
    @myPos = sort {$a <=> $b} @starts;
    %startsHash = ();	# key = start position regardless of strand; value = NET tags (pos - neg)
    %posStartsHash = ();	# key = start position if on pos strand; value = no of pos tags at this position
    %negStartsHash = ();	# key = start position if on neg strand; value = no of neg tags at this pos
    for ($i = 0; $i <= $#starts; $i++)
    {
      if ( !exists($startsHash{$starts[$i]}) )
      {
        $startsHash{$starts[$i]} = 1;
        $posStartsHash{$starts[$i]} = 1;
        $negStartsHash{$starts[$i]} = 0;
      }
      else
      {
        $startsHash{$starts[$i]} += 1;
        $posStartsHash{$starts[$i]} += 1;
      }
    }
  
    ### for every starting position, "un"counting the number of tags in neg strand 
  
    @starts = split /\s/, $bedStartsNeg{$chrom};
    @myNeg = sort {$a <=> $b} @starts;
    for ($i = 0; $i <= $#starts; $i++)
    {
      if ( !exists($startsHash{$starts[$i]}) )
      {
        $startsHash{$starts[$i]} = -1;
        $negStartsHash{$starts[$i]} = 1;
        $posStartsHash{$starts[$i]} = 0 if ( !exists($posStartsHash{$starts[$i]}) );
  					### if clause above is unnecessary
      }
      else
      {
        $startsHash{$starts[$i]} -= 1;
        $negStartsHash{$starts[$i]} += 1;
      }
    }

    ### CONTROL: for every starting position, counting the number of tags in pos strand

    if ( exists($Options{b}) )
    {
      @myPosControl = split /\s/, $bedStartsPosControl{$chrom};
      @myPosControl = sort {$a <=> $b} @myPosControl;
      @myNegControl = split /\s/, $bedStartsNegControl{$chrom};
      @myNegControl = sort {$a <=> $b} @myNegControl;
    }
  
    $tagCnt = 0;
    $posTagCnt = 0;
    $negTagCnt = 0;
    $carryOver = 0;
    $posCarryOver = 0;
    $negCarryOver = 0;
    $anchor = 0; 	### anchor is where the tag difference count will be plotted
    $prevAnchor = 0;
    %windowCount = ();	### these are 20 bp window
    %posWindowCount = ();
    %negWindowCount = ();
    foreach $key (sort { $a <=> $b } keys%startsHash)
    {
      if ($anchor == 0)	### the following block is entered only at the beginning
      {
        $anchor = $key - ($key%$interval) + 1 if ($key%$interval > 0); 
        $anchor = $key - $interval + 1 if ($key%$interval == 0); 
        $tagCnt = 0;
        $posTagCnt = 0;
        $negTagCnt = 0;
      }
  
      if ($key <= $anchor + $interval - 1)	### key < next anchor
      {
        $tagCnt += $startsHash{$key};
        $posTagCnt += $posStartsHash{$key};
        $negTagCnt += $negStartsHash{$key};
      }
      else
      {
        $totalTagCnt = $tagCnt;
        $posTotalTagCnt = $posTagCnt;
        $negTotalTagCnt = $negTagCnt;
        $totalTagCnt += $carryOver if ( $prevAnchor + $interval == $anchor );
        $posTotalTagCnt += $posCarryOver if ( $prevAnchor + $interval == $anchor );
        $negTotalTagCnt += $negCarryOver if ( $prevAnchor + $interval == $anchor );
        $windowCount{$anchor} = $totalTagCnt;
        $posWindowCount{$anchor} = $posTotalTagCnt;
        $negWindowCount{$anchor} = $negTotalTagCnt;
        $carryOver = $tagCnt;
        $posCarryOver = $posTagCnt;
        $negCarryOver = $negTagCnt;
        $prevAnchor = $anchor;
  
        $anchor = $key - ($key%$interval) + 1 if ($key%$interval > 0); 
        $anchor = $key - $interval + 1 if ($key%$interval == 0); 
        if ($anchor > $prevAnchor + $interval)
        {
          $nextAnchor = $prevAnchor + $interval;
          $windowCount{$nextAnchor} = $tagCnt;		### tagCnt or carryOver: both same
          $posWindowCount{$nextAnchor} = $posTagCnt;	### posTagCnt or posCarryOver: both same
          $negWindowCount{$nextAnchor} = $negTagCnt;	### negTagCnt or negCarryOver: both same
        }
        $tagCnt = $startsHash{$key};	
        $posTagCnt = $posStartsHash{$key};
        $negTagCnt = $negStartsHash{$key};
      }
    }
  
    $totalTagCnt = $tagCnt;
    $posTotalTagCnt = $posTagCnt;
    $negTotalTagCnt = $negTagCnt;
    $totalTagCnt += $carryOver if ( $prevAnchor + $interval == $anchor );
    $posTotalTagCnt += $posCarryOver if ( $prevAnchor + $interval == $anchor );
    $negTotalTagCnt += $negCarryOver if ( $prevAnchor + $interval == $anchor );
    $windowCount{$anchor} = $totalTagCnt;
    $posWindowCount{$anchor} = $posTotalTagCnt;
    $negWindowCount{$anchor} = $negTotalTagCnt;
    $nextAnchor = $anchor + $interval;
    $windowCount{$nextAnchor} = $tagCnt;
    $posWindowCount{$nextAnchor} = $posTagCnt;
    $negWindowCount{$nextAnchor} = $negTagCnt;
   
    #############################################################
    ### PRINTING OUT THE BINDING SITES FOUND ON THIS CHROMOSOME
    #############################################################
    $prevAnchor = -1;
    $startAnchor = -1;
    $stopAnchor = -1;
    $prevBar = 0;	### 1 if pos, 0 if zero, -1 if neg 
    foreach $anchor (sort {$a <=> $b} keys%windowCount)
    {
      if ($windowCount{$anchor} > 0)	### positive NET tag count
      {
        if ($prevBar == 0)  	### enters block for the first time
        {
          $prevBar = 1;
  	  $startAnchor = $anchor;
        }
        elsif ($prevBar == 1)	### pos to pos transition
        {
  	  $startAnchor = $anchor;
        }
        elsif ($prevBar == -1)	### neg to pos transition
        {
          $bindingSiteLocation = $startAnchor + ($stopAnchor - $startAnchor)/2;	### avg: same as (s+e)/2
          $bindingSiteStart = $startAnchor - ($interval/2);
          $bindingSiteEnd = $stopAnchor + ($interval/2);

          ### COUNTING EVIDENCE  

  	  $numPosCount = 0; ### number of evidences on the positive (left) side of bsite
  	  $numNegCount = 0; ### number of evidences on the negative (right) side of bsite
  
          $rightEdgeCoord= $bindingSiteLocation + $windowLookout;
          $leftIndex = binarySearch(1, $bindingSiteLocation, \@myNeg);
          $rightIndex = binarySearch(0, $rightEdgeCoord, \@myNeg);
          $numNegCount = $rightIndex - $leftIndex + 1;

          $leftEdgeCoord = $bindingSiteLocation - $windowLookout;
          $leftIndex = binarySearch(1, $leftEdgeCoord, \@myPos);
          $rightIndex = binarySearch(0, $bindingSiteLocation, \@myPos);
          $numPosCount = $rightIndex - $leftIndex + 1;
  
 	  $evid = $numPosCount+$numNegCount;
	
	  if ($numPosCount >= $numReadsOnEachSide && $numNegCount >= $numReadsOnEachSide && $evid >= $evidNecessary)
  	  {
            if ( exists($Options{b}) )
	    {
  	      $numPosCountControl = 0; ### number of evidences on the positive (left) side of bsite coordinate in the control set
  	      $numNegCountControl = 0; ### number of evidences on the negative (right) side of bsite coordinate in the control set
  
              $rightEdgeCoordControl = $bindingSiteLocation + $windowLookout;
              $leftIndexControl = binarySearch(1, $bindingSiteLocation, \@myNegControl);
              $rightIndexControl = binarySearch(0, $rightEdgeCoordControl, \@myNegControl);
              $numNegCountControl = $rightIndexControl - $leftIndexControl + 1;
    
              $leftEdgeCoordControl = $bindingSiteLocation - $windowLookout;
              $leftIndexControl = binarySearch(1, $leftEdgeCoordControl, \@myPosControl);
              $rightIndexControl = binarySearch(0, $bindingSiteLocation, \@myPosControl);
              $numPosCountControl = $rightIndexControl - $leftIndexControl + 1;
  
              $evidControl = $numPosCountControl+$numNegCountControl;
              if ($evidControl != 0)
	      {
	        $foldEnrichment = ($evid/$evidControl) * ($controlTagsSelected/$tagsSelected);
	      }
	      else
	      {
	        $foldEnrichment = $evid * ($controlTagsSelected/$tagsSelected);
	      }

              print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evid\t$foldEnrichment\n"; 
	    }
	    else
	    {
              print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evid\n"; 
	    }
  	  }
	  $prevBar = 1;
	  $startAnchor = $anchor;
        }
      }
      elsif ($windowCount{$anchor} < 0)
      {
        if ($prevBar == 1)	### pos to neg transition
        {
          $stopAnchor = $anchor;
          $prevBar = -1;
        }
        elsif ($prevBar == -1)	### neg to neg transition
        {
        }
      }
    }

    if ($prevBar == -1)		### ended on a negative bar, print out last binding site
    {
      $bindingSiteLocation = $startAnchor + ($stopAnchor - $startAnchor)/2;
      $bindingSiteStart = $startAnchor - ($interval/2);
      $bindingSiteEnd = $stopAnchor + ($interval/2);
  
      $numPosCount = 0;
      $numNegCount = 0;
  
      $rightEdgeCoord= $bindingSiteLocation + $windowLookout;
      $leftIndex = binarySearch(1, $bindingSiteLocation, \@myNeg);
      $rightIndex = binarySearch(0, $rightEdgeCoord, \@myNeg);
      $numNegCount = $rightIndex - $leftIndex + 1;
  
      $leftEdgeCoord = $bindingSiteLocation - $windowLookout;
      $leftIndex = binarySearch(1, $leftEdgeCoord, \@myPos);
      $rightIndex = binarySearch(0, $bindingSiteLocation, \@myPos);
      $numPosCount = $rightIndex - $leftIndex + 1;
  
      $evid = $numPosCount + $numNegCount;
      if ($numPosCount >= $numReadsOnEachSide && $numNegCount >= $numReadsOnEachSide && $evid >= $evidNecessary)
      {
        if ( exists($Options{b}) )
        {
          $numPosCountControl = 0; ### number of evidences on the positive (left) side of bsite coordinate in the control set
          $numNegCountControl = 0; ### number of evidences on the negative (right) side of bsite coordinate in the control set

          $rightEdgeCoordControl = $bindingSiteLocation + $windowLookout;
          $leftIndexControl = binarySearch(1, $bindingSiteLocation, \@myNegControl);
          $rightIndexControl = binarySearch(0, $rightEdgeCoordControl, \@myNegControl);
          $numNegCountControl = $rightIndexControl - $leftIndexControl + 1;

          $leftEdgeCoordControl = $bindingSiteLocation - $windowLookout;
          $leftIndexControl = binarySearch(1, $leftEdgeCoordControl, \@myPosControl);
          $rightIndexControl = binarySearch(0, $bindingSiteLocation, \@myPosControl);
          $numPosCountControl = $rightIndexControl - $leftIndexControl + 1;

          $evidControl = $numPosCountControl+$numNegCountControl;
          if ($evidControl != 0)
	  {
	    $foldEnrichment = ($evid/$evidControl) * ($controlTagsSelected/$tagsSelected);
	  }
	  else
	  {
	    $foldEnrichment = $evid * ($controlTagsSelected/$tagsSelected);
	  }

          print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evid\t$foldEnrichment\n";
        }
        else
        {
          print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evid\n";
        }
      }
    }
  }
  close BED;
  
  ############################################################################
  ### Finding binding sites with tags only in one direction
  ### which would have been missed by above two-side evid procedure
  ############################################################################
  
  if ($oneside == 1)
  {
  
    open(BED, ">>$outputFile.tmp") || die "can't open file $outputFile.tmp\n";
    $bsitesMissedPos = 0;
    foreach $chrom (sort keys%bedStartsPos)
    {
      @posTags = split /\s/, $bedStartsPos{$chrom};
      @posTags = sort {$a <=> $b} @posTags;
      @negTags = split /\s/, $bedStartsNeg{$chrom};
      @negTags = sort {$a <=> $b} @negTags;
      $firstTagIndex = 0;
      $evidCnt = 0;
      for ($i = 0; $i <= $#posTags; $i++)
      {
        $nextNegTagIndex = binarySearch(1, $posTags[$i], \@negTags);
        if ($negTags[$nextNegTagIndex] - $posTags[$i] > 2*$fragmentLength) 
		### next NEG tag is at least two fragment length more away
        {
          $evidCnt++;
          while($posTags[$i] - $posTags[$firstTagIndex] + 1 > $fragmentLength)
          {
            $firstTagIndex++;
            $evidCnt--;
          }
    
          if (  ($i<$#posTags && $posTags[$i+1]-$posTags[$i] > $fragmentLength) || ($i==$#posTags)  )       
          	### if next POS tag is at least one fragment length away
          {
            if ($evidCnt >= $evidNecessary && $posTags[$i] != $posTags[$firstTagIndex]) ### missed bsite
            {
              $bsitesMissedPos++;
    	      $bindingSiteLocation = $posTags[$i] + (int($fragmentLength/2));
	      $bindingSiteStart = $posTags[$i];
	      $bindingSiteEnd = $posTags[$i] + $fragmentLength;
              if ( exists($Options{b}) )
              {
                $numPosCountControl = 0; ### number of evidences on the positive (left) side of bsite coordinate in the control set
                $numNegCountControl = 0; ### number of evidences on the negative (right) side of bsite coordinate in the control set

                $rightEdgeCoordControl = $posTags[$i] + $windowLookout;
                $leftIndexControl = binarySearch(1, $posTags[$i], \@myNegControl);
                $rightIndexControl = binarySearch(0, $rightEdgeCoordControl, \@myNegControl);
                $numNegCountControl = $rightIndexControl - $leftIndexControl + 1;

                $leftEdgeCoordControl = $posTags[$i] - $windowLookout;
                $leftIndexControl = binarySearch(1, $leftEdgeCoordControl, \@myPosControl);
                $rightIndexControl = binarySearch(0, $posTags[$i], \@myPosControl);
                $numPosCountControl = $rightIndexControl - $leftIndexControl + 1;

                $evidControl = $numPosCountControl+$numNegCountControl;
                if ($evidControl != 0)
	        {
	          $foldEnrichment = ($evidCnt/$evidControl) * ($controlTagsSelected/$tagsSelected);
	        }
	        else
	        {
	          $foldEnrichment = $evidCnt * ($controlTagsSelected/$tagsSelected);
	        }

                print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evidCnt\t$foldEnrichment\n";
              }
	      else
	      {
                #print BED "$chrom\t$bindingSiteLocation\t$evidCnt\n";
                print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evidCnt\n";
	      }
            }
            else
            {
              $firstTagIndex = $i+1;
              $evidCnt = 0;
            }
          }
        }
        else
        {
          $firstTagIndex = $i+1;
          $evidCnt = 0;
        }
      }
    }
  
  
    $bsitesMissedNeg = 0;
    foreach $chrom (sort keys%bedStartsNeg)
    {
      @posTags = split /\s/, $bedStartsPos{$chrom};
      @posTags = sort {$b <=> $a} @posTags;	### reverse sorting
      @negTags = split /\s/, $bedStartsNeg{$chrom};
      @negTags = sort {$b <=> $a} @negTags;	### reverse sorting
      $firstTagIndex = 0;
      $evidCnt = 0;
      for ($i = 0; $i <= $#negTags; $i++)
      {
        $nextPosTagIndex = binarySearchReverseSorted(1, $negTags[$i], \@posTags);
        if ($negTags[$i] - $posTags[$nextPosTagIndex] > 2*$fragmentLength) ### next POS tag is 2*200bp or more away
        {
          $evidCnt++;
          while($negTags[$firstTagIndex] - $negTags[$i] + 1 > $fragmentLength)
          {
            $firstTagIndex++;
            $evidCnt--;
          }
          if (  ($i<$#negTags && $negTags[$i]-$negTags[$i+1] > $fragmentLength) || ($i==$#negTags)  )
                                    ### if next NEG tag is one fragment length away
          {
            if ($evidCnt >= $evidNecessary && $negTags[$i] != $negTags[$firstTagIndex]) ### missed bsite
            {
              $bsitesMissedNeg++;
              $bindingSiteLocation = $negTags[$i] - (int($fragmentLength/2));
	      $bindingSiteStart = $negTags[$i] - $fragmentLength;
	      $bindingSiteEnd = $negTags[$i];

              if ( exists($Options{b}) )
              {
                $numPosCountControl = 0; ### number of evidences on the positive (left) side of bsite coordinate in the control set
                $numNegCountControl = 0; ### number of evidences on the negative (right) side of bsite coordinate in the control set

                $rightEdgeCoordControl = $negTags[$i] + $windowLookout;
                $leftIndexControl = binarySearch(1, $negTags[$i], \@myNegControl);
                $rightIndexControl = binarySearch(0, $rightEdgeCoordControl, \@myNegControl);
                $numNegCountControl = $rightIndexControl - $leftIndexControl + 1;

                $leftEdgeCoordControl = $negTags[$i] - $windowLookout;
                $leftIndexControl = binarySearch(1, $leftEdgeCoordControl, \@myPosControl);
                $rightIndexControl = binarySearch(0, $negTags[$i], \@myPosControl);
                $numPosCountControl = $rightIndexControl - $leftIndexControl + 1;

                $evidControl = $numPosCountControl+$numNegCountControl;
                if ($evidControl != 0)
	        {
	          $foldEnrichment = ($evidCnt/$evidControl) * ($controlTagsSelected/$tagsSelected);
	        }
	        else
	        {
	          $foldEnrichment = $evidCnt * ($controlTagsSelected/$tagsSelected);
	        }

                print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evidCnt\t$foldEnrichment\n";
              }
	      else
	      {
                #print BED "$chrom\t$bindingSiteLocation\t$evidCnt\n";
                print BED "$chrom\t$bindingSiteStart\t$bindingSiteEnd\t$evidCnt\n";
	      }
            }
            else
            {
              $firstTagIndex = $i+1;
              $evidCnt = 0;
            }
          }
        }
        else
        {
          $firstTagIndex = $i+1;
          $evidCnt = 0;
        }
      }
    }
    close BED;
  
  }	### end of one-sided tags

  undef @starts; 
  undef @myPos;
  undef @myNeg;
  undef @myPosControl;
  undef @myNegControl;
  undef %startsHash;
  undef %posStartsHash;
  undef %negStartsHash;
  undef %windowCount;
  undef %posWindowCount;
  undef %negWindowCount;
}

sub pValueDistribution
{
  my($t, %chromLengths, %chromStart, %chromEnd, $effectiveGenomeLength, $chrom, @posReads, @negReads, @myPosFore, @myNegFore, @myPosBack, @myNegBack, $trials, $i, $randCoord, $numPosCountFore, $numNegCountFore, $rightEdgeCoordFore, $leftIndexFore, $rightIndexFore, $numNegCountFore, $leftEdgeCoordFore, $leftIndexFore, $rightIndexFore, $numPosCountFore, $evidFore, $numPosCountBack, $numNegCountBack, $rightEdgeCoordBack, $leftIndexBack, $rightIndexBack, $numNegCountBack, $leftEdgeCoordBack, $leftIndexBack, $rightIndexBack, $numPosCountBack, $evidBack, $foldRatio);
  %chromLengths = ();
  %chromStart = ();
  %chromEnd = ();
  $effectiveGenomeLength = 0;
  foreach $chrom (sort keys%bedStartsPos)
  {
    @posReads = split /\s/, $bedStartsPos{$chrom};
    @posReads = sort {$a <=> $b} @posReads;
    @negReads = split /\s/, $bedStartsNeg{$chrom};
    @negReads = sort {$a <=> $b} @negReads;
    $chromStart{$chrom} = $posReads[0];
    $chromStart{$chrom} = $negReads[0] if ($negReads[0] < $posReads[0]);
    $chromEnd{$chrom} = $posReads[$#posReads];
    $chromEnd{$chrom} = $negReads[$#negReads] if ($negReads[$#negReads] > $posReads[$#posReads]);
    $chromLengths{$chrom} = $chromEnd{$chrom} - $chromStart{$chrom} + 1;
    $effectiveGenomeLength += $chromLengths{$chrom};
  }
  foreach $chrom (sort keys%chromLengths)
  {
    @myPosFore = split /\s/, $bedStartsPos{$chrom};
    @myPosFore = sort {$a <=> $b} @myPosFore;
    @myNegFore = split /\s/, $bedStartsNeg{$chrom};
    @myNegFore = sort {$a <=> $b} @myNegFore;
    @myPosBack = split /\s/, $bedStartsPosControl{$chrom};
    @myPosBack = sort {$a <=> $b} @myPosBack;
    @myNegBack = split /\s/, $bedStartsNegControl{$chrom};
    @myNegBack = sort {$a <=> $b} @myNegBack;

    $trials = $numTrials * $chromLengths{$chrom} / $effectiveGenomeLength;
    for ($i = 0; $i < $trials; $i++)
    {
      $randCoord = int(rand($chromLengths{$chrom})) + $chromStart{$chrom};

      $numPosCountFore = 0; ### number of evidences on the positive (left) side of random coordinate in the control set
      $numNegCountFore = 0; ### number of evidences on the negative (right) side of random coordinate in the control set
      $rightEdgeCoordFore = $randCoord+ $windowLookout;
      $leftIndexFore = binarySearch(1, $randCoord, \@myNegFore);
      $rightIndexFore = binarySearch(0, $rightEdgeCoordFore, \@myNegFore);
      $numNegCountFore = $rightIndexFore - $leftIndexFore + 1;
      $leftEdgeCoordFore = $randCoord- $windowLookout;
      $leftIndexFore = binarySearch(1, $leftEdgeCoordFore, \@myPosFore);
      $rightIndexFore = binarySearch(0, $randCoord, \@myPosFore);
      $numPosCountFore= $rightIndexFore- $leftIndexFore + 1;
      $evidFore = $numPosCountFore+$numNegCountFore;

      $numPosCountBack = 0; ### number of evidences on the positive (left) side of random coordinate in the control set
      $numNegCountBack = 0; ### number of evidences on the negative (right) side of random coordinate in the control set
      $rightEdgeCoordBack = $randCoord+ $windowLookout;
      $leftIndexBack = binarySearch(1, $randCoord, \@myNegBack);
      $rightIndexBack = binarySearch(0, $rightEdgeCoordBack, \@myNegBack);
      $numNegCountBack = $rightIndexBack - $leftIndexBack + 1;
      $leftEdgeCoordBack = $randCoord- $windowLookout;
      $leftIndexBack = binarySearch(1, $leftEdgeCoordBack, \@myPosBack);
      $rightIndexBack = binarySearch(0, $randCoord, \@myPosBack);
      $numPosCountBack= $rightIndexBack- $leftIndexBack + 1;
      $evidBack = $numPosCountBack+$numNegCountBack;

      if ( $evidBack != 0 )
      {
        $foldRatio = $evidFore/$evidBack;
      }
      else
      {
        $foldRatio = $evidFore;
      }
      push(@foldArray, $foldRatio);
    }
  }
  undef %chromLengths;
  undef %chromStart;
  undef %chromEnd;
  undef @posReads;
  undef @negReads;
  undef @myPosFore;
  undef @myNegFore;
  undef @myPosBack;
  undef @myNegBack;
}

sub get_pValue
{
  my($foldArg) = @_;
  my($pValueIndex) = binarySearchReverseSorted(1, $foldArg, \@foldArray);
  my($pValueLocal) = $pValueIndex/$#foldArray;
  return $pValueLocal;
}


sub printSummary
{
  open(OUT, ">$outputFile") || die "can't open file $outputFile\n";

  print OUT "======================================================================\n";
  print OUT "SISSRs: A tool to identify binding sites from ChIP-Seq data\n";
  print OUT "======================================================================\n";
  print OUT "SISSRs version 1.4 (Release date: Mon, 24 November 2008)\n\n";
  print OUT "For further information on how to interpret these results or to get\n";
  print OUT "a recent version of the SISSRs software, please visit\n";
  print OUT "http://sissrs.rajajothi.com\n";
  print OUT "======================================================================\n\n\n";

  print OUT "======================================================================\n";
  print OUT "REFERENCE\n";
  print OUT "======================================================================\n";
  print OUT "If you use this program in your research, please cite:\n\n";
  print OUT "Raja Jothi, Suresh Cuddapah, Artem Barski, Kairong Cui, Keji Zhao\n";
  print OUT "Genome-wide identification of in vivo protein-DNA binding sites\n";
  print OUT "from ChIP-Seq Data. Nucleic Acids Research, 36(16):5221-31 (2008)\n";
  print OUT "======================================================================\n\n\n";


  print OUT "======================================================================\n";
  print OUT "COMMAND LINE SUMMARY & ESTIMATED PARAMETERS\n";
  print OUT "======================================================================\n";
  print "\nData file (i)                          : $bedFileName\n" if ($progress);
  print OUT "Data file (i)                           : $bedFileName\n";
  print "Number of tags in the data file         : $ctr\n" if ($progress);
  print OUT "Number of tags in the data file         : $ctr\n";
  if ( exists($Options{q}) )
  {
    printf "Number of tags selected for analysis\n  after ignoring tags that overlap\n  regions listed in $repeatsFileName        : $tagsSelected (%5.2f%)\n", $fracSelected if ($progress);
    printf OUT "Number of tags selected for analysis\n  after ignoring tags that overlap\n  regions listed in $repeatsFileName        : $tagsSelected (%5.2f%)\n", $fracSelected;
  }
  else
  {
    printf "Number of tags selected for analysis    : $tagsSelected (%5.2f%)\n", $fracSelected if ($progress);
    printf OUT "Number of tags selected for analysis    : $tagsSelected (%5.2f%)\n", $fracSelected;
  }
  printf "Tags mapped to sense strand             : $numPositive (%5.2f%)\n", $posFrac if ($progress);
  printf OUT "Tags mapped to sense strand             : $numPositive (%5.2f%)\n", $posFrac;
  printf "Tags mapped to anti-sense strand        : $numNegative (%5.2f%)\n", $negFrac if ($progress);
  printf OUT "Tags mapped to anti-sense strand        : $numNegative (%5.2f%)\n", $negFrac;
  print "Background model (Negative Control)     : $Options{b}\n" if ($progress && exists($Options{b}) );
  print OUT "Background model (Negative Control)     : $Options{b}\n" if ( exists($Options{b}) );
  print "Number of tags in the control file      : $controlTagsSelected\n" if ($progress && exists($Options{b}) );
  print OUT "Number of tags in the control file      : $controlTagsSelected\n" if ( exists($Options{b}) );
  print "Background model (Negative Control)     : Poisson\n" if ($progress && !exists($Options{b}) );
  print OUT "Background model (Negative Control)     : Poisson\n" if ( !exists($Options{b}) );
  print "Genome length (s)                       : $Options{s}\n" if ($progress);
  print OUT "Genome length (s)                       : $Options{s}\n";
  printf "Fraction of genome mappable by reads (m): %3.2f\n", $mappableFraction if ( !exists($Options{q})  && $progress);
  printf OUT "Fraction of genome mappable by reads (m): %3.2f\n", $mappableFraction if ( !exists($Options{q}));
  printf "Fraction of genome mappable by reads (m): %3.2f-%8.6f\n", $mappableFraction, $ignoreFraction if ( exists($Options{q}) && $progress);
  printf OUT "Fraction of genome mappable by reads (m): %3.2f-%8.6f\n", $mappableFraction, $ignoreFraction if ( exists($Options{q}));
  print "Effective genome length (s*m)           : $genomeLength\n" if ($progress);
  print OUT "Effective Genome length  (s*m)          : $genomeLength\n";
  print "E-value (e)                             : $eValue\n" if ($progress && exists($Options{b}) );
  print OUT "E-value (e)                             : $eValue\n" if ( exists($Options{b}) );
  print "P-value (p)                             : $fold_pValue\n" if ($progress && exists($Options{b}) );
  print OUT "P-value (p)                             : $fold_pValue\n" if ( exists($Options{b}) );
  print "False discovery rate (d)                : $fdr\n" if ($progress && !exists($Options{b}) );
  print OUT "False discovery rate (d)                : $fdr\n" if ( !exists($Options{b}) );
  print "Scanning window size (w)                : $interval\n" if ($progress);
  print OUT "Scanning window size (w)                : $interval\n";
  print "Average DNA fragment length (f)         : $fragmentLength\n" if ($progress);
  print OUT "Average DNA fragment length (f)         : $fragmentLength\n";
  print "Minimum number of 'directional' tags\n  required on each side of the inferred\n  binding site (E)                      : $numReadsOnEachSide\n" if ($progress);
  print OUT "Minimum number of 'directional' tags\n  required on each side of the inferred\n  binding site (E)                      : $numReadsOnEachSide\n";
  print "Keep one tag per genomic position (a)   : NO\n" if ($keepOneRead == 0 && $progress);
  print "Keep one tag per genomic position (a)   : YES\n" if ($keepOneRead == 1 && $progress);
  print OUT "Keep one tag per genomic position (a)   : NO\n" if ($keepOneRead == 0);
  print OUT "Keep one tag per genomic position (a)   : YES\n" if ($keepOneRead == 1);
  print "Also reports binding sites supported\n  only by reads mapped to either sense\n  or anti-sense strand (u)              : NO\n" if ($oneside == 0 && $progress);
  print OUT "Also reports binding sites supported\n  only by reads mapped to either sense\n  or anti-sense strand (u)              : NO\n" if ($oneside == 0);
  print "Also reports binding sites supported\n  only by reads mapped to either sense\n  or anti-sense strand (u)              : YES\n" if ($oneside == 1 && $progress);
  print OUT "Also reports binding sites supported\n  only by reads mapped to either sense\n  or anti-sense strand (u)              : YES\n" if ($oneside == 1);

  print OUT "======================================================================\n\n\n";

  close OUT;

}

sub printUsage
{
  print "\n=======================================================================\n"; 
  print "\nUSAGE:\n\n sissrs.pl -i <input-file> -o <output-file> -s <genome-size> [OPTIONS]\n\n"; 
  print " [-i <file>]\tinput file containing tags/reads in BED format\n";
  print " [-o <file>]\toutput file into which results will be stored\n";
  print " [-s <int>]\tgenome size (number of bases or nucleotides)\n";
  print " [-a]\t\tonly one read is kept if multiple reads align to the\n\t\tsame genomic coordinate (minimizes amplification bias)\n";
  print " [-F <int>]\taverage length of DNA fragments that were sequenced\n\t\t(default: estimated from reads)\n";
  print " [-D <real>]\tfalse discovery rate (default: 0.001) if random\n\t\tbackground model based on Poisson probabilities need\n\t\tto be used as control (also check option -b below)\n";
  print " [-b <real>]\tbackground file containing tags in BED format to be\n\t\tused as control; -e and -p can be set to desired\n\t\tvalues to control specificity and specificity, resp.\n";
  print " [-e <real>]\te-value (>=0); it is the number of binding sites\n\t\tone might expect to infer by chance (default: 10);\n\t\tthis option is irrelevant if -b option is NOT used\n";
  print " [-p <real>]\tp-value threshold for fold enrichment of ChIP tags\n\t\tat a binding site location compared to that at the\n\t\tsame location in the control data (default: 0.001);\n\t\tthis option is irrelevant if -b option is NOT used\n";
  print " [-m <int>]\tfraction of genome mappable by reads (default: 0.8\n\t\tfor hg18, assuming ELAND was used to map the reads;\n\t\tcould be different for different organisms and other\n\t\talignment algorithms)\n";
  print " [-w <int>]\tscanning window size (even number > 1), which\n\t\tcontrols for noise (default: 20)\n";
  print " [-E <int>]\tmin number of 'directional' reads required on each\n\t\tside of the inferred binding site (>0); (default: 2)\n";
  print " [-L <int>]\tupper-bound on the DNA fragment length (default: 500)\n";
  print " [-q <file>]\tfile containing genomic regions to exclude; reads\n\t\tmapped to these regions will be ignored; file\n\t\tformat: 'chr startCoord endCoord'\n";
  print " [-t]\treports each binding site as a single genomic\n\t\tcoordinate (transition point t in Fig 1A)\n";
  print " [-r]\t\treports each binding site as an X-bp binding region\n\t\tcentered on inferred binding coordinate; X denotes\n\t\tthe distance from the start of the right-most red\n\t\tbar (see Fig 1A in the manuscript; lower-left)\n\t\tto the end of the left-most blue bar surrounding the\n\t\tactual binding site (transition point t in Fig 1A)\n";
  print " [-c]\t\tsame as the -r option, except that it reports binding\n\t\tsites that are clustered within F bp of each other as\n\t\ta single binding region; this is the default option\n";
  print " [-u]\t\t(also) reports binding sites supported only by reads\n\t\tmapped to either sense or anti-sense strand; this\n\t\toption will recover binding sites whose sense or\n\t\tanti-sense reads were not mapped for some reason\n\t\t(e.g., falls in unmappable/repetitive regions)\n";
  print " [-x]\t\tdo not print progress report (default: prints report)\n";
 print  "\n=======================================================================\n\n"; 
}

