#!/usr/bin/env perl

# Author: Karl Gemayel
# Created: November 30, 2016
# 
# Some options edited bt AL.
#
# Run the GeneMarkS-2 gene-finder.

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;

# get script name
my $scriptName = basename($0);

# get path of current script
my $scriptPath = abs_path($0);
$scriptPath =~ /^(.*)\//;
$scriptPath = $1;

my $trainer = "$scriptPath/biogem";        # training of parameters 
my $predictor = "$scriptPath/gmhmmp2";      # predicting genes

my $comparePrediction = "$scriptPath/compp";    # compare prediction files to check for convergence

# ------------------------------ #
#      Modes for iterations      #
# ------------------------------ #
my $modeNoMotif = "no-motif"; 
my $modeGroupAStep1  = "group-a1";
my $modeGroupAStep2  = "group-a2";
my $modeGroupB  = "group-b";
my $modeGroupC  = "group-c";
my $modeGroupD  = "group-d";
my $modeGroupE  = "group-e";
my @validIterationModes = ($modeNoMotif, $modeGroupAStep1, $modeGroupAStep2, $modeGroupB, $modeGroupC, $modeGroupD, $modeGroupE);




# ------------------------------ #
#    Default Variable Values     #
# ------------------------------ # 
my $D_GENETIC_CODE                      = 11                                ;
my $D_FNOUTPUT                          = "gms2.lst"                        ;
my $D_FORMAT_OUTPUT                     = "lst"                             ;
my $D_MGMTYPE                           = "auto"                            ;
my $D_PROM_WIDTH_A                      = 12                                ;
my $D_PROM_WIDTH_B                      = 6                                 ;
my $D_RBS_WIDTH                         = 6                                 ;
my $D_PROM_UPSTR_LEN_A                  = 40                                ;
my $D_PROM_UPSTR_LEN_B                  = 20                                ;
my $D_RBS_UPSTR_LEN                     = 20                                ;
my $D_SPACER_SCORE_THRESH_A             = 0.1                               ;
my $D_SPACER_SCORE_THRESH_B             = 0.25                              ;
my $D_SPACER_DIST_THRESH                = 14                                ;
my $D_SPACER_WINDOW_SIZE                = 1                                 ;
my $D_16S                               = "TAAGGAGGTGA"                     ;
my $D_MIN_MATCH_16S                     = 4                                 ;
my $D_MIN_MATCH_RBS_PROM                = 3                                 ;
my $D_MIN_FRAC_RBS_16S_MATCH            = 0.5                               ;
my $D_UPSTR_SIG_LENGTH                  = 35                                ;
my $D_UPSTR_SIG_ORDER                   = 2                                 ;
my $D_MAX_ITER                          = 10                                ;
my $D_CONV_THRESH                       = 0.99                              ;
my $D_COD_ORDER                         = 5                                 ;
my $D_NONCOD_ORDER                      = 2                                 ;
my $D_START_CONTEXT_ORDER               = 2                                 ;
my $D_FGIO_DIST_THRESH                  = 25                                ;

# ------------------------------ #
#    Command-line variables      #
# ------------------------------ # 
my $fn_genome                                                               ;       # Name of file containing genome sequence
my $genomeType                                                              ;       # Type of genome: Options: archaea, bacteria, auto
my $geneticCode                         = $D_GENETIC_CODE                   ;       # Genetic code
my $fnoutput                            = $D_FNOUTPUT                       ;       # Name of final output file
my $formatOutput                        = $D_FORMAT_OUTPUT                  ;       # Format for output file
my $fnAA                                                                    ;       # amino acid sequences
my $fnNN                                                                    ;       # nucleotide sequences
my $gid                                 = ''                                ;       # change gene id labeling to contig_id with id local to contig
my $ncbi                                = ''                                ;       # parse genetic code from definition line of FASTA file

# Group-A
my $groupA_widthPromoter                = $D_PROM_WIDTH_A                   ;
my $groupA_widthRBS                     = $D_RBS_WIDTH                      ;
my $groupA_promoterUpstreamLength       = $D_PROM_UPSTR_LEN_A               ;
my $groupA_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $groupA_spacerScoreThresh            = $D_SPACER_SCORE_THRESH_A          ;
my $groupA_spacerDistThresh             = $D_SPACER_DIST_THRESH             ;
my $groupA_spacerWindowSize             = $D_SPACER_WINDOW_SIZE             ;

# Group-B
my $groupB_widthPromoter                = $D_PROM_WIDTH_B                   ;
my $groupB_widthRBS                     = $D_RBS_WIDTH                      ;
my $groupB_promoterUpstreamLength       = $D_PROM_UPSTR_LEN_B               ;
my $groupB_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $groupB_spacerScoreThresh            = $D_SPACER_SCORE_THRESH_B          ;
my $groupB_spacerWindowSize             = $D_SPACER_WINDOW_SIZE             ;
my $groupB_spacerDistThresh             = $D_SPACER_DIST_THRESH             ;
my $groupB_tail16S                      = $D_16S                            ;
my $groupB_minMatchToTail               = $D_MIN_MATCH_16S                  ;

# Group-C
my $groupC_widthRBS                     = $D_RBS_WIDTH                      ;
my $groupC_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $groupC_minMatchPromoterRBS          = $D_MIN_MATCH_RBS_PROM             ;
my $groupC_minMatchRBS16S               = $D_MIN_MATCH_16S                  ;

# Group-D
my $groupD_widthRBS                     = $D_RBS_WIDTH                      ;
my $groupD_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $groupD_percentMatchRBS              = $D_MIN_FRAC_RBS_16S_MATCH         ;
my $groupD_minMatchRBS16S               = $D_MIN_MATCH_16S                  ;
my $groupD_tail16S                      = $D_16S                            ;

# Group-E
my $groupE_widthRBS                     = $D_RBS_WIDTH                      ;
my $groupE_rbsUpstreamLength            = $D_RBS_UPSTR_LEN                  ;
my $groupE_upstreamSignatureLength      = $D_UPSTR_SIG_LENGTH               ;
my $groupE_upstreamSignatureOrder       = $D_UPSTR_SIG_ORDER                ;
my $groupE_tail16S                      = $D_16S                            ;

# Iteration control
my $MAX_ITER                            = $D_MAX_ITER                       ;       # number of max iterations in main cycle
my $CONV_THRESH                         = $D_CONV_THRESH                    ;       # convergence threshold
my $numIterWithoutRBS                   = 1                                 ;

# Model Hyperparameters
my $orderCod                            = $D_COD_ORDER                      ;       # order for coding model
my $orderNon                            = $D_NONCOD_ORDER                   ;       # order for noncoding model
my $scOrder                             = $D_START_CONTEXT_ORDER            ;       # start context order
my $fgioDistThresh                      = $D_FGIO_DIST_THRESH               ;
#
# Misc Variables
my $toMgmProb                           = 0.15                              ;
my $toNativeProb                        = 0.85                              ;
my $fixedNativeAtypicalProb                                                 ;
my $trainNonCodingOnFullGenome                                              ;
my $minAtypicalProb                     = 0.02                              ;
my $runMFinderWithoutSpacer                                                 ;
my $showAdvancedOptions                                                     ;            
my $mgmType = $D_MGMTYPE                                                    ;       # Type of MGM model: options: "bac, arc, auto"
my $verbose                                                                 ;       # verbose mode
my $keepAllFiles                                                            ;
my $forceGroup                                                              ;
my $fn_external                                                             ;       # External evidence in GFF format

# Parse command-line options
GetOptions (
    'seq=s'                                 =>  \$fn_genome,
    'genome-type=s'                         =>  \$genomeType,
    'gcode=i'                               =>  \$geneticCode,
    'output=s'                              =>  \$fnoutput,
    'format=s'                              =>  \$formatOutput,
    'faa=s'                                 =>  \$fnAA,
    'fnn=s'                                 =>  \$fnNN,
    'gid'                                   =>  \$gid,
    'ncbi'                                  =>  \$ncbi,
    # Group-A
    'group-a-width-promoter=i'              =>  \$groupA_widthPromoter,
    'group-a-width-rbs=i'                   =>  \$groupA_widthRBS,
    'group-a-promoter-upstream-length=i'    =>  \$groupA_promoterUpstreamLength,
    'group-a-rbs-upstream-length=i'         =>  \$groupA_rbsUpstreamLength,
    'group-a-spacer-score-thresh=f'         =>  \$groupA_spacerScoreThresh,
    'group-a-spacer-dist-thresh=i'          =>  \$groupA_spacerDistThresh,
    'group-a-spacer-window-size=i'          =>  \$groupA_spacerWindowSize,
    # Group-B
    'group-b-width-promoter=i'              =>  \$groupB_widthPromoter,
    'group-b-width-rbs=i'                   =>  \$groupB_widthRBS,
    'group-b-promoter-upstream-length=i'    =>  \$groupB_promoterUpstreamLength,
    'group-b-rbs-upstream-length=i'         =>  \$groupB_rbsUpstreamLength,
    'group-b-spacer-score-thresh=f'         =>  \$groupB_spacerScoreThresh,
    'group-b-spacer-window-size=i'          =>  \$groupB_spacerWindowSize,
    'group-b-tail-16s=s'                    =>  \$groupB_tail16S,
    'group-b-min-match-to-tail=i'           =>  \$groupB_minMatchToTail,
    # Group-C
    'group-c-width-rbs=i'                   =>  \$groupC_widthRBS,
    'group-c-rbs-upstream-length=i'         =>  \$groupC_rbsUpstreamLength,
    'group-c-min-match-promoter-rbs=i'      =>  \$groupC_minMatchPromoterRBS,
    # Group-D
    'group-d-width-rbs=i'                   =>  \$groupD_widthRBS,
    'group-d-rbs-upstream-length=i'         =>  \$groupD_rbsUpstreamLength,
    'group-d-percent-match-rbs=f'           =>  \$groupD_percentMatchRBS,
    # Group-E
    'group-e-width-rbs=i'                   =>  \$groupE_widthRBS,
    'group-e-rbs-upstream-length=i'         =>  \$groupE_rbsUpstreamLength,
    'group-e-upstream-signature-length=i'   =>  \$groupE_upstreamSignatureLength,
    'group-e-upstream-signature-order=i'    =>  \$groupE_upstreamSignatureOrder,
    'group-e-tail-16s=s'                    =>  \$groupE_tail16S,
    # Iteration control
    'max-iter=i'                            =>  \$MAX_ITER,
    'conv-thresh=f'                         =>  \$CONV_THRESH,
    # Model Hyperparameters: Orders
    'order-cod=i'                           =>  \$orderCod,
    'order-non=i'                           =>  \$orderNon,
    'order-sc=i'                            =>  \$scOrder,
    # Model Hyperparameters: lengths
    'fgio-dist-thresh=i'                    =>  \$fgioDistThresh,
    # Misc
    'fixed-native-atypical-prob'            =>  \$fixedNativeAtypicalProb,
    'train-noncoding-on-full-genome'        =>  \$trainNonCodingOnFullGenome,
    'min-atypical-prob=f'                   =>  \$minAtypicalProb,
    'run-mfinder-without-spacer'            =>  \$runMFinderWithoutSpacer,
    'v'                                     =>  \$verbose,
    'advanced-options'                      =>  \$showAdvancedOptions,
    'mgm-type=s'                            =>  \$mgmType,
    'keep-all-files'                        =>  \$keepAllFiles,
    'force-group=s'                         =>  \$forceGroup,
    'ext=s'                                 =>  \$fn_external,
);

Usage($scriptName) if (!defined $fn_genome or !defined $genomeType or !isValidGenomeType($genomeType));


if ($ncbi)
{
	$geneticCode = GetGeneticCodeFromFile( $fn_genome);
}

# variable that's set if group B is tested and Promoter and RBS matched
my $testGroupB_PromoterMatchedRBS;

# setup temporary file collection
my @tempFiles;

# create "single fasta format" from multifasta file
my $fnseq = "tmpseq.fna";
MultiToSingleFASTA($fn_genome, $fnseq);

# add temporary files
push @tempFiles, ($fnseq) unless $keepAllFiles;

my $mgmMod = "$scriptPath/mgm_$geneticCode.mod";        # name of MGM mod file (based on genetic code)
my $modForFinalPred = "tmp.mod";                        # used to keep a version of the model at every iteration 

my $alignmentInMFinder = "right";
if (defined $runMFinderWithoutSpacer) {
    $alignmentInMFinder = "none";
}


my $testGroupA = ($genomeType eq "archaea"  or $genomeType eq "auto");
my $testGroupB = ($genomeType eq "bacteria" or $genomeType eq "auto");

my $testGroupAStep1 = $testGroupA;
my $testGroupAStep2 = $testGroupA;


#----------------------------------------
# Run initial MGM prediction
#----------------------------------------
my $mgmPred = CreatePredFileName("0");                  # create a prediction filename for iteration 0
#run("$scriptPath/gmhmmp2 -M $mgmMod -s $fnseq -o $mgmPred --mgm_type $mgmType ");       # Run MGM
run("$predictor -M $mgmMod -s $fnseq -o $mgmPred --mgm_type $mgmType ");       # Run MGM

# add temporary files
push @tempFiles, ($mgmPred) unless $keepAllFiles;


# Compute probability of bacteria #bac/total; add probability to native mod file
my ($bacProb, $arcProb) = EstimateBacArc($mgmPred);

#----------------------------------------
# Main Cycle: Run X iterations 
#----------------------------------------
my $prevPred = $mgmPred;        # Previous iteration prediction: start with MGM predictions
my $prevMod = $mgmMod;          # Previous iteration model:      start with MGM model

# Run iterations without start motif model
my $iterBegin = 1;
my $iterEnd = $numIterWithoutRBS;
my $prevIter = RunIterations( { "mode" => $modeNoMotif, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );


if (defined $forceGroup) {
    my $forceMode;
    $forceMode = $modeGroupAStep1 if ($forceGroup eq "A");
    $forceMode = $modeGroupAStep2 if ($forceGroup eq "A2");
    $forceMode = $modeGroupB if ($forceGroup eq "B");
    $forceMode = $modeGroupC if ($forceGroup eq "C");
    $forceMode = $modeGroupD if ($forceGroup eq "D");
    $forceMode = $modeGroupE if ($forceGroup eq "E");

    ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);
    
    $prevIter = RunIterations( { "mode" => $forceMode, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
}
else {

    # Group A - step 1: If Group-A testing enabled, run single iteration to test for Group-A membership (step 1)
    if ($testGroupAStep1) {
        $iterBegin = $prevIter + 1;
        $iterEnd   = $iterBegin;            

        print "Testing Step-1 of Group-A membership...\n" if defined $verbose;

        $prevIter = RunIterations( { "mode" => $modeGroupAStep1, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
    }

    if ($testGroupAStep1 and IsGroupA($prevIter)) {
        print "Group-A membership: successful (by step 1).\n" if defined $verbose;

        ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);
        
        $prevIter = RunIterations( { "mode" => $modeGroupAStep1, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
    }
    else {

        # If Group A (step-1) was tested, revert iteration count and move model file
        if ($testGroupAStep1) {
            print "Group-A membership: failed.\n" if defined $verbose;
            MoveFilesFromIteration($prevIter, "groupA1");                              # revert files of failed iteration
            $prevIter -= 1;                                                 # decrement iteration counter by 1
        }


        # Group A: if Group-A testing enabled, run single iteration to test for Group-A membership
        if ($testGroupAStep2) {
            $iterBegin = $prevIter + 1;
            $iterEnd   = $iterBegin;            

            print "Testing Group-A membership...\n" if defined $verbose;

            $prevIter = RunIterations( { "mode" => $modeGroupAStep2, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
        }

        # Group-A: if membership satisfied, run remaining iterations until convergence
        if ($testGroupAStep2 and IsGroupA($prevIter)) {

            print "Group-A membership: successful.\n" if defined $verbose;

            ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);
            
            $prevIter = RunIterations( { "mode" => $modeGroupAStep2, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
        }
        # Group A: If membership not satisfied, move on to group B
        else {
            # If Group A was tested, revert iteration count and move model file
            if ($testGroupAStep2) {
                print "Group-A membership: failed.\n" if defined $verbose;
                MoveFilesFromIteration($prevIter, "groupA2");                              # revert files of failed iteration
                $prevIter -= 1;                                                 # decrement iteration counter by 1
            }

            # # If Group A was tested, revert iteration count and move model file
            # if ($testGroupAStep1) {
            #     print "Group-A membership: failed.\n" if defined $verbose;
            #     MoveFilesFromIteration($prevIter, "groupA1");                              # revert files of failed iteration
            #     $prevIter -= 1;                                                 # decrement iteration counter by 1
            # }

            # Group B: single iteration to test for Group-B membership
            if ($testGroupB) {
                $iterBegin = $prevIter + 1;
                $iterEnd = $iterBegin;          
                $prevIter = RunIterations( { "mode" => $modeGroupB, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
            }

            # Group-B: if membership satisfied, run remaining iterations until convergence
            if ($testGroupB and IsGroupB($prevIter)) {
                
                print "Group-B membership: successful.\n" if defined $verbose;

                ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);

                $prevIter = RunIterations( { "mode" => $modeGroupB, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
            }
            # Group B: If membership not satisfied, move on to group C
            else {

                MoveFilesFromIteration($prevIter, "groupB") if ($testGroupB);
                MoveFilesFromIteration($prevIter, "groupA2") if ($testGroupAStep2 and !$testGroupB);
                MoveFilesFromIteration($prevIter, "groupA1") if ($testGroupAStep1 and not !$testGroupAStep2 and !$testGroupB);
                $prevIter -= 1;

                # Group C: single iteration to test for Group-B membership
                $iterBegin = $prevIter + 1;
                $iterEnd = $iterBegin;
                $prevIter = RunIterations( { "mode" => $modeGroupC, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );

                # Group-C: if membership satisfied, run remaining iterations until convergence
                if (IsGroupC($prevIter)) {

                    print "Group-C membership: successful.\n" if defined $verbose;

                    ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);

                    $prevIter = RunIterations( { "mode" => $modeGroupC, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
                }
                # Group C: If membership not satisfied, move on to group D
                else {

                    # go back one iteration (to cancel group C)
                    MoveFilesFromIteration($prevIter, "groupC");
                    $prevIter -= 1;
                    
                    # Group D: single iteration to test for Group-B membership
                    $iterBegin = $prevIter + 1;
                    $iterEnd = $iterBegin;
                    $prevIter = RunIterations( { "mode" => $modeGroupD, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd  } );
                    
                    # Group-D: if membership satisfied, run remaining iterations until convergence
                    if (IsGroupD($prevIter)) {
                        print "Group-D membership: successful.\n" if defined $verbose;

                        ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);

                        $prevIter = RunIterations( { "mode" => $modeGroupD, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
                    }
                    # Group D: If membership not satisfied, move on to group E
                    else {
                        # go back one iteration (to cancel group D)
                        MoveFilesFromIteration($prevIter, 'groupD');
                        $prevIter -= 1;

                        ($iterBegin, $iterEnd) = GetBeginEndIterations($prevIter);
                        $prevIter = RunIterations( { "mode" => $modeGroupE, "iteration-begin" => $iterBegin, "iteration-end" => $iterEnd } );
                    } 
                }
            }
        }
    }
}


$prevPred = CreatePredFileName($prevIter);       # Prediction file: get name of previous iteration
$prevMod  = CreateModFileName($prevIter);        # Model file: get name of previous iteration



#----------------------------------------
# Clean up and get scores
#----------------------------------------


my $finalPred = $fnoutput;
my $finalMod = "GMS2.mod";
my $finalMGM = $mgmMod;  
run("cp $prevMod $finalMod");



# add bacteria and archaea probability to modfile
AddToModel($finalMod, "TO_ATYPICAL_FIRST_BACTERIA", $bacProb);
AddToModel($finalMod, "TO_ATYPICAL_SECOND_ARCHAEA", $arcProb);

if (not $fixedNativeAtypicalProb) {
    ($toNativeProb, $toMgmProb) = EstimateNativeAtypical($prevPred);
}
# add mgm and native probabilities to modfile
AddToModel($finalMod, "TO_MGM", $toMgmProb);
AddToModel($finalMod, "TO_NATIVE", $toNativeProb);
# AddToModel($finalMod, "TO_NATIVE", 15);

# check for fnn and faa output options
my $extraOutput = "";
$extraOutput .= "--AA $fnAA " if defined $fnAA;
$extraOutput .= "--NT $fnNN " if defined $fnNN;
$extraOutput .= "--gid_per_contig " if $gid;
$extraOutput .= "--defline_parse " if $ncbi;
$extraOutput .= "-e $fn_external " if $fn_external;

run("$predictor -m $finalMod -M $finalMGM -s $fn_genome -o $finalPred --format $formatOutput $extraOutput");

run ("rm -f @tempFiles");


# 
sub GetGeneticCodeFromFile
{
	my $fname = shift;

	my $genetic_code = 0;
	my %h;

	open( my $IN, $fname) or die"error on open file $fname: $!\n";
	while( my $line = <$IN> )
	{
		next if ( $line !~ /^\s*>/ );

		%h = $line =~ m/\[\s*(\S+)\s*=\s*(\S.*?)s*\]/g;

		if ( exists $h{ 'gcode' } )
		{
			if ( $genetic_code and ( $genetic_code ne $h{ 'gcode' } ))
			{
				die "error: different genetic codes are specified for the same input file: $genetic_code and $h{ 'gcode' }\n";
			}
			else
			{
				$genetic_code = $h{ 'gcode' };
			}
		}	
	}
	close $IN;

	die "error, genetic code information not found in defenition line: $_" if (!$genetic_code);

	return $genetic_code;
}

# Calculate number of iterations remaining (until we reach max number of allowed iterations)
sub NumOfIterRemaining {
    my ($prevIter, $maxIter) = @_;

    return $maxIter - $prevIter;
}

# Return true if the FGIO that don't match 16S have a localized signal located before the distance threshold
sub FGIONotMatching16SHaveSignalBeforeThresh {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    $prevPred = "itr_$prevIter.lst";
    $prevMod = "itr_$prevIter.mod";
    #my $isBacteriaProm = run("$trainer experiment promoter-is-valid-for-bacteria --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh");
    my $isBacteriaProm = run("$trainer experiment promoter-is-valid-for-bacteria --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize --min-leaderless-percent 11 --min-leaderless-count 100 --fnlabels $prevPred --fnseq $fnseq");

    chomp $isBacteriaProm;
    return $isBacteriaProm eq "yes";
}

sub RBSSignalLocalized {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $rbsIsLocalized = run("$trainer experiment rbs-is-localized --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize");
    chomp $rbsIsLocalized;

    return $rbsIsLocalized eq "yes";
}

# Return true if the FGIO have a localized motif signal located further than a distance threshold
sub FGIOHaveSignalAfterThresh {
    my ($prevIter, $distThresh, $scoreThresh, $windowSize) = @_;

    $prevMod = "itr_$prevIter.mod";
    #my $isArchaea = run("$trainer experiment promoter-is-valid-for-archaea --fnmod $prevMod --dist-thresh $distThresh");
    my $isArchaea = run("$trainer experiment promoter-is-valid-for-archaea --fnmod $prevMod --dist-thresh $distThresh --score-thresh $scoreThresh --window-size $windowSize");

    chomp $isArchaea;

    return $isArchaea eq "yes";
}

# Return true if the fraction of predicted RBS is greater than the threshold
sub PredictedRBSMatch16S {
    my ($fnpred, $seq16S, $minMatch) = @_;

    my $rbsMatchedOutput = run("$trainer experiment match-rbs-to-16s --fnlabels $fnpred --match-to $seq16S --min-match $minMatch --allow-ag-sub");

    my @matchInfo = split(' ', $rbsMatchedOutput);
    my $percentMatched = $matchInfo[1] / $matchInfo[0];

    # my $denom=system('cat $fnpred | grep -E "(native|atypical)[[:space:]]+[ACGT]+[[:space:]]+[[:digit:]]+[[:space:]]+1[[:space:]]*" | wc -l');
    # $percentMatched = $matchInfo[1] / $denom;

    print "Percent of matched RBS: $percentMatched\n" if defined $verbose;

    return $percentMatched >= $groupD_percentMatchRBS;
}

# Return true if the Promoter and RBS model consensus sequences match each other
sub PromoterAndRBSConsensusMatch {
    my ($prevIter, $minMatch) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $isGroupC = run("$trainer experiment promoter-and-rbs-match --fnmod $fnmod --match-thresh $minMatch");
    chomp $isGroupC;

    return $isGroupC eq "yes";
}

# Return true if the RBS model consensus matches the 16S tail
sub RBSConsensusAnd16SMatch {
    my ($prevIter, $minMatch) = @_;

    my $fnmod = "itr_$prevIter.mod";

    my $isMatched = run("$trainer experiment rbs-consensus-and-16s-match --fnmod $fnmod --allow-ag-sub");
    chomp $isMatched;

    return $isMatched eq "yes";
}

# Run GMS2 iterations in a particular mode
sub RunIterations {
    my %params = %{ $_[0] };

    my $mode        = $params{"mode"};
    my $iterBegin   = $params{"iteration-begin"};
    my $iterEnd     = $params{"iteration-end"};

    my $iter = $iterBegin;

    if ($iter <= 0) {
        print "Cannot run training for iteration <= 0.";
        return;
    }

    while ($iter <= $iterEnd) {

        print "Mode $mode: Entering iteration $iter...\n" if defined $verbose;

        # train on native genes
        my $nativeOnly = 0;
        if ($iter > 1) {
            $nativeOnly = 1;
        }

        my $currMod  = CreateModFileName($iter);        # model file for current iteration
        my $currPred = CreatePredFileName($iter);       # prediction file for current iteration
        my $prevPred = CreatePredFileName($iter-1);

        # Training step: use prediction of previous iteration
        my $trainingCommand = GetTrainingCommand($iter, $mode);          # construct training command
        run("$trainingCommand");                                         # run training command

        # add bacteria and archaea probability to model file
        AddToModel($currMod, "TO_ATYPICAL_FIRST_BACTERIA", $bacProb);
        AddToModel($currMod, "TO_ATYPICAL_SECOND_ARCHAEA", $arcProb);

        if (not $fixedNativeAtypicalProb and $iter > 1) {
            ($toNativeProb, $toMgmProb) = EstimateNativeAtypical($prevPred);
        }
        # add mgm and native probabilities to modfile
        AddToModel($currMod, "TO_MGM", $toMgmProb);
        AddToModel($currMod, "TO_NATIVE", $toNativeProb);

        # Prediction step: using current model file
        my $errCode = run("$predictor -m $currMod -M $mgmMod -s $fnseq -o $currPred --format train");

        # Check for convergence
        my $similarity = run("$comparePrediction -n -a $prevPred -b $currPred -G");

        print "Iteration : $similarity\n" if defined $verbose;


        # add temporary files
        push @tempFiles, ($currMod, $currPred) unless $keepAllFiles;
        

        if ( ($similarity > 99 && $iter > 2) ) {
            print "Converged at iteration $iter\n" if defined $verbose;
            return $iter;
        }

        # set previous prediction (before exiting loop)
        $prevPred = $currPred;
        $prevMod = $currMod;


        $iter++;
    }

    return $iter-1;
}

# Run a system command and log it
sub run {
    my $command = shift;
    open(FILE, ">>log");
    print FILE $command . "\n";
    my $value = `$command`;
    chomp($value);
    return $value;
}

# Estimate bacteria and archaea probabilities based on the counts in the prediction file
sub EstimateBacArc {
    my $fname = shift;
    my $counts_all = 0;
    my $counts_bac = 0;

    my $min_gene_length_bac_arc = 600;

    open( my $IN , $fname ) or die "error on open file $fname: $!\n";
    while( my $line = <$IN>)
    {
        next if ( $line =~ /^\s*$/ );
        next if ( $line =~ /^#/ );
        next if ( $line =~ /SequenceID:/ );

        if ( $line =~ /^\s*\d+\s+[+-]\s+\S+\s+\S+\s+(\d+)\s+bac\s*/ )
        {
            if ( $1 >= $min_gene_length_bac_arc )
            {
                ++$counts_bac;
                ++$counts_all;
            }
        }
        elsif ( $line =~ /^\s*\d+\s+[+-]\s+\S+\s+\S+\s+(\d+)\s+arc\s*/ )
        {
            if ( $1 >= $min_gene_length_bac_arc )
            {
                ++$counts_all;
            }
        }
        else {die;}
    }
    close $IN;

    if (!$counts_all) {print "error, unexpected format foundin file: $fname"; exit 1; }

    if (defined $verbose) {
        my $counts_arc = $counts_all - $counts_bac;
        my $bacProb = $counts_bac / $counts_all;
        my $arcProb = $counts_arc / $counts_all;
        print "NumBac = $counts_bac\n";
        print "NumArc = $counts_arc\n";
        print "Bacteria Probability: $bacProb\n";
        print "Archaea Probability: $arcProb\n";
    }

    return ( sprintf( "%.5f", $counts_bac/$counts_all ), sprintf( "%.5f", ($counts_all - $counts_bac)/$counts_all ) );
}

# Esitmate native and atypical probabilities based on the counts in the prediction file
sub EstimateNativeAtypical {
    my $fname = shift;
    my $counts_all = 0;
    my $counts_native = 0;

    my $min_gene_length_native_atypical = 600;

    open( my $IN , $fname ) or die "error on open file $fname: $!\n";
    while( my $line = <$IN>)
    {
        next if ( $line =~ /^\s*$/ );
        next if ( $line =~ /^#/ );
        next if ( $line =~ /SequenceID:/ );

        my $currLength = -1;
        if ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s*/ ) {
            $currLength = $1;
        }

        # if ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+native\s*/ )
        if ( $line =~ /\s+native\s*$/ ) 
        {
            if ( $currLength >= $min_gene_length_native_atypical )
            {
                ++$counts_native;
                ++$counts_all;
            }
        }
        #elsif ( $line =~ /^\s*\d+\s+[+-]\s+\d+\s+\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+atypical\s*/ )
        elsif ( $line =~ /\s+atypical\s*$/ ) 
        {
            if ( $currLength >= $min_gene_length_native_atypical )
            {
                ++$counts_all;
            }
        }
        #else {die;}
    }
    close $IN;

    if (!$counts_all) {print "error, unexpected format foundin file: $fname"; exit 1; }

    if (defined $verbose) {
        my $counts_atypical = $counts_all - $counts_native;
        my $nativeProb = $counts_native / $counts_all;
        my $atypicalProb = $counts_atypical / $counts_all;
        print "NumNative = $counts_native\n";
        print "NumAtypical = $counts_atypical\n";
        print "Native Probability: $nativeProb\n";
        print "Atypical Probability: $atypicalProb\n";
    }

    my $probNative = $counts_native/$counts_all;
    my $probAtypical = ($counts_all - $counts_native)/$counts_all;
    if ($probAtypical < $minAtypicalProb) {
        $probAtypical = $minAtypicalProb;
        $probNative = 1 - $probAtypical;
    }

    return ( sprintf( "%.5f", $probNative ), sprintf( "%.5f", $probAtypical) );
}

# Converts a multifasta file to single fasta by concatenating all sequences into one
sub MultiToSingleFASTA {
    # input and output filenames
    my ($fnin, $fnout) = @_;
    
    run ("echo '>anydef' > $fnout");
    run ("grep -v '>' $fnin | tr '[:lower:]' '[:upper:]' >> $fnout");
    return;
}

# Add label/value pair to a model file
sub AddToModel {
    my ( $fname, $label, $value ) = @_;
    open (my $fout, ">>", $fname) or die "Error: Could not open file $fname\n";

    print $fout "\$" . $label . " $value\n";
    close $fout;
}

# Create name for model file based on iteration number
sub CreateModFileName {
    my $iter = $_[0];
    return "itr_$iter.mod";
}

# Create name for prediction file based on iteration number
sub CreatePredFileName {
    my $iter = $_[0];
    return "itr_$iter.lst";
}

# Create name for Motif Finder Output file based on iteration number
sub CreateMFinderResultFileName {
    my $iter = $_[0];
    return "itr_$iter.mfinderresult";
}

# Returns true if the genome type is valid: archaea, bacteria, auto
sub isValidGenomeType {
    my $gt = $_[0];
    return ($gt eq "archaea" or $gt eq "bacteria" or $gt eq "auto");
}





##############################
#                            #
#   Group Membership Tests   #
#                            #
##############################




sub IsGroupA {
    my $iter = $_[0];

    my $testResult = FGIOHaveSignalAfterThresh($iter, $groupA_spacerDistThresh, $groupA_spacerScoreThresh, $groupA_spacerWindowSize);

    return $testResult;
}

sub IsGroupB {
    my $iter = $_[0];

    my $test1 = FGIONotMatching16SHaveSignalBeforeThresh($iter, $groupB_spacerDistThresh, $groupB_spacerScoreThresh, $groupB_spacerWindowSize);
    my $test2 = PromoterAndRBSConsensusMatch($iter, $groupC_minMatchPromoterRBS);

    $testGroupB_PromoterMatchedRBS = $test2;

    if ($test1 and not $test2) {
        return 1;
    }
    else {
        return undef;
    }
}

sub IsGroupC {
    my $iter = $_[0];

    my $test = RBSConsensusAnd16SMatch($iter, $groupC_minMatchRBS16S);

    if (!$test) {

        if (RBSSignalLocalized($iter, 14, 0.15, 1)) {
            return 1;
        }
        return 0;
        # return 1;

        if ($testGroupB) {
            
            if ($testGroupB_PromoterMatchedRBS) {
                return 1;
            }
            else {
                return undef;
            }
        }
        # if group B wasn't tested for
        else {
            return 1;
        }
    }
    else {
        return undef;
    }
}

sub IsGroupD {
    my $iter = $_[0];

    my $fnpred = CreatePredFileName($iter);

    my $test = PredictedRBSMatch16S($fnpred, $groupD_tail16S, $groupD_minMatchRBS16S);

    return $test;

}


##############################
#                            #
#   Auxilliary Functions     #
#                            #
##############################


sub GetTrainingCommand {
    my ($currIter, $mode) = @_;

    if ($currIter == 0) {
        print "Cannot construct training model at iteration 0";
        return undef;
    }

    my $prevIter = $currIter - 1;

    my $currMod  = CreateModFileName($currIter);        # model file for current iteration
    my $prevPred = CreatePredFileName($prevIter);       # prediction file of previous iteration


    # train on native genes
    my $nativeOnly = 0;
    if ($currIter > 1) {
        $nativeOnly = 1;
    }
    

    # Training step: use prediction of previous iteration
    my $trainingCommand = "$trainer gms2-training -s $fnseq -l $prevPred -m $currMod --order-coding $orderCod --order-noncoding $orderNon --only-train-on-native $nativeOnly --genetic-code $geneticCode --order-start-context $scOrder --fgio-dist-thr $fgioDistThresh";

    #$trainingCommand .= " --len-start-context 6 --margin-start-context -3 ";

    if ($mode eq $modeNoMotif) {
        $trainingCommand .= " --run-motif-search false";
        $trainingCommand .= " --genome-group D";        # FIXME: remove requirement from training
    }
    elsif ($mode eq $modeGroupAStep1) {
        $trainingCommand .= " --genome-group A --ga-upstr-len-rbs $groupA_rbsUpstreamLength --align $alignmentInMFinder --ga-width-rbs $groupA_widthRBS --ga-upstr-len-prom $groupA_promoterUpstreamLength --ga-width-prom $groupA_widthPromoter";            
    }
    elsif ($mode eq $modeGroupAStep2) {
        $trainingCommand .= " --genome-group A2 --ga-upstr-len-rbs $groupB_rbsUpstreamLength --align $alignmentInMFinder --ga-width-rbs $groupB_widthRBS --ga-upstr-len-prom $groupA_promoterUpstreamLength --ga-width-prom $groupA_widthPromoter --ga-extended-sd $groupB_tail16S";            
    }
    elsif ($mode eq $modeGroupB) {
        $trainingCommand .= " --genome-group B --gb-upstr-len-rbs $groupB_rbsUpstreamLength --align $alignmentInMFinder --gb-width-rbs $groupB_widthRBS --gb-upstr-len-prom $groupB_promoterUpstreamLength --gb-width-prom $groupB_widthPromoter --gb-extended-sd $groupB_tail16S";
    }
    elsif ($mode eq $modeGroupC) {
        $trainingCommand .= " --genome-group C --gc-upstr-len-rbs $groupC_rbsUpstreamLength --align $alignmentInMFinder --gc-width-rbs $groupC_widthRBS"; #  --gc-upstr-reg-3-prime 3";  # 
        # $trainingCommand .= " --genome-group C2 --align $alignmentInMFinder ";  # --gc-upstr-reg-3-prime 3
    }
    elsif ($mode eq $modeGroupD) {
        $trainingCommand .= " --genome-group D --gd-upstr-len-rbs $groupD_rbsUpstreamLength --align $alignmentInMFinder --gd-width-rbs $groupD_widthRBS";
    }
    elsif ($mode eq $modeGroupE) {
        $trainingCommand .= " --genome-group E --ge-upstr-len-rbs $groupE_rbsUpstreamLength --align $alignmentInMFinder --ge-width-rbs $groupE_widthRBS --ge-len-upstr-sig $groupE_upstreamSignatureLength --ge-order-upstr-sig $groupE_upstreamSignatureOrder --ge-extended-sd $groupE_tail16S";
    }
    else {
        die "Mode invalid: should not reach this point";
    }

    return $trainingCommand;
}

sub MoveFilesFromIteration {
    my $iter = $_[0];
    my $name = $_[1];

    my $fnmod  = CreateModFileName  ($iter);
    my $fnpred = CreatePredFileName ($iter);

    print "Move files from iteration $iter to name '$name'\n" if defined $verbose;

    run("mv $fnmod  $name.mod");
    run("mv $fnpred $name.lst");
}


sub GetBeginEndIterations {
    my ($prevIter) = @_;

    my $numIterRemain = NumOfIterRemaining($prevIter, $MAX_ITER);

    $iterBegin = $prevIter + 1;
    $iterEnd = $iterBegin + $numIterRemain - 1;

    return ($iterBegin, $iterEnd);
}

# FIXME: fnn, faa, check gcode 11,4, keep-all-files




# Usage function: print usage message and exit script
sub Usage {
    my $name = $_[0];
    print "Usage: $name --seq SEQ --genome-type TYPE";
    print
"
Basic Options: 
--seq                                   File containing genome sequence in FASTA format
--genome-type                           Type of genome: archaea, bacteria, auto 
--gcode                                 The genetic code number (default: $D_GENETIC_CODE. Choices: 11 and 4)
--output                                Name of output file (default: $D_FNOUTPUT)
--format                                Format of output file (default: $D_FORMAT_OUTPUT)
--ext                                   Name of file with external information in GFF format (PLUS mode of GMS2)
--fnn                                   Name of output file that will hold nucleotide sequences of predicted genes
--faa                                   Name of output file that will hold protein sequences of predicted genes
--gid                                   Change gene ID format
--advanced-options                      Show the advanced options

Version: 1.07
";

    if (defined $showAdvancedOptions) {
        print 
"
Advanced Options:
# Iteration control
--max-iter                              Number of max iterations in the main cycle (default: $D_MAX_ITER)
--conv-thresh                           The convergence threshold (in range [0,1]) (default: $D_CONV_THRESH)

# Misc
fixed-native-atypical-prob              Fix the native and atypical prior probabilities
train-noncoding-on-full-genome          Train the non-coding model on the full genome
min-atypical-prob                       Set minimum prior probability for atypical genes
run-mfinder-without-spacer              Disable the \"location distribution\" in the motif finders
mgm-type                                Type of genome model to use for MGM predictions
                                        Option: bac, arc, auto. Default: (default: $D_MGMTYPE)
keep-all-files                          Keep all intermediary files 
fgio-dist-thresh                        Distance threshold for FGIO identification

# Group-A
group-a-width-promoter                  Width of the promoter motif model (default: $D_PROM_WIDTH_A)
group-a-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
group-a-promoter-upstream-length        Upstream length for promoter training (default: $D_PROM_UPSTR_LEN_A)
group-a-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
group-a-spacer-score-thresh             Minimum peak threshold for the spacer distribution (default: $D_SPACER_SCORE_THRESH_A)
group-a-spacer-dist-thresh              Minimum distance threshold for the spacer distribution (default: $D_SPACER_DIST_THRESH)
group-a-spacer-window-size              Window size for calculating the \"peak value\" to compare
                                        against the score threshold. (default: $D_SPACER_WINDOW_SIZE)
# Group-B
group-b-width-promoter                  Width of the promoter motif model (default: $D_PROM_WIDTH_B)
group-b-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
group-b-promoter-upstream-length        Upstream length for promoter training (default: $D_PROM_UPSTR_LEN_B)
group-b-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
group-b-spacer-score-thresh             Minimum peak threshold for the spacer distribution (default: $D_SPACER_SCORE_THRESH_B)
group-b-spacer-window-size              Window size for calculating the \"peak value\" to compare
                                        against the score threshold (default: $D_SPACER_WINDOW_SIZE)
group-b-tail-16s                        The 16S rRNA tail used for selecting training sequences for
                                        the promoter model (default: $D_16S)
group-b-min-match-to-tail               Minimum number of consecutive nucleotide matches to the 16S (default: $D_MIN_MATCH_16S)

# Group-C
group-c-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
group-c-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
group-c-min-match-promoter-rbs          Minimum number of consecutive nucleotide matches between the 
                                        promoter and RBS (default: $D_MIN_MATCH_RBS_PROM)

# Group-D
group-d-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
group-d-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
group-d-percent-match-rbs               Minimum percentage of predicted RBS sequences that match to 16S (default: $D_MIN_FRAC_RBS_16S_MATCH)

# Group-E
group-e-width-rbs                       Width of the rbs motif model (default: $D_RBS_WIDTH)
group-e-rbs-upstream-length             Upstream length for rbs training (default: $D_RBS_UPSTR_LEN)
group-e-upstream-signature-length       Length of the upstream-signature Nonuniform Markov model (default: $D_UPSTR_SIG_LENGTH)
group-e-upstream-signature-order        Order of the upstream-signature Nonuniform Markov model (default: $D_UPSTR_SIG_ORDER)
group-e-tail-16s                        The 16S rRNA tail used for selecting training sequences for
                                        the RBS model (default: $D_16S)
";
    }

    exit;
}
