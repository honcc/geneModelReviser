#!/usr/bin/perl/ -w
$|++;
use strict;
use Math::Random; #---For generating normal distribution random numbers
use File::Path; #---for removing tmp GNUPLOT dat files
use Time::HiRes qw (sleep);
######################################################################################################################################################
#
#	Description
#		This is a perl script to revise the gene model (provided as a gff) using the NGS data (generated from geneModelBulider.pl/cufflinks as a gtf). It will 
#	first find the overlapping between the existing gene model (refGff) and the NGS data (NGSGtf, normally generated from geneModelBulider.pl/cufflinks), then
#	it will revise the gene models using the refGFF as the guide, and combines the splcing junction information.
#
#	Input
#		--refGff=				file path; the exising gene model Gff, i.e. reference;
#		--NGSGtf=				file path; the Gff for NGS data;
#		--junctBEB=				file path; the bed file that contains all junctions; To pickup the junctions that is absent in major and minor isoforms;
#		--pileupPath=			file path; the pileup file generated from pileup counter, used to determine the ration of splice and unspliced read flanking the splicing jucntions;
#		--intronBoundCovPath=	"no" or the path of a file that contains the intron bound coverage of the junctions. If "no", the coverage will be calculated from the pileu file;
#		--fastaSeq=				file path; a file that contains the reference sequence fasta;
#		--junctionBedType=		"HMMSplicer" or "tophat"; as the bed is a bit different between the two, mainly on supporting read numnber; "HMMSplicer" fits the output from HMMSplicerBEDToSAMParser, "tophat" fits for junctions.bed from tophat;
#		--strndSpecificNGS=		'yes' or 'no'; [no]; the NGS data is strand specific or not;
#		--outDir=				dir; output directory;
#
#	Output
#
#	Usage
#		perl geneModelReviser_v0.6.pl --refGff=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/all.revised.withRepElmntgff --NGSGtf=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/geneModel/finalGTF/allIsoform.gtf --outDir=./s2.s3.s4.s6.s7.s8.pooled/ --junctBEB=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/HMMSplicerBEDToSAMParser/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed --pileupPath=no --intronBoundCovPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviser/v0.5/s2.s3.s4.s6.s7.s8.pooled/refIntronBoundCov.30.txt --fastaSeq=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa
#
#	Version History
#		v0.2
#			-readRefGff subroutine rewritten;
#
#		v0.3
#			-completely rewritten: originally designed for output of transfragDiscoverer and old-style reference GFF (changed to new style in v0.2). Now it is rewritten for output of gene model reviser and new-style of reference gff.
#
#		v0.4
#			-corrected for the output exon ranges instead of gene range in the NGS unpair file;
#
#		v0.5
#			-will ouput info of novel junctions;
#
#		v0.6
#			-added the functionality to scan for proximity of N region and cntg ends, quick and dirty code inherited from motifProximityScanner_v0.1.pl
#
#		v0.7
#		- added junctionBedType option;
#		- modified checkOverlapAndProximity subroutnine so that it recognize the non-strand specific GTF from transfragDiscoverer;
#		- will summerize gene based info;
#
#		v0.8
#		- added --strndSpecificNGS= option;
#
#		v0.9
#		- relaxed 'passed' criteria, coverage down to 70pct and cancelled the extended transfrag comment
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#

#1----------Read the parameters----------#
use vars qw ($refGff $NGSGtf $junctBEB $pileupPath $intronBoundCovPath $fastaSeq $junctionBedType $strndSpecificNGS $outDir);
($refGff, $NGSGtf, $junctBEB, $pileupPath, $intronBoundCovPath, $fastaSeq, $junctionBedType, $strndSpecificNGS, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#2----------Read the Gff----------#
my ($refRngSSHsh_ref, $refCtgryByGeneHsh_ref, $refTrnscptStrndHsh_ref, $refRngXSHsh_ref, $refExonRngByGeneHsh_ref, $refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref) = readRefGff();
my ($NGSRngSSHsh_ref, $NGSStrndHsh_ref, $NGSTrnscptRngXSHsh_ref, $NGSGeneByTrnscptHSH_ref, $NGSExonRngByTrnsptHsh_ref, $NGSIntronRngByGeneHsh_ref, $NGSJunctStrByIntronIDHsh_ref, $NGSIntronRngXS_ref, $NGSIntronStrndHsh_ref, $NGSCntgByTrnscptHsh_ref, $uniqueNGSJunctStrByIntronIDHsh_ref) = readNGSGtf();

my $refSeqHsh_ref = readMultiFasta($fastaSeq);
my ($ATRichRngXSHsh_ref, $NRegionRngXSHsh_ref, $cntgEndRngXSHsh_ref) = getMotifRng($refSeqHsh_ref);
my $ATRichStrndHsh_ref = assignDummyStrand($ATRichRngXSHsh_ref);
my $NRegionStrndHsh_ref = assignDummyStrand($NRegionRngXSHsh_ref);
my $cntgEndStrndHsh_ref = assignDummyStrand($cntgEndRngXSHsh_ref);

my ($dummy, $NRegionNGSHitHsh_ref, $NRegionNGSPrxmtyHsh_ref, $cntgEndNGSHitHsh_ref, $cntgEndNGSPrxmtyHsh_ref, $NRegionRefHitHsh_ref, $NRegionRefPrxmtyHsh_ref, $cntgEndRefHitHsh_ref, $cntgEndRefPrxmtyHsh_ref);
#---------check NGS with N region with NGS---------#
($dummy, $dummy, $NRegionNGSHitHsh_ref, $dummy, $dummy, $dummy, $NRegionNGSPrxmtyHsh_ref, $dummy) = checkOverlapAndProximity($NGSTrnscptRngXSHsh_ref, $NRegionRngXSHsh_ref, $NGSStrndHsh_ref, $NRegionStrndHsh_ref, "Finding NGS overlapping and proximty with N regions ", "yes");

#---------check NGS with cntgEnd with NGS---------#
($dummy, $dummy, $cntgEndNGSHitHsh_ref, $dummy, $dummy, $dummy, $cntgEndNGSPrxmtyHsh_ref, $dummy) = checkOverlapAndProximity($NGSTrnscptRngXSHsh_ref, $cntgEndRngXSHsh_ref, $NGSStrndHsh_ref, $cntgEndStrndHsh_ref, "Finding NGS overlapping and proximty with cntg Ends ", "yes");

#---------check NGS with N region with NGS---------#
($dummy, $dummy, $NRegionRefHitHsh_ref, $dummy, $dummy, $dummy, $NRegionRefPrxmtyHsh_ref, $dummy) = checkOverlapAndProximity($refRngXSHsh_ref, $NRegionRngXSHsh_ref, $refTrnscptStrndHsh_ref, $NRegionStrndHsh_ref, "Finding Ref overlapping and proximty with N regions ", "yes");

#---------check NGS with cntgEnd with NGS---------#
($dummy, $dummy, $cntgEndRefHitHsh_ref, $dummy, $dummy, $dummy, $cntgEndRefPrxmtyHsh_ref, $dummy) = checkOverlapAndProximity($refRngXSHsh_ref, $cntgEndRngXSHsh_ref, $refTrnscptStrndHsh_ref, $cntgEndStrndHsh_ref, "Finding Ref overlapping and proximty with cntg Ends ", "yes");

my ($NGSPrmntTrnscptRngXSHsh_ref, $NGSPrmntStrndHsh_ref) = filterProminentIsoform($NGSTrnscptRngXSHsh_ref, $NGSStrndHsh_ref);

#----calculate the bound coverage of ref intron
my $refJunctIntronBoundCovHsh_ref= scanJunctionInnerOutterBoundXSCoverage($refIntronRngXS_ref, 30);

#3---------find overlap---------#
my ($SSHitByRefTrnscptHsh_ref, $SSHitByNGSTrnscptHsh_ref, $XSHitByRefTrnscptHsh_ref, $XSHitByNGSTrnscptHsh_ref, $SSPrxmtyByRefTrnscptHsh_ref, $SSPrxmtyByNGSTrnscptHsh_ref, $XSPrxmtyByRefTrnscptHsh_ref, $XSPrxmtyByNGSTrnscptHsh_ref) = checkOverlapAndProximity($refRngXSHsh_ref, $NGSPrmntTrnscptRngXSHsh_ref, $refTrnscptStrndHsh_ref, $NGSPrmntStrndHsh_ref, "Finding overlapping and proximal ref and prominent NGS transcripts", "yes");

#4------find the hit and unpair----------#
my ($refXSUnpairLenHsh_ref, $NGSXSUnpairLenHsh_ref, $refXSHitCountHsh_ref, $NGSXSHitCountHsh_ref, $refBasedNGSHitInfoHsh_ref) = findHitAndUnpair($SSHitByRefTrnscptHsh_ref, $SSHitByNGSTrnscptHsh_ref, $XSHitByRefTrnscptHsh_ref, $XSHitByNGSTrnscptHsh_ref, $refTrnscptStrndHsh_ref, $NGSPrmntStrndHsh_ref, $refRngXSHsh_ref, $NGSPrmntTrnscptRngXSHsh_ref, $NGSExonRngByTrnsptHsh_ref, $strndSpecificNGS);

#5---summarize the hit and unpair
summarizeHitAndUnpair($refXSUnpairLenHsh_ref, $NGSXSUnpairLenHsh_ref, $refXSHitCountHsh_ref, $NGSXSHitCountHsh_ref, $refCtgryByGeneHsh_ref, $SSPrxmtyByRefTrnscptHsh_ref, $SSPrxmtyByNGSTrnscptHsh_ref, $XSPrxmtyByRefTrnscptHsh_ref, $XSPrxmtyByNGSTrnscptHsh_ref, $NGSExonRngByTrnsptHsh_ref, $NGSPrmntTrnscptRngXSHsh_ref, $NGSCntgByTrnscptHsh_ref, $NGSStrndHsh_ref, $NRegionNGSHitHsh_ref, $NRegionNGSPrxmtyHsh_ref, $cntgEndNGSHitHsh_ref, $cntgEndNGSPrxmtyHsh_ref);

#6------check ref and NGS overlapping junction
my ($SSNGSIntronHitByRefHsh_ref, $SSRefIntronHitByNGSHsh_ref, $XSNGSIntronHitByRefHsh_ref, $XSRefIntronHitByNGSHsh_ref) = checkOverlapAndProximity($refIntronRngXS_ref, $NGSIntronRngXS_ref, $refStrndByIntronHsh_ref, $NGSIntronStrndHsh_ref, "Finding overlapping ref and NGS introns", "no");

#7----------Read the BED----------#
my ($BEDStrndByIntronHsh_ref, $BEDJunctRngHsh_ref, $BEDJunctCntgHsh_ref, $BEDJunctScoreHsh_ref, $BEDJunctReadNumHsh_ref, $BEDIntronRngXS_ref) = readJunctBED($junctBEB);

#8------check ref and BED overlapping junction
my ($SSBEDIntronHitByRefHsh_ref, $SSBEDIntronHitByBEDHsh_ref, $XSBEDIntronHitByRefHsh_ref, $XSBEDIntronHitByBEDHsh_ref) = checkOverlapAndProximity($refIntronRngXS_ref, $BEDIntronRngXS_ref, $refStrndByIntronHsh_ref, $BEDStrndByIntronHsh_ref, "Finding overlapping ref and BED introns", "no");

#9-----summarize Junction Overlapping Scenarios
my $allScenarioByRefIntronIDHsh_ref = summarizeRefJunctOverlapScenarios($SSNGSIntronHitByRefHsh_ref, $refJunctStrByIntronIDHsh_ref, $NGSJunctStrByIntronIDHsh_ref, $SSBEDIntronHitByRefHsh_ref, $refJunctIntronBoundCovHsh_ref);

#10---Checking overlapping and proximity between NGS junctions and refModels
my ($SSNGSIntronHitByRefTrnscptHsh_ref, $SSRefTrnscptHitByNGSIntronHsh_ref, $XSNGSIntronHitByRefTrnscptHsh_ref, $XSRefTrnscptHitByNGSIntronHsh_ref, $SSNGSIntronPrxmtyByRefTrnscptHsh_ref, $SSRefTrnscptPrxmtyByNGSIntronHsh_ref, $XSNGSIntronPrxmtyByRefTrnscptHsh_ref, $XSRefTrnscptPrxmtyByNGSIntronHsh_ref) = checkOverlapAndProximity($refRngXSHsh_ref, $NGSIntronRngXS_ref, $refTrnscptStrndHsh_ref, $NGSIntronStrndHsh_ref, "Checking overlapping and proximity between NGS junctions and refModels", "yes");

my $NGSJRefPrxmtyLimit = 0; #----minimal distance between and NGS intron to refRng to be regarded as in proximity, use "zero" to disable

my $NGSJunctOnRefCommentHsh_ref = summarizeNGSJunctOverlapScenarios($SSNGSIntronHitByRefTrnscptHsh_ref, $SSNGSIntronPrxmtyByRefTrnscptHsh_ref, $SSRefTrnscptPrxmtyByNGSIntronHsh_ref, $SSRefTrnscptHitByNGSIntronHsh_ref, $SSRefIntronHitByNGSHsh_ref,  $refJunctStrByIntronIDHsh_ref, $uniqueNGSJunctStrByIntronIDHsh_ref, $allScenarioByRefIntronIDHsh_ref, $refCtgryByGeneHsh_ref, $refTrnscptStrndHsh_ref, $NGSIntronStrndHsh_ref, $XSRefTrnscptHitByNGSIntronHsh_ref, $NGSJRefPrxmtyLimit);

my $cutoffCovPct = 60;
$refBasedNGSHitInfoHsh_ref = summerizeRefTranscriptBasedNGSHitInfo($allScenarioByRefIntronIDHsh_ref, $refBasedNGSHitInfoHsh_ref, $refIntronRngByGeneHsh_ref, $NRegionRefHitHsh_ref, $NRegionRefPrxmtyHsh_ref, $cntgEndRefHitHsh_ref, $cntgEndRefPrxmtyHsh_ref, $NGSJunctOnRefCommentHsh_ref, $cutoffCovPct);

printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$pileupPath= "no";
	$intronBoundCovPath = "no";
	$junctionBedType = "tophat";
	$strndSpecificNGS = 'no';
	
	foreach my $param (@ARGV) {

		if ($param =~ m/--refGff=/) {$refGff = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--NGSGtf=/) {$NGSGtf = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--junctBEB=/) {$junctBEB = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--pileupPath=/) {$pileupPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--intronBoundCovPath=/) {$intronBoundCovPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--fastaSeq=/) {$fastaSeq = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--junctionBedType=/) {$junctionBedType = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--strndSpecificNGS=/) {$strndSpecificNGS = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);}
	}

	system ("mkdir -p -m 777 $outDir");
	system ("mkdir -p -m 777 $outDir/unpariedTransfrag/");
	system ("mkdir -p -m 777 $outDir/RefJunctions/");
	system ("mkdir -p -m 777 $outDir/NGSJunctions/");
	system ("mkdir -p -m 777 $outDir/hitlog/");
	system ("mkdir -p -m 777 $outDir/GFF/");

	return ($refGff, $NGSGtf, $junctBEB, $pileupPath, $intronBoundCovPath, $fastaSeq, $junctionBedType, $strndSpecificNGS, $outDir);
}
########################################################################## readRefGff
sub readRefGff {

	#---the whole subrountine was inherited from pileupCounter_v0.8 so it looks a bit redundant. Will come back later to clean it up.

	#---variables to retun
	my (%strndByGeneHsh, %cntgByGeneHsh, %exonRngByGeneHsh, %exonNumByCntgHsh, %geneExonLenHsh, %geneCDSLenHsh, %ctgryReadCountHsh, %CDSRngByGeneHsh, %geneByCtgryHsh, %ctgryByGeneHsh);
	my (%geneByRNAHsh, %CDSCountHsh, %exonCountHsh, %geneExonLocationHsh);
	my (%SSRngByCntgByWholeGeneHsh, %XStrndRngByCntgByWholeGeneHsh);

	#---read the gff
	open (INFILE, $refGff) || die "Cannot open $refGff";
	print "Reading $refGff for categorizing the features.\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
			my @theLineSplt = split (/\t/, $theLine);
			my $seq = $theLineSplt[0];

			my $geneCategory = $theLineSplt[2];
			my $featureStart = $theLineSplt[3];
			my $featureEnd = $theLineSplt[4];
			my $geneStrnd = $theLineSplt[6];
			my $dscrptns = $theLineSplt[8];
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}
	
			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				$strndByGeneHsh{$geneID} = $geneStrnd;
				$cntgByGeneHsh{$geneID} = $seq;
				${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$geneID}}{"start"} = $featureStart;
				${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$geneID}}{"end"} = $featureEnd;

				${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$geneID}}{"start"} = $featureStart;
				${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$geneID}}{"end"} = $featureEnd;

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				# The CDS is ignored at the moment, until it reaches the point that we are looking at UTRs
				#
				#my $mRNAID = $parent;
				#my $geneID = $geneByRNAHsh{$mRNAID};
				#$CDSCountHsh{$geneID}++;
				#my $CDSCount = $CDSCountHsh{$geneID};
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"} = $featureStart;
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} = $featureEnd;
			 	#$geneCDSLenHsh{$geneID} = 0 if $CDSCount == 1; #---define the length hashfor the 1st time
			 	#$geneCDSLenHsh{$geneID} += ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} - ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"};
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $exonID = $unqID;
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh{$RNAID};
				my $locationTag = $seq.":".$featureStart.":".$featureEnd;

				${$geneExonLocationHsh{$locationTag}}{$geneID}++;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} = $featureEnd;
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh{$RNAID} = $geneID;
				$ctgryByGeneHsh{$geneID} = $geneCategory;
				$geneByCtgryHsh{$geneCategory} = $geneID;
				
				if (not(exists $ctgryReadCountHsh{$geneCategory})) {#---initialize the $geneCategory for all category 
					${$ctgryReadCountHsh{$geneCategory}}{"s"} = 0;
					${$ctgryReadCountHsh{$geneCategory}}{"a"} = 0;
				}
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;

	my ($refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref) = getIntronFromExonRng(\%exonRngByGeneHsh, \%cntgByGeneHsh, \%strndByGeneHsh, "Getting the intron boundaries on refGff", "ref");
	
	return (\%SSRngByCntgByWholeGeneHsh, \%ctgryByGeneHsh, \%strndByGeneHsh, \%XStrndRngByCntgByWholeGeneHsh, \%exonRngByGeneHsh, $refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref);
}
########################################################################## readNGSGtf
sub readNGSGtf {

#DS571327	geneModelBuilder	transcript	24695	25646	1000	+	.	gene_id "JT_1004"; transcript_id "JT_1004.0";
#DS571327	geneModelBuilder	exon	24695	25213	1000	+	.	gene_id "JT_1004"; transcript_id "JT_1004.0"; exon_number "1";
#DS571327	geneModelBuilder	exon	25277	25646	1000	+	.	gene_id "JT_1004"; transcript_id "JT_1004.0"; exon_number "2";
#DS571287	geneModelBuilder	transcript	12395	14104	1000	+	.	gene_id "JT_1005"; transcript_id "JT_1005.0";
#DS571287	geneModelBuilder	exon	12395	12924	1000	+	.	gene_id "JT_1005"; transcript_id "JT_1005.0"; exon_number "1";
#DS571287	geneModelBuilder	exon	13100	13243	1000	+	.	gene_id "JT_1005"; transcript_id "JT_1005.0"; exon_number "2";
#DS571287	geneModelBuilder	exon	13397	14104	1000	+	.	gene_id "JT_1005"; transcript_id "JT_1005.0"; exon_number "3";


	#---the whole subrountine was inherited from pileupCounter_v0.8 so it looks a bit redundant. Will come back later to clean it up.

	#---variables to retun
	my (%strndByTrnscptHsh, %cntgByTrnscptHsh, %exonRngByGeneHsh, %exonNumByCntgHsh, %geneExonLenHsh, %geneCDSLenHsh, %ctgryReadCountHsh, %CDSRngByGeneHsh, %geneByCtgryHsh, %ctgryByGeneHsh, %geneByTrnscptHsh);
	my (%geneByRNAHsh, %CDSCountHsh, %exonCountHsh, %geneExonLocationHsh);
	my (%SSRngByCntgByWholeGeneHsh, %XStrndRngByCntgByWholeGeneHsh);

	#---read the gff
	open (INFILE, $NGSGtf) || die "Cannot open $NGSGtf";
	print "Reading $NGSGtf for categorizing the features.\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		
		my @theLineSplt = split (/\t/, $theLine);
		my $seq = $theLineSplt[0];

		my $geneCategory = $theLineSplt[2];
		my $featureStart = $theLineSplt[3];
		my $featureEnd = $theLineSplt[4];
		my $geneStrnd = $theLineSplt[6];
		my $dscrptns = $theLineSplt[8];
		my @dscrptnsSplt = split /;/, $dscrptns;
		my ($geneID, $trnscptID, $exonNum);
		
		foreach my $theDscptn (@dscrptnsSplt) {
			my @theDscptnSplt = split /\"/, $theDscptn;
			if ($theDscptn =~ m/gene_id/) {$geneID = $theDscptnSplt[1];}
			if ($theDscptn =~ m/transcript_id/) {$trnscptID = $theDscptnSplt[1];}
			if ($theDscptn =~ m/exon_number/) {$exonNum = $theDscptnSplt[1];}
		}
	
		if ($geneCategory eq "transcript") {#---transcript range
			
			$strndByTrnscptHsh{$trnscptID} = $geneStrnd;
			$cntgByTrnscptHsh{$trnscptID} = $seq;
			$geneByTrnscptHsh{$trnscptID}= $geneID;
			
			${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$trnscptID}}{"start"} = $featureStart;
			${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$trnscptID}}{"end"} = $featureEnd;

			${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$trnscptID}}{"start"} = $featureStart;
			${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$trnscptID}}{"end"} = $featureEnd;

		} elsif ($geneCategory eq "exon") {#---exon range

			${${$exonRngByGeneHsh{$trnscptID}}{$exonNum}}{"start"} = $featureStart;
			${${$exonRngByGeneHsh{$trnscptID}}{$exonNum}}{"end"} = $featureEnd;
			
		}
		
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;

	my ($NGSIntronRngByGeneHsh_ref, $NGSJunctStrByIntronIDHsh_ref, $NGSIntronRngXS_ref, $NGSIntronStrndHsh_ref) = getIntronFromExonRng(\%exonRngByGeneHsh, \%cntgByTrnscptHsh, \%strndByTrnscptHsh, "Getting the intron boundaries on NGSGtf", "NGS");
	
	#---remove the duplicated NGS introns on the minor isoforms, as many minor isoforms have exactly the same intron as in the prominent isoforms
	my %NGSJunctStrByIntronIDHsh = %{$NGSJunctStrByIntronIDHsh_ref};
	my %tmpJunctStrHsh;
	foreach my $NGSIntronID (keys %NGSJunctStrByIntronIDHsh) {
		my $NGSJunctStr = $NGSJunctStrByIntronIDHsh{$NGSIntronID};
		my @NGSJunctStrSplt = split /:/, $NGSJunctStr;
		my $NGSJunctStrNoPM = join ":", @NGSJunctStrSplt[0..2];
		push @{${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"junctStr"}}, $NGSJunctStr;
		push @{${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"ID"}}, $NGSIntronID;
	}
	
	#---go through each unique location
	my %uniqueNGSJunctStrByIntronIDHsh;
	foreach my $NGSJunctStrNoPM (keys %tmpJunctStrHsh) {
		my $targetIndex = 0;
		for (my $i=0; $i < @{${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"junctStr"}}; $i++) {
			$targetIndex = $i if (${${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"junctStr"}}[$i] =~ m/\:P$/); #---get the prominent isoform index
		}
		my $targetJunctStr = ${${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"junctStr"}}[$targetIndex];
		my $targetID = ${${$tmpJunctStrHsh{$NGSJunctStrNoPM}}{"ID"}}[$targetIndex];
		$uniqueNGSJunctStrByIntronIDHsh{$targetID} = $targetJunctStr;
	}
	
	return (\%SSRngByCntgByWholeGeneHsh, \%strndByTrnscptHsh, \%XStrndRngByCntgByWholeGeneHsh, \%geneByTrnscptHsh, \%exonRngByGeneHsh, $NGSIntronRngByGeneHsh_ref, $NGSJunctStrByIntronIDHsh_ref, $NGSIntronRngXS_ref, $NGSIntronStrndHsh_ref, \%cntgByTrnscptHsh, \%uniqueNGSJunctStrByIntronIDHsh);
}
########################################################################## filterProminentIsoform
sub filterProminentIsoform {

	my %NGSTrnscptRngXSHsh = %{$_[0]};
	my %NGSStrndHsh = %{$_[1]};
	
	my (%NGSPrmntTrnscptRngXSHsh, %NGSPrmntStrndHsh);
	
	foreach my $cntg (keys %NGSTrnscptRngXSHsh) {
		foreach my $NGSFtur (keys %{$NGSTrnscptRngXSHsh{$cntg}}) {
			if ($NGSFtur =~ m/\.0$/) {#---all prominent isoform ends with a 0
				${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"start"} = ${${$NGSTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"start"};
				${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"end"} = ${${$NGSTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"end"};
				$NGSPrmntStrndHsh{$NGSFtur} = $NGSStrndHsh{$NGSFtur};
			}
		}
	}
	
	return (\%NGSPrmntTrnscptRngXSHsh, \%NGSPrmntStrndHsh);
}
########################################################################## findHitAndUnpair
sub findHitAndUnpair {
	
	my %SSHitByRefTrnscptHsh = %{$_[0]};
	my %SSHitByNGSTrnscptHsh = %{$_[1]};
	my %XSHitByRefTrnscptHsh = %{$_[2]};
	my %XSHitByNGSTrnscptHsh = %{$_[3]};
	my %refTrnscptStrndHsh = %{$_[4]};
	my %NGSStrndHsh = %{$_[5]};
	my %refRngXSHsh = %{$_[6]};
	my %NGSPrmntTrnscptRngXSHsh = %{$_[7]};
	my %NGSExonRngByTrnsptHsh = %{$_[8]};
	my $strndSpecificNGS = $_[9];

	my (%refXSUnpairLenHsh, %NGSXSUnpairLenHsh, %refXSHitCountHsh, %NGSXSHitCountHsh, %refBasedNGSHitInfoHsh);

	open (REFXSHITLOG, ">$outDir/hitlog/refHitLogXS.tsv");
	print REFXSHITLOG "\@Hits not nesscessarily on the same strand\n";
	print REFXSHITLOG "\@refFtur\tNGSFtur\toverlapScenerio\tcntg\tStrnd.refStart.refEnd\tNGSStart.NGSEnd\n";

	open (REFSSHITLOG, ">$outDir/hitlog/refHitLogStrndSpecific.tsv");
	print REFSSHITLOG "\@Hits on the same strand\n";
	print REFSSHITLOG "\@refFtur\tNGSFtur\toverlapScenerio\tcntg\tStrnd.refStart.refEnd\tNGSStart.NGSEnd\n";

	open (REFSSUNPAIRGFF, ">$outDir/hitlog/refUnpairStrndSpecific.gff");#---feature that doesnt have hit on other strand but have hit on the same strand
	open (REFXSUNPAIRGFF, ">$outDir/GFF/refUnpairXS.gff"); 
	
	
	foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {
		foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {#--- all ftur on $cntg of refGff
			my $refStart = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"start"};
			my $refEnd = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"end"};
			my $refStrnd = $refTrnscptStrndHsh{$refFtur};
			my $NGSTransfragHit = 0;
			
			my %tmpRefHitPosHsh;
			my @tmpNGSTransfragPosAry;
			
			foreach my $pos ($refStart..($refEnd-1)) {
				$tmpRefHitPosHsh{$pos} = 0;
			}

			if (exists $XSHitByRefTrnscptHsh{$refFtur}) {
				foreach my $NGSFtur (sort {$a cmp $b} keys %{$XSHitByRefTrnscptHsh{$refFtur}}) {
					my $NGSStrnd = $NGSStrndHsh{$NGSFtur};
					my $NGSStart = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"start"};
					my $NGSEnd = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"end"};
					print REFXSHITLOG $refFtur."\t".$NGSFtur."\tscence_${$XSHitByRefTrnscptHsh{$refFtur}}{$NGSFtur}\t$cntg\t$refStrnd.$refStart.$refEnd\t$NGSStrnd.$NGSStart.$NGSEnd\n";
					
					$refXSHitCountHsh{$refFtur}++;
					
					if ($strndSpecificNGS eq 'yes') {
					
						if (exists ${$SSHitByRefTrnscptHsh{$refFtur}}{$NGSFtur}) {

							print REFSSHITLOG $refFtur."\t".$NGSFtur."\tscence_${$SSHitByRefTrnscptHsh{$refFtur}}{$NGSFtur}\t$cntg\t$refStrnd.$refStart.$refEnd\t$NGSStrnd.$NGSStart.$NGSEnd\n";

							$NGSTransfragHit++;

							#----store the NGS coverage info
							foreach my $pos ($NGSStart..$NGSEnd) {
								push @tmpNGSTransfragPosAry, $pos;
								$tmpRefHitPosHsh{$pos} = 1;
							}

						} else {
					
							#DS571233	BCP	gene	434	676	.	-	.	ID=EhSINE3_9
							#DS571233	BCP	SINE	434	676	.	-	.	ID=rna_EhSINE3_9-1;Parent=EhSINE3_9
							#DS571233	BCP	exon	434	676	.	-	.	ID=exon_EhSINE3_9-1;Parent=rna_EhSINE3_9-1

							print REFSSUNPAIRGFF "$cntg\tgeneModelReviser\tgene\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=$refFtur\n";
							print REFSSUNPAIRGFF "$cntg\tgeneModelReviser\trefFrag\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=rna_$refFtur\-1;Parent=$refFtur\n";
							print REFSSUNPAIRGFF "$cntg\tgeneModelReviser\texon\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=exon_$refFtur\-1;Parent=rna_$refFtur\-1\n";
						}
					
					} elsif ($strndSpecificNGS eq 'no') {
						
						$NGSTransfragHit++;

						#----store the NGS coverage info
						foreach my $pos ($NGSStart..$NGSEnd) {
							push @tmpNGSTransfragPosAry, $pos;
							$tmpRefHitPosHsh{$pos} = 1;
						}
					} else {
						die "strndSpecificNGS option invalid: it has to be 'yes' or 'no'\n";
					}
				}
				
			} else {#---unpaired refFutr
				
				my $tmpLength = $refEnd - $refStart;
				$refXSUnpairLenHsh{$refFtur} = $tmpLength;

				print REFXSUNPAIRGFF "$cntg\tgeneModelReviser\tgene\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=$refFtur\n";
				print REFXSUNPAIRGFF "$cntg\tgeneModelReviser\trefFrag\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=rna_$refFtur\-1;Parent=$refFtur\n";
				print REFXSUNPAIRGFF "$cntg\tgeneModelReviser\texon\t$refStart\t$refEnd\t\.\t$refStrnd\t\.\tID=exon_$refFtur\-1;Parent=rna_$refFtur\-1\n";
			}
			
			#----summarize the coverage info
			my $coveredPos = 0;
			my $refLength = $refEnd-$refStart; 
			foreach my $pos ($refStart..($refEnd-1)) {
				$coveredPos += $tmpRefHitPosHsh{$pos};
			}

			my $NGSTransfragCovPct = sprintf "%.2f", 100*$coveredPos/$refLength;

			#----summarize the transfrag extension info
			my $end5TransfragDiff = my $end3TransfragDiff = 'na';
			if (@tmpNGSTransfragPosAry > 0) {
				@tmpNGSTransfragPosAry = sort {$a <=> $b} @tmpNGSTransfragPosAry;
				$end5TransfragDiff = $refStart - $tmpNGSTransfragPosAry[0];
				$end3TransfragDiff = $tmpNGSTransfragPosAry[-1] - $refEnd;
			}

			($end5TransfragDiff, $end3TransfragDiff) = ($end3TransfragDiff, $end5TransfragDiff) if $refStrnd eq '-';

			${$refBasedNGSHitInfoHsh{$refFtur}}{"NGSTransfragHit"} = $NGSTransfragHit;
			${$refBasedNGSHitInfoHsh{$refFtur}}{"NGSTransfragCovPct"} = $NGSTransfragCovPct;
			${$refBasedNGSHitInfoHsh{$refFtur}}{"3EndTransfragDiff"} = $end3TransfragDiff;
			${$refBasedNGSHitInfoHsh{$refFtur}}{"5EndTransfragDiff"} = $end5TransfragDiff;
			${$refBasedNGSHitInfoHsh{$refFtur}}{"refLength"} = $refLength;
			${$refBasedNGSHitInfoHsh{$refFtur}}{"coveredPos"} = $coveredPos;
		}
	}

	close REFSSHITLOG;
	close REFXSHITLOG;
	close REFSSUNPAIRGFF;
	close REFXSUNPAIRGFF;
	
	open (NGSXSHITLOG, ">$outDir/hitlog/NGSHitLogXS.tsv");
	print NGSXSHITLOG "\@Hits not nesscessarily on the same strand, but include on hit on the same strand\n";
	print NGSXSHITLOG "\@refFtur\tNGSFtur\toverlapScenerio\tcntg\tStrnd.refStart.refEnd\tNGSStart.NGSEnd\n";

	open (NGSSSHITLOG, ">$outDir/hitlog/NGSHitLogStrndSpecific.tsv");
	print NGSSSHITLOG "\@Hits on the same strand\n";
	print NGSSSHITLOG "\@refFtur\tNGSFtur\toverlapScenerio\tcntg\tStrnd.refStart.refEnd\tNGSStart.NGSEnd\n";

	open (NGSSSUNPAIRGFF, ">$outDir/GFF/NGSUnpairStrndSpecific.gff");#---feature that doesnt have hit on other strand but have hit on the same strand
	open (NGSXSUNPAIRGFF, ">$outDir/GFF/NGSUnpairXS.gff"); 
	
	foreach my $cntg (sort {$a cmp $b} keys %NGSPrmntTrnscptRngXSHsh) {
		foreach my $NGSFtur (sort {$a cmp $b} keys %{$NGSPrmntTrnscptRngXSHsh{$cntg}}) {#--- all ftur on $cntg of NGSGtf
			my $NGSStart = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"start"};
			my $NGSEnd = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"end"};
			my $NGSStrnd = $NGSStrndHsh{$NGSFtur};

			if (exists $XSHitByNGSTrnscptHsh{$NGSFtur}) {
				foreach my $refFtur (sort {$a cmp $b} keys %{$XSHitByNGSTrnscptHsh{$NGSFtur}}) {

					$NGSXSHitCountHsh{$NGSFtur}++;
					
					my $refStrnd = $refTrnscptStrndHsh{$refFtur};
					my $refStart = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"start"};
					my $refEnd = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"end"};
					print NGSXSHITLOG $refFtur."\t".$NGSFtur."\tscence_${$XSHitByNGSTrnscptHsh{$NGSFtur}}{$refFtur}\t$cntg\t$refStrnd.$refStart.$refEnd\t$NGSStrnd.$NGSStart.$NGSEnd\n";
					
					if (exists ${$SSHitByNGSTrnscptHsh{$NGSFtur}}{$refFtur}) {
						 print NGSSSHITLOG $refFtur."\t".$NGSFtur."\tscence_${$SSHitByNGSTrnscptHsh{$NGSFtur}}{$refFtur}\t$cntg\t$refStrnd.$refStart.$refEnd\t$NGSStrnd.$NGSStart.$NGSEnd\n";
					} else {
						print NGSSSUNPAIRGFF "$cntg\tgeneModelReviser\tgene\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=$NGSFtur\n";
						print NGSSSUNPAIRGFF "$cntg\tgeneModelReviser\ttransfrag\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=rna_$NGSFtur\-1;Parent=$NGSFtur\n";
						foreach my $exonNum (sort {$a <=> $b} keys %{$NGSExonRngByTrnsptHsh{$NGSFtur}}) {
							my $exonStart = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"start"};
							my $exonEnd = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"end"};
							print NGSXSUNPAIRGFF "$cntg\tgeneModelReviser\texon\t$exonStart\t$exonEnd\t\.\t$NGSStrnd\t\.\tID=exon_$NGSFtur\-$exonNum;Parent=rna_$NGSFtur\-1\n";
						}
					}
				}

			} else {
				
				print NGSXSUNPAIRGFF "$cntg\tgeneModelReviser\tgene\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=$NGSFtur\n";
				print NGSXSUNPAIRGFF "$cntg\tgeneModelReviser\ttransfrag\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=rna_$NGSFtur\-1;Parent=$NGSFtur\n";
				foreach my $exonNum (sort {$a <=> $b} keys %{$NGSExonRngByTrnsptHsh{$NGSFtur}}) {
					my $exonStart = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"start"};
					my $exonEnd = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"end"};
					print NGSXSUNPAIRGFF "$cntg\tgeneModelReviser\texon\t$exonStart\t$exonEnd\t\.\t$NGSStrnd\t\.\tID=exon_$NGSFtur\-$exonNum;Parent=rna_$NGSFtur\-1\n";
				}
				my $tmpLength = $NGSEnd - $NGSStart;
				$NGSXSUnpairLenHsh{$NGSFtur} = $tmpLength;
			}
		}
	}

	close NGSSSHITLOG;
	close NGSXSHITLOG;
	close NGSSSUNPAIRGFF;
	close NGSXSUNPAIRGFF;

	return (\%refXSUnpairLenHsh, \%NGSXSUnpairLenHsh, \%refXSHitCountHsh, \%NGSXSHitCountHsh, \%refBasedNGSHitInfoHsh); 
}
########################################################################## summarizeHitAndUnpair
sub summarizeHitAndUnpair {

	my %refXSUnpairLenHsh = %{$_[0]};
	my %NGSXSUnpairLenHsh = %{$_[1]};
	my %refXSHitCountHsh = %{$_[2]};
	my %NGSXSHitCountHsh = %{$_[3]};
	my %refCtgryByGeneHsh = %{$_[4]};
	my %SSPrxmtyByRefTrnscptHsh = %{$_[5]};
	my %SSPrxmtyByNGSTrnscptHsh = %{$_[6]};
	my %XSPrxmtyByRefTrnscptHsh = %{$_[7]};
	my %XSPrxmtyByNGSTrnscptHsh = %{$_[8]};
	my %NGSExonRngByTrnsptHsh = %{$_[9]};
	my %NGSPrmntTrnscptRngXSHsh = %{$_[10]};
	my %NGSCntgByTrnscptHsh	= %{$_[11]};
	my %NGSStrndHsh = %{$_[12]};
	my %NRegionNGSHitHsh = %{$_[13]};
	my %NRegionNGSPrxmtyHsh = %{$_[14]};
	my %cntgEndNGSHitHsh = %{$_[15]};
	my %cntgEndNGSPrxmtyHsh = %{$_[16]};

	#---print the Length info of the non-strand specific unpair of NGS transfrags
	my $minDistance = 100;
	
	my (%allNGSXSUnpairLenFreqHsh, %isolatedNGSXSUnpairLenFreqHsh);

	open (ALLNGSXSUNPAIRLENLOG, ">$outDir/unpariedTransfrag/all.NGSUnpairXS.len.tsv"); 
	open (ISONGSXSUNPAIRLENLOG, ">$outDir/unpariedTransfrag/isolated.$minDistance.NGSUnpairXS.len.tsv"); 
	open (PROXIMALNGSXSUNPAIRLENLOG, ">$outDir/unpariedTransfrag/proximal.$minDistance.NGSUnpairXS.len.tsv"); 

	open (NGSISOLATEGFF, ">$outDir/unpariedTransfrag/isolated.$minDistance.NGSUnpairXS.gff"); 
	open (NGSPROXIMALEGFF, ">$outDir/unpariedTransfrag/proximal.$minDistance.NGSUnpairXS.gff"); 
	open (NGSISONREGCNTGGFF, ">$outDir/unpariedTransfrag/isolated.$minDistance.NReg.Cntg.NGSUnpairXS.gff"); 
	
	open (NGSISONREGCNTGLOG, ">$outDir/unpariedTransfrag/isolated.$minDistance.NReg.Cntg.NGSUnpairXS.len.tsv"); 

	print ALLNGSXSUNPAIRLENLOG join "\t", ("NGSFtur", "exonNumber", "length", "headPrxmty", "tailPrxmty\n");
	print ISONGSXSUNPAIRLENLOG join "\t", ("NGSFtur", "exonNumber", "length", "headPrxmty", "tailPrxmty\n");
	print PROXIMALNGSXSUNPAIRLENLOG join "\t", ("NGSFtur", "exonNumber", "length", "headPrxmty", "tailPrxmty\n");
	print NGSISONREGCNTGLOG "NGSFtur\ttotalExonLength\n";

	my $isolatedRefNRegCntgEnd = 0;

	foreach my $NGSFtur (sort {$NGSXSUnpairLenHsh{$b} <=> $NGSXSUnpairLenHsh{$a}} keys %NGSXSUnpairLenHsh) {
		
		my $exonNumber = keys %{$NGSExonRngByTrnsptHsh{$NGSFtur}};
		
		#---check the proximity of this NGSFutr
		my $headPrxmty = ${${$XSPrxmtyByNGSTrnscptHsh{$NGSFtur}}{"H"}}[0];
		my $tailPrxmty = ${${$XSPrxmtyByNGSTrnscptHsh{$NGSFtur}}{"T"}}[0];

		print ALLNGSXSUNPAIRLENLOG join "\t", ($NGSFtur, $exonNumber, $NGSXSUnpairLenHsh{$NGSFtur}, $headPrxmty, $tailPrxmty."\n");
		$allNGSXSUnpairLenFreqHsh{$NGSXSUnpairLenHsh{$NGSFtur}}++;
		
		#---check for N region and cntg ends
		my $NRegionWithin = my $cntgEndWithin = "no";
		$NRegionWithin = "yes" if (exists $NRegionNGSHitHsh{$NGSFtur});
		$cntgEndWithin = "yes" if (exists $cntgEndNGSHitHsh{$NGSFtur});
		my $NRegionHeadPrmxty = ${${$NRegionNGSPrxmtyHsh{$NGSFtur}}{"H"}}[0];
		my $NRegionTailPrmxty = ${${$NRegionNGSPrxmtyHsh{$NGSFtur}}{"T"}}[0];
		my $cntgEndHeadPrmxty = ${${$cntgEndNGSPrxmtyHsh{$NGSFtur}}{"H"}}[0];
		my $cntgEndTailPrmxty = ${${$cntgEndNGSPrxmtyHsh{$NGSFtur}}{"T"}}[0];
		
		my $NRegCntgFiltered = "no";
		my $comment = "no";
		if (($NRegionWithin eq "yes") or ($cntgEndWithin eq "yes") or (($NRegionHeadPrmxty <= $minDistance) and ($NRegionHeadPrmxty != "-999")) or (($NRegionTailPrmxty <= $minDistance) and ($NRegionTailPrmxty != "-999")) or (($cntgEndHeadPrmxty <= $minDistance) and ($cntgEndHeadPrmxty != "-999")) or (($cntgEndTailPrmxty <= $minDistance) and ($cntgEndTailPrmxty != "-999"))) {
			$NRegCntgFiltered = "yes";
			$comment = "";
			$comment .= "N region within;" if ($NRegionWithin eq "yes");
			$comment .= "at contig end;" if ($cntgEndWithin eq "yes");
			$comment .= "adjacent to contig end;" if ((($cntgEndHeadPrmxty <= $minDistance) and ($cntgEndHeadPrmxty != "-999")) or (($cntgEndTailPrmxty <= $minDistance) and ($cntgEndTailPrmxty != "-999")));
			$comment .= "adjacent to N region;" if ((($NRegionHeadPrmxty <= $minDistance) and ($NRegionHeadPrmxty != "-999")) or (($NRegionTailPrmxty <= $minDistance) and ($NRegionTailPrmxty != "-999")));
		}

		my $cntg = $NGSCntgByTrnscptHsh{$NGSFtur};
		my $NGSStart = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"start"};
		my $NGSEnd = ${${$NGSPrmntTrnscptRngXSHsh{$cntg}}{$NGSFtur}}{"end"};
		my $NGSStrnd = $NGSStrndHsh{$NGSFtur};
		

		if ((($headPrxmty >= $minDistance) or ($headPrxmty == -999)) and (($tailPrxmty >= $minDistance)or ($tailPrxmty == -999))) {#---head and tail are isolated from refFtur at least at $minDistance
			print ISONGSXSUNPAIRLENLOG join "\t", ($NGSFtur, $exonNumber, $NGSXSUnpairLenHsh{$NGSFtur}, $headPrxmty, $tailPrxmty."\n");
			$isolatedNGSXSUnpairLenFreqHsh{$NGSXSUnpairLenHsh{$NGSFtur}}++;
			if ($NRegCntgFiltered eq "no") {
				$isolatedRefNRegCntgEnd++;
				print NGSISONREGCNTGGFF "$cntg\tgeneModelReviser\tgene\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=$NGSFtur\n";
				print NGSISONREGCNTGGFF "$cntg\tgeneModelReviser\ttransfrag\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=rna_$NGSFtur\-1;Parent=$NGSFtur\n";
			}
			
			print NGSISOLATEGFF "$cntg\tgeneModelReviser\tgene\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=$NGSFtur\n";
			print NGSISOLATEGFF "$cntg\tgeneModelReviser\ttransfrag\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=rna_$NGSFtur\-1;Parent=$NGSFtur\n";
			my $totalExonLength = 0;
			foreach my $exonNum (sort {$a <=> $b} keys %{$NGSExonRngByTrnsptHsh{$NGSFtur}}) {
				my $exonStart = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"start"};
				my $exonEnd = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"end"};
				my $exonLength = $exonEnd - $exonStart;
				$totalExonLength += $exonLength;
				print NGSISOLATEGFF "$cntg\tgeneModelReviser\texon\t$exonStart\t$exonEnd\t\.\t$NGSStrnd\t\.\tID=exon_$NGSFtur\-$exonNum;Parent=rna_$NGSFtur\-1\n";
				if ($NRegCntgFiltered eq "no") {
					print NGSISONREGCNTGGFF "$cntg\tgeneModelReviser\texon\t$exonStart\t$exonEnd\t\.\t$NGSStrnd\t\.\tID=exon_$NGSFtur\-$exonNum;Parent=rna_$NGSFtur\-1\n";
				}
			}

			print NGSISONREGCNTGLOG $NGSFtur."\t".$totalExonLength."\n";

		} else {

			print PROXIMALNGSXSUNPAIRLENLOG join "\t", ($NGSFtur, $exonNumber, $NGSXSUnpairLenHsh{$NGSFtur}, $headPrxmty, $tailPrxmty."\n");

			print NGSPROXIMALEGFF "$cntg\tgeneModelReviser\tgene\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=$NGSFtur\n";
			print NGSPROXIMALEGFF "$cntg\tgeneModelReviser\ttransfrag\t$NGSStart\t$NGSEnd\t\.\t$NGSStrnd\t\.\tID=rna_$NGSFtur\-1;Parent=$NGSFtur\n";
			foreach my $exonNum (sort {$a <=> $b} keys %{$NGSExonRngByTrnsptHsh{$NGSFtur}}) {
				my $exonStart = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"start"};
				my $exonEnd = ${${$NGSExonRngByTrnsptHsh{$NGSFtur}}{$exonNum}}{"end"};
				print NGSPROXIMALEGFF "$cntg\tgeneModelReviser\texon\t$exonStart\t$exonEnd\t\.\t$NGSStrnd\t\.\tID=exon_$NGSFtur\-$exonNum;Parent=rna_$NGSFtur\-1\n";
			}
		}
	}
	
	close ALLNGSXSUNPAIRLENLOG;
	close ISONGSXSUNPAIRLENLOG;
	close NGSISOLATEGFF;
	close NGSPROXIMALEGFF;
	close PROXIMALNGSXSUNPAIRLENLOG;
	close NGSISONREGCNTGGFF;

	print $isolatedRefNRegCntgEnd." transfrags were found to be isolated from all ref anno, NReg and Cntg ends.\n";

	GNUPlotXYScatterWithLines(\%allNGSXSUnpairLenFreqHsh, "$outDir/unpariedTransfrag/all.NGSUnpairXS.len.pdf", "$outDir/unpariedTransfrag/all.NGSUnpairXS.len.dat", "Lenght (nt)", "Frequency", "linear", "linear", "Length unpaired NGS features");
	GNUPlotXYScatterWithLines(\%allNGSXSUnpairLenFreqHsh, "$outDir/unpariedTransfrag/isolated.$minDistance.NGSUnpairXS.len.pdf", "$outDir/unpariedTransfrag/isolated.$minDistance.NGSUnpairXS.len.dat", "Lenght (nt)", "Frequency", "linear", "linear", "Length unpaired NGS features");
	
}
########################################################################## checkOverlapAndProximity
sub checkOverlapAndProximity {
#		
#	dependence: printProgressScale, updateProgressBar
#
#					The 7 scenes of overlapping and proximity 
#
#
#     case 0: complete overlapp (($refStart == $qryStart) && ($refEnd == $qryEnd))
#			
#     case 1: overlapHead         case 2: overlapTail	        case 3: cover		     case 4: within		case 5: prxmtyTail	     case 6: prxmtyHead
#
#ref  |--------|		        |---------|	            |-------------|	              |-----|					|-----|				  			     	   |-------|
#Qry  	<=========>	         <==========>		          <=========>	            <==========>						<==========>	      <==========>
#
#     ($refStart<$qryStart)&&	   ($refStart>=$qryStart)&&	  ($refStart<$qryStart)&&   ($refStart>$qryStart)&& ($refEnd<=$qryStart)&&		($refStart>$qryStart)&&
#     ($refEnd>=$qryStart)&&	   ($refStart<=$qryEnd)&&	  ($refEnd>$qryEnd)	    ($refEnd<$qryEnd)			($refEnd<$qryEnd)			($refStart>=$qryEnd)
#     ($refEnd<=$qryEnd)	   ($refEnd>$qryEnd)												 
#
	#---incoming variables
	my %refRngXSHsh = %{$_[0]}; 
	my %qryRngXSHsh = %{$_[1]}; 
	my %refTrnscptStrndHsh = %{$_[2]};
	my %qryStrndHsh = %{$_[3]};
	my $strToPrint = $_[4];
	my $checkPrxmty = $_[5];
	
	#---outgoing variables
	my (%SSHitByRefTrnscptHsh, %SSHitByQryHsh, %XSHitByRefTrnscptHsh, %XSHitByQryHsh, %SSPrxmtyByRefTrnscptHsh, %SSPrxmtyByQryHsh, %XSPrxmtyByRefTrnscptHsh, %XSPrxmtyByQryHsh);

	#---on screen progress scale
	printProgressScale("\n$strToPrint", 50, 10);
	my $progressCntgNum = 0;

	#---make a tmpHsh to contain all cntgs the have either ref and qry
	my %tmpCntgHsh;
	foreach my $cntg (keys %refRngXSHsh) {$tmpCntgHsh{$cntg}++;}
	foreach my $cntg (keys %qryRngXSHsh) {$tmpCntgHsh{$cntg}++;}

	my $totalCntgNum = keys %tmpCntgHsh;

	foreach my $cntg (sort {$a cmp $b} keys %tmpCntgHsh) {

		#---update the on-screen progress
		$progressCntgNum++;
		updateProgressBar("$cntg", $progressCntgNum, $totalCntgNum, 50, 10);
		
		my (%tmpSSPrxmtyByRefTrnscptHsh, %tmpSSPrxmtyByQryHsh, %tmpXSPrxmtyByRefTrnscptHsh, %tmpXSPrxmtyByQryHsh);

		if ((exists $qryRngXSHsh{$cntg}) and (exists $refRngXSHsh{$cntg})) {#---if there are both ref and qry can both be found on cntg
			foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
				my $refStart = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"start"};
				my $refEnd = ${${$refRngXSHsh{$cntg}}{$refFtur}}{"end"};
				foreach  my $qryFtur (sort {$a cmp $b} keys %{$qryRngXSHsh{$cntg}}) {#--- all ftur on the $strnd of $cntg of QryGtf
					
					my $sameStrnd = "no";
					if (($refTrnscptStrndHsh{$refFtur} eq $qryStrndHsh{$qryFtur}) or ($qryStrndHsh{$qryFtur} eq ".") or ($refTrnscptStrndHsh{$refFtur} eq ".")){
						$sameStrnd = "yes";
					}
					
					my $qryStart = ${${$qryRngXSHsh{$cntg}}{$qryFtur}}{"start"};
					my $qryEnd = ${${$qryRngXSHsh{$cntg}}{$qryFtur}}{"end"};
					
					if  (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0

						${$XSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 0;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 0;

						if ($sameStrnd eq "yes") {
							${$SSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 0;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 0;
						}
						
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1		

						${$XSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 1;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 1;

						if ($sameStrnd eq "yes") {
							${$SSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 1;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 1;
						}
											
					} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2					

						${$XSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 2;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 2;

						if ($sameStrnd eq "yes") {
							${$SSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 2;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 2;
						}
					
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3		

						${$XSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 3;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 3;

						if ($sameStrnd eq "yes") {
							${$SSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 3;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 3;
						}
					
					} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4						

						${$XSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 4;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 4;

						if ($sameStrnd eq "yes") {
							${$SSHitByRefTrnscptHsh{$refFtur}}{$qryFtur} = 4;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 4;
						}
					
					#------Proximity with ref's tail proximal to qry's head
					} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

						if ($checkPrxmty eq "yes") {
							my $tmpPrmxty = $qryStart - $refEnd;
							${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$qryFtur} = $tmpPrmxty;
							${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur} = $tmpPrmxty;

							if ($sameStrnd eq "yes") {	
								${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$qryFtur} = $tmpPrmxty;
								${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur} = $tmpPrmxty;
							}
						}

					#------Proximity with ref's head proximal to qry's tail
					} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

						if ($checkPrxmty eq "yes") {
							my $tmpPrmxty = $refStart - $qryEnd;
							${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$qryFtur} = $tmpPrmxty;
							${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur} = $tmpPrmxty;

							if ($sameStrnd eq "yes") {	
								${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$qryFtur} = $tmpPrmxty;
								${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur} = $tmpPrmxty;
							}
						}
						
					} else {#---BUG! possibly other scene?
						print "refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";
						die "Unexpected overlapping scene between $refFtur and $qryFtur. It's a Bug. Program qutting.\n";
					}
				}
			}#---end of foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
		} #---end of if (exists $qryRngXSHsh{$cntg}) {
		
		#---find the closest proximity for all refs
		if ($checkPrxmty eq "yes") {
		
			#---for all ref based info
			foreach my $refFtur (keys %{$refRngXSHsh{$cntg}}) {
	
				#---in cases if the proximity are edges
				${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}); 
				${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}); 
				${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}); 
				${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}); 
	
				#---for all XS heads
				foreach my $qryFtur (sort {${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$a} <=> ${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$b}} keys %{${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}) {
					@{${$XSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}} = (${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all XS tails
				foreach my $qryFtur (sort {${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$a} <=> ${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$b}} keys %{${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}) {
					@{${$XSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}} = (${${$tmpXSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all SS heads
				foreach my $qryFtur (sort {${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$a} <=> ${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$b}} keys %{${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}) {
					@{${$SSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}} = (${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"H"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all SS tails
				foreach my $qryFtur (sort {${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$a} <=> ${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$b}} keys %{${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}) {
					@{${$SSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}} = (${${$tmpSSPrxmtyByRefTrnscptHsh{$refFtur}}{"T"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}
				
			}#---end of foreach my $refFtur (keys %refRngXSHsh)
			
			#---for all qry based info
			foreach my $qryFtur (keys %{$qryRngXSHsh{$cntg}}) {
	
				#---in cases if the proximity are edges
				${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}); 
				${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}); 
				${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}); 
				${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}); 
	
				#---for all XS heads
				foreach my $refFtur (sort {${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$a} <=> ${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$b}} keys %{${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}) {
					@{${$XSPrxmtyByQryHsh{$qryFtur}}{"H"}} = (${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all XS tails
				foreach my $refFtur (sort {${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$a} <=> ${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$b}} keys %{${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}) {
					@{${$XSPrxmtyByQryHsh{$qryFtur}}{"T"}} = (${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all SS heads
				foreach my $refFtur (sort {${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$a} <=> ${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$b}} keys %{${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}) {
					@{${$SSPrxmtyByQryHsh{$qryFtur}}{"H"}} = (${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all SS tails
				foreach my $refFtur (sort {${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$a} <=> ${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$b}} keys %{${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}) {
					@{${$SSPrxmtyByQryHsh{$qryFtur}}{"T"}} = (${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}
				
			}#---end of foreach my $qryFtur (keys %qryRngXSHsh) {
		}
		
	}#---end foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {

	print "\n\n";
	
	return (\%SSHitByRefTrnscptHsh, \%SSHitByQryHsh, \%XSHitByRefTrnscptHsh, \%XSHitByQryHsh, \%SSPrxmtyByRefTrnscptHsh, \%SSPrxmtyByQryHsh, \%XSPrxmtyByRefTrnscptHsh, \%XSPrxmtyByQryHsh);
}
########################################################################## GNUPLOTAryHistogram
sub GNUPLOTAryHistogram {

	my @numAry = @{$_[0]};
	my $fileName = $_[1];
	
	print "Running GNUPLOTAryHistogram for $fileName.\n";
	
	#---calculate the optimal bin
	@numAry = sort {$a <=> $b} @numAry;
	my $tail5PctPos = int ($#numAry*0.05);
	my $upper95PctPos = $#numAry - $tail5PctPos;
	my $lower95PctPos = $tail5PctPos;
	my $binWidth = int (($numAry[$upper95PctPos] - $numAry[$lower95PctPos])/100);
	$binWidth = 1 if ($binWidth < 1);

	#---creat a tmp file
	open (TMPFILE, ">tmp.dat");
	for (@numAry) {print TMPFILE $_."\n";}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	bin_width = $binWidth
	bin_number(x) = floor(x/bin_width)
	rounded(x) = bin_width * (bin_number(x) + 0.5)
	UNITY = 1
	set output "| ps2pdf - $fileName";
	set xlabel "value";
	set ylabel "frequency";
   	plot 'tmp.dat' u (rounded(\$1)):(UNITY) t 'bin at $binWidth' smooth frequency w histeps
EOPLOT
	close(GNUPLOT);
	rmtree(['./tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## GNUPlotXYScatterWithLines
sub GNUPlotXYScatterWithLines {

	my %XYHsh = %{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $xscale = $_[5];
	my $yscale = $_[6];
	my $title = $_[7];
	
	my $GNULogXCmd = "";
	$GNULogXCmd = "set logscale x" if ($xscale eq "log");
	my $GNULogYCmd = "";
	$GNULogYCmd = "set logscale y" if ($yscale eq "log");

	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotXYScatterWithLines for $fileName.\n";
	
	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for my $x (sort {$a <=> $b} keys %XYHsh) {
		print TMPFILE $x."\t".$XYHsh{$x}."\n";
	}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath";
	unset logscale x; 
	unset logscale y; 
	$GNULogXCmd;
	$GNULogYCmd;
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
	set nokey;
   	plot '$plotDataPath' using 1:2 with lines;
EOPLOT
	close(GNUPLOT);
	#rmtree(['tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## getIntronFromExonRng
sub getIntronFromExonRng {

	#---incoming variables
	my %exonRngByGeneHsh = %{$_[0]};
	my %cntgByGeneHsh = %{$_[1]};
	my %strndByGeneHsh = %{$_[2]};
	my $strToPrint = $_[3];
	my $refOrNGS = $_[4];
	
	#---outgoing variables
	my (%intronRngByGeneHsh, %junctStrByIntronIDHsh, %XStrndRngByCntgByIntronID, %strndByIntronIDHsh);
	
	#---on screen progress scale
	printProgressScale("\n$strToPrint", 50, 10);
	my $totalGeneNum = keys %exonRngByGeneHsh;
	my $progressGeneNum = 0;
	
	#---go through each gene, see if there's intron, and store the ranges
	foreach my $geneID (sort {$a cmp $b} keys %exonRngByGeneHsh) {

		#---update the on-screen progress
		$progressGeneNum++;
		updateProgressBar("$geneID", $progressGeneNum, $totalGeneNum, 50, 10);
		
		#---get the cntg and strnd
		my $cntg = $cntgByGeneHsh{$geneID};
		my $strnd = $strndByGeneHsh{$geneID};
		
		#---get all exon bounds
		my @tmpBoundAry;
		foreach my $exon (keys %{$exonRngByGeneHsh{$geneID}}) {
			push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"start"};
			push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"end"};
		}
		
		#---if more than one exon
		if (@tmpBoundAry > 2) { 
			#---sort the bounds
			my @sortedTmpBoundAry = sort {$a <=> $b} @tmpBoundAry;
			#---get the intronBounds, go to all odd number indexes and take itself and the +1 index values, i.e. $i = 1, 3, 5 if there're 0,1,2,3,4,5,6,7  
			my $intronNum = 0;
			for (my $i=1; $i<(@sortedTmpBoundAry-1); $i=$i+2) {
				$intronNum++;
				my $intronStart = $sortedTmpBoundAry[$i]+1;
				my $intronEnd = $sortedTmpBoundAry[$i+1]-1;
				my $intronID = $geneID.":".$intronNum;
				
				#---determine the junct is on ref, prominent or minor isoform
				my $prominentMinorRef;
				if ($refOrNGS eq "ref") {
					$prominentMinorRef = "R";
				} else {
					my @geneIDSplt = split /\./, $geneID;
					my $isoformNum = $geneIDSplt[-1];
					if ($isoformNum == 0) {#---prominent isoform
						$prominentMinorRef = "P";
					} else {
						$prominentMinorRef = "M";
					}
				}
				my $junctStr = join ":", ($cntg, $intronStart, $intronEnd, $prominentMinorRef);
				${${$intronRngByGeneHsh{$geneID}}{$intronNum}}{"start"} = $intronStart;
				${${$intronRngByGeneHsh{$geneID}}{$intronNum}}{"end"} = $intronEnd;
				${${$XStrndRngByCntgByIntronID{$cntg}}{$intronID}}{"start"} = $intronStart;
				${${$XStrndRngByCntgByIntronID{$cntg}}{$intronID}}{"end"} = $intronEnd;
				$strndByIntronIDHsh{$intronID} = $strnd;
				$junctStrByIntronIDHsh{$intronID} = $junctStr;
			}
		}
	}

	print "\n\n";

	return (\%intronRngByGeneHsh, \%junctStrByIntronIDHsh, \%XStrndRngByCntgByIntronID, \%strndByIntronIDHsh);
	
}
########################################################################## updateProgressBar
sub updateProgressBar {
	
	my $strToPrint = $_[0];
	my $progressCount = $_[1];
	my $totalCount = $_[2];
	my $scaleMax = $_[3];
	my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
	my $progressPct = int (($progressCount/$totalCount)*$scaleMax);

	my $progressBar = "|";
	for my $i (1..$progressPct) {$progressBar .= ">";}
	for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
	$progressBar .= "|";

	my $extraWhiteSpaceStr = "";
	for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
	print $progressBar.$strToPrint.$extraWhiteSpaceStr."\r";

}
########################################################################## printProgressScale
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
########################################################################## summarizeRefJunctOverlapScenarios
sub summarizeRefJunctOverlapScenarios {
	
	#---incoming vars
	my %SSNGSIntronHitByRefHsh = %{$_[0]};
	my %refJunctStrByIntronIDHsh = %{$_[1]};
	my %NGSJunctStrByIntronIDHsh = %{$_[2]};
	my %SSBEDIntronHitByRefHsh = %{$_[3]};
	my %refJunctIntronBoundCovHsh = %{$_[4]};
	
	#---outgoing vars
	my (%allScenarioByRefIntronIDHsh);

	#--var intitialization
	my $refIntronNum = keys %refJunctStrByIntronIDHsh;
	my $NGSIntronNum = keys %NGSJunctStrByIntronIDHsh;

	#---summarize the ref hits
	open (CONFIRM, ">$outDir/RefJunctions/confirmed.RefJunctions.tsv");
	open (MODIFY, ">$outDir/RefJunctions/modified.RefJunctions.tsv");
	print CONFIRM join '', ((join "\t", ("scenario", "scenarioNumber", "startDiff", "endDiff", "refIntronID", "NGSIntronID", "refJunctStr", "NGSJunctStr")), "\n");
	print MODIFY join '', ((join "\t", ("scenario", "scenarioNumber", "startDiff", "endDiff", "refIntronID", "NGSIntronID", "refJunctStr", "NGSJunctStr")), "\n");
	print "Printing overlapping introns\n";

	foreach my $refIntronID (sort {$a cmp $b} keys %SSNGSIntronHitByRefHsh) {
		my $refJunctStr = $refJunctStrByIntronIDHsh{$refIntronID};
		my @refJunctStrSplt = split /:/, $refJunctStr;
		my $refStart = $refJunctStrSplt[1];
		my $refEnd = $refJunctStrSplt[2];
		
		foreach my $NGSIntronID (sort {${$SSNGSIntronHitByRefHsh{$refIntronID}}{$a} <=> ${$SSNGSIntronHitByRefHsh{$refIntronID}}{$b}} keys %{$SSNGSIntronHitByRefHsh{$refIntronID}}) {
			my $NGSJunctStr = $NGSJunctStrByIntronIDHsh{$NGSIntronID};
			my @NGSJunctStrSplt = split /:/, $NGSJunctStr;
			my $NGSStart = $NGSJunctStrSplt[1];
			my $NGSEnd = $NGSJunctStrSplt[2];

			my $startDiff = $refStart-$NGSStart;
			my $endDiff = $refEnd-$NGSEnd;
			
			my $scenario = ${$SSNGSIntronHitByRefHsh{$refIntronID}}{$NGSIntronID};
			my $printStr = join "\t", ($scenario, $startDiff, $endDiff, $refIntronID, $NGSIntronID, $refJunctStr, $NGSJunctStr."\n");
			
			if ($NGSJunctStr =~ m/P$/) {#---NGS is a prominent Jucntion
				if ($scenario == 0) {#---complete exact overlap

					if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
						$allScenarioByRefIntronIDHsh{$refIntronID} = "confirm_P";
					} else {
						print "confirm_P: ".$refIntronID." mapped to multiple scenarios\n";
						$allScenarioByRefIntronIDHsh{$refIntronID} .= ":confirm_P";
					}
									
					print CONFIRM $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
	
				} else {
	
					if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
						$allScenarioByRefIntronIDHsh{$refIntronID} = "modify_P";
					} elsif ($allScenarioByRefIntronIDHsh{$refIntronID}) {
						print "modify_P: ".$refIntronID." mapped to multiple scenarios\n";
						$allScenarioByRefIntronIDHsh{$refIntronID} .= ":modify_P";
					}
	
					print MODIFY $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
				}
			}
		}
	}

	#---check the minor isoforms and check BED overlaping, if not overlapping with anything, put it to unknown.
	open (UNKNOWN, ">$outDir/RefJunctions/unknown.RefJunctions.tsv");
	print UNKNOWN "scenario\trefIntronID\trefJunctStr\n";

	foreach my $refIntronID (keys %refJunctStrByIntronIDHsh) {
	
		if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {#---no record in the prominent isoforms

			my $refJunctStr = $refJunctStrByIntronIDHsh{$refIntronID};
			my @refJunctStrSplt = split /:/, $refJunctStr;
			my $refStart = $refJunctStrSplt[1];
			my $refEnd = $refJunctStrSplt[2];

			if (exists $SSNGSIntronHitByRefHsh{$refIntronID}) {#----no record in the prominent isoforms but have overlapping, i.e. minor isoforms
			
				foreach my $NGSIntronID (sort {$a cmp $b} keys %{$SSNGSIntronHitByRefHsh{$refIntronID}}) {
				
					my $NGSJunctStr = $NGSJunctStrByIntronIDHsh{$NGSIntronID};
					my @NGSJunctStrSplt = split /:/, $NGSJunctStr;
					my $NGSStart = $NGSJunctStrSplt[1];
					my $NGSEnd = $NGSJunctStrSplt[2];

					my $startDiff = $refStart-$NGSStart;
					my $endDiff = $refEnd-$NGSEnd;
			
					my $scenario = ${$SSNGSIntronHitByRefHsh{$refIntronID}}{$NGSIntronID};
					my $printStr = join "\t", ($scenario, $startDiff, $endDiff, $refIntronID, $NGSIntronID, $refJunctStr, $NGSJunctStr."\n");
			
					die "unexpected scenario, no reason to be prominent, quitting.\n" if ($NGSJunctStr =~ m/P$/);
			
					if ($NGSJunctStr =~ m/M$/) {#---NGS is a prominent Jucntion
						if ($scenario == 0) {#---complete exact overlap

							if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
								$allScenarioByRefIntronIDHsh{$refIntronID} = "confirm_M";
							} else {
								print "confirm_M: ".$refIntronID." mapped to multiple scenarios\n";
								$allScenarioByRefIntronIDHsh{$refIntronID} .= ":confirm_M";
							}
									
							print CONFIRM $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
	
						} else {
	
							if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
								$allScenarioByRefIntronIDHsh{$refIntronID} = "modify_M";
							} else {
								print "modify_M: ".$refIntronID." mapped to multiple scenarios\n";
								$allScenarioByRefIntronIDHsh{$refIntronID} .= ":modify_M";
							}
	
							print MODIFY $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
						}
					}
				}

			} elsif (exists $SSBEDIntronHitByRefHsh{$refIntronID}) {#----no with all isoforms, but may still overlap with the the BED junctions
			
				foreach my $BEDJunctStr (sort {$a cmp $b} keys %{$SSBEDIntronHitByRefHsh{$refIntronID}}) {
					
					my $BEDIntronID = $BEDJunctStr;
					my @BEDJunctStrSplt = split /:/, $BEDJunctStr;
					my $BEDStart = $BEDJunctStrSplt[1];
					my $BEDEnd = $BEDJunctStrSplt[2];

					my $startDiff = $refStart-$BEDStart;
					my $endDiff = $refEnd-$BEDEnd;
			
					my $scenario = ${$SSBEDIntronHitByRefHsh{$refIntronID}}{$BEDIntronID};
					
					my $printStr = join "\t", ($scenario, $startDiff, $endDiff, $refIntronID, $BEDIntronID, $refJunctStr, $BEDJunctStr."\n");
			
					die "unexpected scenario, no reason to be prominent, quitting.\n" if ($BEDJunctStr !~ m/B$/);
			
					if ($scenario == 0) {#---complete exact overlap
					
						if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
							$allScenarioByRefIntronIDHsh{$refIntronID} = "confirm_B";
						} else {
							print "confirm_B: ".$refIntronID." mapped to multiple scenarios\n";
							$allScenarioByRefIntronIDHsh{$refIntronID} .= ":confirm_B";
						}
								
						print CONFIRM $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
	
					} else {
	
						#if (not exists $allScenarioByRefIntronIDHsh{$refIntronID}) {
						#	$allScenarioByRefIntronIDHsh{$refIntronID} = "modify_B";
						#} else {
						#	print "modify_B: ".$refIntronID." mapped to multiple scenarios\n";
						#	$allScenarioByRefIntronIDHsh{$refIntronID} .= ":modify_B";
						#}
	
						#print MODIFY $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$printStr;
				
						$allScenarioByRefIntronIDHsh{$refIntronID} = "unknown";
						print UNKNOWN $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$refIntronID."\t".$refJunctStr."\n";

					}
				}
			
			} else {#---no hit on both prominent and minor isoforms
				$allScenarioByRefIntronIDHsh{$refIntronID} = "unknown";
				print UNKNOWN $allScenarioByRefIntronIDHsh{$refIntronID}."\t".$refIntronID."\t".$refJunctStr."\n";
			}
		}
	}

	close UNKNOWN;
	close CONFIRM;
	close MODIFY;
	
	my %scenarioCountHsh;
	open (OVERALL, ">$outDir/RefJunctions/all.RefJunctions.tsv");
	print OVERALL join "\t", ("scenario", "refIntronID", "junctStr", "startRatio", "endRatio", "startOutterCov", "endOutterCov\n");
	foreach my $refIntronID (sort {$allScenarioByRefIntronIDHsh{$a} cmp $allScenarioByRefIntronIDHsh{$b}} keys %allScenarioByRefIntronIDHsh) {
		my $scenario = $allScenarioByRefIntronIDHsh{$refIntronID};
		my $junctStr = $refJunctStrByIntronIDHsh{$refIntronID};
		my $startOutterCov = my $startInnerCov = my $endOutterCov = my $endInnerCov = 0;
		if (exists $refJunctIntronBoundCovHsh{$refIntronID}) {
			$startOutterCov = ${${$refJunctIntronBoundCovHsh{$refIntronID}}{"start"}}{"out"};
			$startInnerCov = ${${$refJunctIntronBoundCovHsh{$refIntronID}}{"start"}}{"in"};
			$endOutterCov = ${${$refJunctIntronBoundCovHsh{$refIntronID}}{"end"}}{"out"};
			$endInnerCov = ${${$refJunctIntronBoundCovHsh{$refIntronID}}{"end"}}{"in"};
		}
		
		my $startRatio = my $endRatio = 9999;
		$startRatio = sprintf '%.2f', ($startOutterCov/$startInnerCov) if ($startInnerCov > 0);
		$endRatio = sprintf '%.2f', ($endOutterCov/$endInnerCov) if ($endInnerCov > 0);
		print OVERALL join "\t", ($scenario, $refIntronID, $junctStr, $startRatio, $endRatio, $startOutterCov, $endOutterCov."\n");

		$scenarioCountHsh{$scenario}++;
	}
	close OVERALL;
	
	#---go through all the scenarios
	my $totalRefIntronNum = keys %allScenarioByRefIntronIDHsh;
	open (STATS, ">$outDir/stats.log.txt");
	print STATS ("Statistics on NGS and Ref intron overlapping\n");
	print STATS join "\t", ("percetage", "scenario", "count\n");

	print ("\nStatistics on NGS and Ref intron overlapping\n");
	print join "\t", ("percetage", "scenario", "count\n");
	
	foreach my $scenario (sort {$a cmp $b} keys %scenarioCountHsh) {
		my $count = $scenarioCountHsh{$scenario};
		my $tmpPct = sprintf '%.2f', ($count/$totalRefIntronNum)*100;

		print STATS join "\t", ($tmpPct, $scenario, $count."\n");
		print join "\t", ($tmpPct, $scenario, $count."\n");
	}
	close STATS;

	return \%allScenarioByRefIntronIDHsh;
}
########################################################################## summarizeNGSJunctOverlapScenarios
sub summarizeNGSJunctOverlapScenarios {

	my %SSNGSIntronHitByRefTrnscptHsh = %{$_[0]};
	my %SSNGSIntronPrxmtyByRefTrnscptHsh = %{$_[1]};
	my %SSRefTrnscptPrxmtyByNGSIntronHsh = %{$_[2]};
	my %SSRefTrnscptHitByNGSIntronHsh = %{$_[3]};
	my %SSRefIntronHitByNGSHsh = %{$_[4]};
	my %refJunctStrByIntronIDHsh = %{$_[5]};
	my %uniqueNGSJunctStrByIntronIDHsh = %{$_[6]};
	my %allScenarioByRefIntronIDHsh = %{$_[7]};
	my %refCtgryByGeneHsh = %{$_[8]};
	my %refTrnscptStrndHsh = %{$_[9]};
	my %NGSIntronStrndHsh = %{$_[10]};
	my %XSRefTrnscptHitByNGSIntronHsh = %{$_[11]};
	my $NGSJRefPrxmtyLimit = $_[12];

	my $test = keys %SSNGSIntronPrxmtyByRefTrnscptHsh;

	printProgressScale("Summarizing NGS junctions hit and proximity", 50);
	my $totalNGSIntronNum = keys %uniqueNGSJunctStrByIntronIDHsh;
	my $procNGSIntronNum = 0;
	
	my %NGSJunctOnRefCommentHsh;
	
	open (NGSINTRONINFO, ">$outDir/NGSJunctions/all.NGSJunctions.tsv");
	print NGSINTRONINFO join "", (join "\t", ("NGSIntronID", "NGSJunctStr", "comment", "prominentMinor", "refIntronHitID", "refIntronHitScenerio", "refTrnscptSenseHitID", "refTrnscptSenseHitType", "refTrnscptAntiSenseHitID", "refTrnscptAntiSenseHitType", "refTrnscptHeadPrxmtyID", "refTrnscptHeadPrxmtyDistance", "refTrnscptTailPrxmtyID", "refTrnscptTailPrxmtyDistance")), "\n";

	my %scenarioCountHsh;
	
	#---go through each NGSIntronID and print all info about each NGS intron
	foreach my $NGSIntronID (sort {$uniqueNGSJunctStrByIntronIDHsh{$a} cmp $uniqueNGSJunctStrByIntronIDHsh{$b}} keys %uniqueNGSJunctStrByIntronIDHsh) {
		
		$procNGSIntronNum++;
		updateProgressBar($NGSIntronID, $procNGSIntronNum, $totalNGSIntronNum, 50, 10);
		
		my $NGSJunctStr = $uniqueNGSJunctStrByIntronIDHsh{$NGSIntronID};
		my $prominentMinor = "minor";
		$prominentMinor = "prominent" if ($NGSJunctStr =~ m/\:P$/);
		
		#---check if overlapping with ref Intron
		my $refIntronHitID = "NA";
		my $refIntronHitScenerio = "NA";
		if (exists $SSRefIntronHitByNGSHsh{$NGSIntronID}) {#---in cases of overlapping ref ftur, use ";" delimited string
			$refIntronHitID = "";
			$refIntronHitScenerio = "";
			foreach my $refIntronID (keys %{$SSRefIntronHitByNGSHsh{$NGSIntronID}}) {
				$refIntronHitID .= $refJunctStrByIntronIDHsh{$refIntronID}.";";
				$refIntronHitScenerio .= $allScenarioByRefIntronIDHsh{$refIntronID}.";";
			}
		}
		
		#---check if SENSE overlapping with ref Trnscrpt
		my $refTrnscptSenseHitID = "NA";
		my $refTrnscptSenseHitType = "NA";
		if (exists $SSRefTrnscptHitByNGSIntronHsh{$NGSIntronID}) {#---in cases of overlapping ref ftur, use ";" delimited string
			$refTrnscptSenseHitID = "";
			$refTrnscptSenseHitType = "";
			foreach my $refTrnscptID (keys %{$SSRefTrnscptHitByNGSIntronHsh{$NGSIntronID}}) {
				$refTrnscptSenseHitType .= $refCtgryByGeneHsh{$refTrnscptID}.";";
				$refTrnscptSenseHitID .= $refTrnscptID.";";
				push @{$NGSJunctOnRefCommentHsh{$refTrnscptID}}, "prominentJunctInORF" if (($prominentMinor eq "prominent") and ($refIntronHitScenerio eq "NA")); #---hit on transcript but not hit on intron
			}
		}
		
		#---check if ANTISENSE overlapping with ref Trnscrpt
		my $refTrnscptAntiSenseHitID = "NA";
		my $refTrnscptAntiSenseHitType = "NA";
		if (exists $XSRefTrnscptHitByNGSIntronHsh{$NGSIntronID}) {#---in cases of overlapping ref ftur, use ";" delimited string
			my @tmpRefIDAry;
			foreach my $refTrnscptID (keys %{$XSRefTrnscptHitByNGSIntronHsh{$NGSIntronID}}) {
				push @tmpRefIDAry, $refTrnscptID if  ($refTrnscptStrndHsh{$refTrnscptID} ne $NGSIntronStrndHsh{$NGSIntronID});
			}
			if (@tmpRefIDAry > 0) {
				$refTrnscptAntiSenseHitID = "";
				$refTrnscptAntiSenseHitType = "";
				foreach my $refTrnscptID (@tmpRefIDAry) {
					$refTrnscptAntiSenseHitType .= $refCtgryByGeneHsh{$refTrnscptID}.";";
					$refTrnscptAntiSenseHitID .= $refTrnscptID.";";
				}
			}
		}

		#---check th proximity with ref transcript
		my $refTrnscptHeadPrxmtyID = "NA";
		my $refTrnscptHeadPrxmtyDistance = "NA";
		my $refTrnscptTailPrxmtyID = "NA";
		my $refTrnscptTailPrxmtyDistance = "NA";
		my $withinPrxmty = "no";
		my $prxmtyType = "";
		
		if (exists ${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"H"}) {
			$refTrnscptHeadPrxmtyDistance = ${${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"H"}}[0];
			$refTrnscptHeadPrxmtyID = ${${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"H"}}[1];
			if (($refTrnscptHeadPrxmtyDistance > 0) and ($refTrnscptHeadPrxmtyDistance <= $NGSJRefPrxmtyLimit)) {
				$withinPrxmty = "yes";
				$prxmtyType .= $refCtgryByGeneHsh{$refTrnscptHeadPrxmtyID}.";";
				push @{$NGSJunctOnRefCommentHsh{$refTrnscptHeadPrxmtyID}}, "prominentJunctAt3End$NGSJRefPrxmtyLimit" if (($prominentMinor eq "prominent") and ($refTrnscptSenseHitID eq "NA") and ($refTrnscptAntiSenseHitID eq "NA") and ($refIntronHitID eq "NA"));
			}
		}
		
		if (exists ${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"T"}) {
			$refTrnscptTailPrxmtyDistance = ${${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"T"}}[0];
			$refTrnscptTailPrxmtyID = ${${$SSRefTrnscptPrxmtyByNGSIntronHsh{$NGSIntronID}}{"T"}}[1];
			if (($refTrnscptTailPrxmtyDistance > 0) and ($refTrnscptTailPrxmtyDistance <= $NGSJRefPrxmtyLimit)) {
				$withinPrxmty = "yes";
				$prxmtyType .= $refCtgryByGeneHsh{$refTrnscptTailPrxmtyID}.";";
				push @{$NGSJunctOnRefCommentHsh{$refTrnscptTailPrxmtyID}}, "prominentJunctAt5End$NGSJRefPrxmtyLimit" if (($prominentMinor eq "prominent") and ($refTrnscptSenseHitID eq "NA") and ($refTrnscptAntiSenseHitID eq "NA") and ($refIntronHitID eq "NA"));
			}
		}
		
		my $comment = "others";
		
		if (($prominentMinor eq "prominent") and ($refIntronHitID ne "NA")) {$comment = "prominent_refIntron_Hit";}
		elsif  (($prominentMinor eq "minor") and ($refIntronHitID ne "NA")) {$comment = "minor_refIntron_Hit";}
		elsif  (($prominentMinor eq "prominent") and ($refTrnscptSenseHitID ne "NA")) {$comment = "prominent_senseTrnscpt_Hit_$refTrnscptSenseHitType";}
		elsif  (($prominentMinor eq "minor") and ($refTrnscptSenseHitID ne "NA")) {$comment = "minor_senseTrnscpt_Hit_$refTrnscptSenseHitType";}
		elsif  (($prominentMinor eq "prominent") and ($withinPrxmty eq "yes")) {$comment = "prominent_proximity_Hit_$prxmtyType";}
		elsif  (($prominentMinor eq "minor") and ($withinPrxmty eq "yes")) {$comment = "minor_proximity_Hit_$prxmtyType";}
		elsif  (($prominentMinor eq "prominent") and ($refTrnscptAntiSenseHitID ne "NA")) {$comment = "prominent_antisenseTrnscpt_Hit_$refTrnscptAntiSenseHitType";}
		elsif  (($prominentMinor eq "minor") and ($refTrnscptAntiSenseHitID ne "NA")) {$comment = "minor_antisenseTrnscpt_Hit_$refTrnscptAntiSenseHitType";}
		elsif  ($prominentMinor eq "prominent") {$comment = "prominent_unknown";}
		elsif  ($prominentMinor eq "minor") {$comment = "minor_unknown";}
		
		$scenarioCountHsh{$comment}++;
		
		print NGSINTRONINFO join "", (join "\t", ($NGSIntronID, $NGSJunctStr, $comment, $prominentMinor, $refIntronHitID, $refIntronHitScenerio, $refTrnscptSenseHitID, $refTrnscptSenseHitType, $refTrnscptAntiSenseHitID, $refTrnscptAntiSenseHitType, $refTrnscptHeadPrxmtyID, $refTrnscptHeadPrxmtyDistance, $refTrnscptTailPrxmtyID, $refTrnscptTailPrxmtyDistance)), "\n";
	}
	close NGSINTRONINFO;
	
	open (STATS, ">>$outDir/stats.log.txt");
	print STATS ("\nStatistics on NGS junction scenarios\n");
	print STATS join "\t", ("percetage", "scenario", "count\n");

	print ("\nStatistics on NGS junction scenarios\n");
	print join "\t", ("percetage", "scenario", "count\n");
	
	foreach my $scenario (sort {$a cmp $b} keys %scenarioCountHsh) {
		my $count = $scenarioCountHsh{$scenario};
		my $tmpPct = sprintf '%.2f', ($count/$totalNGSIntronNum)*100;

		print STATS join "\t", ($tmpPct, $scenario, $count."\n");
		print join "\t", ($tmpPct, $scenario, $count."\n");
	}
	close STATS;

	print "\n\n";
	
	return \%NGSJunctOnRefCommentHsh;

}
########################################################################## readJunctBED
sub readJunctBED {

	my $junctBEDPath = $_[0];
	
	my (%junctRngHsh, %junctStrndHsh, %junctCtngHsh, %junctReadNumHsh, %junctScoreHsh, %junctIntronRngHsh);
	
	my $dupJunctNum = my $totalJunctNum = 0;
	
	open (INFILE, "$junctBEDPath");
	open (TMPLOG, ">$outDir/junctionInBedMoreThanOnce.txt");
	print "Reading $junctBEDPath\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		next if ($theLine =~ m/^track name/);
		$totalJunctNum++;
		print "$totalJunctNum junctions have been stored. $dupJunctNum junction appeared in the bed file twice.\r";
		my @theLineSplt = split (/\t/, $theLine);
		my $cntg = $theLineSplt[0];
		my $bedStart = $theLineSplt[1];
		my $bedEnd = $theLineSplt[2];
		my $strd = $theLineSplt[5];
		my @blkSizesSplt = split /,/, $theLineSplt[10];
		my $blk1Size = $blkSizesSplt[0];
		my $blk2Size = $blkSizesSplt[1];

		my $readNum = my $score = 0;

		if ($junctionBedType eq "HMMSplicer") {
			my @readNumAndScoreSplt = split /\|/, $theLineSplt[3];
			$readNum = substr ($readNumAndScoreSplt[0], index ($readNumAndScoreSplt[0], "=")+1);
			$score = substr ($readNumAndScoreSplt[1], index ($readNumAndScoreSplt[1], "=")+1);
		} elsif ($junctionBedType eq "tophat") {
			$readNum = $score = $theLineSplt[4];
		} else {
			die "unspecified junction bed type\n";
		}

		my $intronStart = $bedStart + $blk1Size + 1;
		my $intronEnd = $bedEnd - $blk2Size;

		my $junctStr = join ":", ($cntg, $intronStart, $intronEnd, "B"); #---assumed to be unique, the B refers to BED, analogous to "P" and "M" for prominent and minor

		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"start"} = $intronStart;
		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"end"} = $intronEnd;
		
		#---multiple $junctStr may exist in HMMSplicer as unique and duplicated individual junctions are collapsed seperately
		if (not exists $junctCtngHsh{$junctStr}) { #---most of the case
			$junctCtngHsh{$junctStr} = $cntg;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd;
			
			$junctStrndHsh{$junctStr} = $strd;
			$junctReadNumHsh{$junctStr} = $readNum;
			$junctScoreHsh{$junctStr} = $score;

		} else { #---appeared twice, extend the range 
			
			$dupJunctNum++;
			
			my $storedStart = ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"}; 
			my $storedEnd = ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
			my $storedScore = $junctScoreHsh{$junctStr};
			my $storedReadNum = $junctReadNumHsh{$junctStr};
			
			$junctReadNumHsh{$junctStr} = $readNum + $storedReadNum;
			$junctScoreHsh{$junctStr} = $score if ($score > $storedScore);;
			
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart if ($bedStart < $storedStart);
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd if ($bedEnd > $storedEnd);

			print TMPLOG "$junctStr appeared in the bed file twice. Collpasing the two BED lines. Score = $storedScore|$score, read=$storedReadNum|$readNum\n";
			
		}

	}	
	close INFILE;
	close TMPLOG;
	
	my $junctNum = keys %junctStrndHsh;
	
	print "\n\n";

	return (\%junctStrndHsh, \%junctRngHsh, \%junctCtngHsh, \%junctScoreHsh, \%junctReadNumHsh, \%junctIntronRngHsh);
	
}
########################################################################## scanJunctionInnerOutterBoundXSCoverage
sub scanJunctionInnerOutterBoundXSCoverage {#----scan for the coverage flanking the intron bounds and record the splicing ratio 

	my %junctIntronRngHsh = %{$_[0]};
	my $flankSize = $_[1];
	
	#---var to return
	my %junctIntronBoundCovHsh = ();

	#---read or calculate the intron bound coverage
	if ($intronBoundCovPath ne "no") {#---read
		
		print "Reading $intronBoundCovPath for intron bound coverage\n";
		
		open (INTRNBOUNDCOV, "$intronBoundCovPath")|| die "Can't read $intronBoundCovPath :$!\n";
		
		while (my $theLine = <INTRNBOUNDCOV>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			my $intronID = $theLineSplt[0];
			${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"out"} = $theLineSplt[1];
			${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"in"} = $theLineSplt[2];
			${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"out"} = $theLineSplt[3];
			${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"in"} = $theLineSplt[4];
		}
		close INTRNBOUNDCOV;
		
	} else {#calculate
		
		open (PILEUPFILE, "$pileupPath");
		open (INTRNBOUNDCOV, ">$outDir/refIntronBoundCov.$flankSize.txt");
	
		print "Reading $pileupPath.\n";
	
		my %tmpCovByPosHsh = ();
		
		#--get the first line
		my $theCurntLine;
		while ($theCurntLine = <PILEUPFILE>) {
			if ($theCurntLine !~ m/^[\@|\#]/) {#---ignore comments
				chomp $theCurntLine;
				last;
			}
		}
		
		my $cntgNum = 0;
		
		#--go through the rest of the file
		while (my $theNextLine = <PILEUPFILE>) {
			next if ($theNextLine =~ m/^[\@|\#]/);
			chomp $theNextLine;
		
			my @theNextLineSplt = split /\t/, $theNextLine;
			my $nextCtng = $theNextLineSplt[0];
	
			my @theCurntLineSplt = split /\t/, $theCurntLine;
			my $curntCtng = $theCurntLineSplt[0];
	
			my $pos = $theCurntLineSplt[1];
			my $plusCov = $theCurntLineSplt[3];
			my $minusCov = $theCurntLineSplt[4];
		
			#---store the cov if plus or minus strand > 0;
			if (($plusCov > 0) or ($minusCov > 0)) {
				${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
				${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
			}
			
			if (eof(PILEUPFILE)) {#---last line of the file
				my $pos = $theNextLineSplt[1];
				my $plusCov = $theNextLineSplt[3];
				my $minusCov = $theNextLineSplt[4];
	
				#---store the cov if plus or minus strand > 0;
				if (($plusCov > 0) or ($minusCov > 0)) {
					${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
					${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
				}
			}
			
			#---change contig or end of file
			if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
				$cntgNum++;
				print "Scanning intron bound coverage in contig $curntCtng. $cntgNum contig scanned.\r";
				
				foreach my $intronID (keys %{$junctIntronRngHsh{$curntCtng}}) {
					
					my $intronStart = ${${$junctIntronRngHsh{$curntCtng}}{$intronID}}{"start"};
					my $intronEnd = ${${$junctIntronRngHsh{$curntCtng}}{$intronID}}{"end"};
					
					my $startOutCovSum = my $startInCovSum = my $endOutCovSum = my $endInCovSum = 0;
					
					#---start inner flank region cov
					for my $pos (($intronStart)..($intronStart + $flankSize)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$startInCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$startInCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
						
					#---end inner flank region cov
					for my $pos (($intronEnd - $flankSize)..($intronEnd)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$endInCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$endInCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
					
					#---start outter flank region cov
					for my $pos (($intronStart - $flankSize)..($intronStart)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$startOutCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$startOutCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
						
					#---end outter flank region cov
					for my $pos (($intronEnd)..($intronEnd + $flankSize)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$endOutCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$endOutCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}

					${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"out"} = sprintf "%.02f", $startOutCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"in"} = sprintf "%.02f", $startInCovSum/$flankSize;
					
					${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"out"} = sprintf "%.02f", $endOutCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"in"} = sprintf "%.02f", $endInCovSum/$flankSize;
	
					print INTRNBOUNDCOV $intronID."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"out"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$intronID}}{"start"}}{"in"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"out"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$intronID}}{"end"}}{"in"}."\n";
	
				} #---end of foreach my $intronID (keys %{$junctIntronRngHsh{$curntCtng}})
				
				%tmpCovByPosHsh = (); #---empty the hash
	
			} #---end of if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
			
			$theCurntLine = $theNextLine;
			
		}#---end of while (my $theNextLine = <PILEUPFILE>) {
	
		close PILEUPFILE;
		close INTRNBOUNDCOV;
	}
	
	print "\n";
	
	return (\%junctIntronBoundCovHsh);
	
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];

	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		open (OUTDIRCMDLOG, ">$outDir/run.cmd.log.txt"); #---make the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print OUTDIRCMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close OUTDIRCMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}

}
########################################################################## readMultiFasta
sub readMultiFasta {

	my $refFastaPath = $_[0];
	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	print "Reading $refFastaPath into a hash.\n";
	open (INFILE, $refFastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq = $seq.$nextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;
	
	return (\%fastaHsh);
}
########################################################################## getMotifRng
sub getMotifRng {
	
	my %refSeqHsh = %{$_[0]};
	
	my (%ATRichRngXSHsh, %NRegionRngXSHsh, %cntgEndRngXSHsh);
	
	my $winSize = 50;
	my $stepSize = 10;
	my $ATRichPctCutoff = 90;
	my $winHalf = int ($winSize/2);
	my $NPctCutoff = 100;
	
	printProgressScale("\nScanning for N regions and AT rich regions", 50, 10);
	my $progressCntgNum = 0;

	my $totalCntgNum = keys %refSeqHsh;

	foreach my $cntg (keys %refSeqHsh) {
		my $cntgSeq = $refSeqHsh{$cntg};
		my $cntgLen = length $cntgSeq;

		#---update the on-screen progress
		$progressCntgNum++;
		updateProgressBar("$cntg", $progressCntgNum, $totalCntgNum, 50, 10);

		#----set the window mid-points
		my $winMidMin = $winHalf;
		my $winMidMax = $cntgLen - $winHalf - 1;
		my @winMidPointAry;
		
		for (my $i=$winMidMin; $i <= $winMidMax; $i = $i+$stepSize) {
			push @winMidPointAry, $i;
		}
		
		#---slide the window!
		foreach my $winMidPoint (@winMidPointAry) {

			${${$cntgEndRngXSHsh{$cntg}}{$cntg."_head"}}{"start"} = 0;
			${${$cntgEndRngXSHsh{$cntg}}{$cntg."_head"}}{"end"} = 10;

			${${$cntgEndRngXSHsh{$cntg}}{$cntg."_tail"}}{"start"} = $cntgLen - 10;
			${${$cntgEndRngXSHsh{$cntg}}{$cntg."_tail"}}{"end"} = $cntgLen;
		
			my $winStart = $winMidPoint - $winHalf;  
			my $winSeq = substr $cntgSeq, $winStart, $winSize;
			my $A=($winSeq=~tr/A//);
			my $T=($winSeq=~tr/T//);
			my $G=($winSeq=~tr/G//);
			my $C=($winSeq=~tr/C//);
			my $N=($winSeq=~tr/N//);
			my $motifID = $cntg."_".$winStart;

			#----N window
			if ((100*$N/$winSize) >= $NPctCutoff) {
				${${$NRegionRngXSHsh{$cntg}}{$motifID}}{"start"} = $winStart;
				${${$NRegionRngXSHsh{$cntg}}{$motifID}}{"end"} = $winStart + $winSize;
			}
			
			#----AT rich window
			if ((100*($A+$T)/$winSize) >= $ATRichPctCutoff) {
				${${$ATRichRngXSHsh{$cntg}}{$motifID}}{"start"} = $winStart;
				${${$ATRichRngXSHsh{$cntg}}{$motifID}}{"end"} = $winStart + $winSize;
			}
		}	
	}
	
	print "\n\n";
	
	return (\%ATRichRngXSHsh, \%NRegionRngXSHsh, \%cntgEndRngXSHsh);
	
}
########################################################################## assignDummyStrand
sub assignDummyStrand {
	
	my %rngHsh = %{$_[0]};
	
	my %strndHsh;
	
	foreach my $cntg (keys %rngHsh) {
		foreach my $ID (keys %{$rngHsh{$cntg}}) {
			$strndHsh{$ID} = ".";
		} 
	}
	
	return \%strndHsh;
}
########################################################################## summerizeRefTranscriptBasedNGSHitInfo
sub summerizeRefTranscriptBasedNGSHitInfo {
	
	#---$refBasedNGSHitInfoHsh_ref = summerizeRefTranscriptBasedNGSHitInfo($allScenarioByRefIntronIDHsh_ref, $refBasedNGSHitInfoHsh_ref, $refIntronRngByGeneHsh_ref, $NRegionRefHitHsh_ref, $NRegionRefPrxmtyHsh_ref, $cntgEndRefHitHsh_ref, $cntgEndRefPrxmtyHsh_ref);
	my %allScenarioByRefIntronIDHsh = %{$_[0]};
	my %refBasedNGSHitInfoHsh = %{$_[1]};
	my %refIntronRngByGeneHsh = %{$_[2]};
	my %NRegionRefHitHsh = %{$_[3]};
	my %NRegionRefPrxmtyHsh = %{$_[4]};
	my %cntgEndRefHitHsh = %{$_[5]};
	my %cntgEndRefPrxmtyHsh = %{$_[6]};
	my %NGSJunctOnRefCommentHsh = %{$_[7]};
	my $cutoffCovPct = $_[8];
	
	open (GENEBAESDINFO, ">$outDir/geneBasedNGSHitInfo.tsv");

	print GENEBAESDINFO join "", ((join "\t", ("geneID", "inherentCriterion", "coverageCriterion", "junctionCriterion", "coveredPos", "refLength", "NGSTransfragHit", "NGSTransfragCovPct", "end3TransfragDiff", "end5TransfragDiff", "totalIntronNum", "confirmedIntron", "modifiedIntron", "unknownIntron", "intronScenerio", "end5CntgEndPrxmty", "end3CntgEndPrxmty", "end5NRegPrxmty", "end3NRegPrxmty", "NRegionHit", "cntgEndHit",  "newJunctionComment", "comment")), "\n");

	#---go through each gene
	foreach my $geneID (keys %refBasedNGSHitInfoHsh) {

		my $inherentCriterion = my $coverageCriterion = my $junctionCriterion = 'passed';
		my @commentAry = ();
		
		#----get the potentially new junction comments
		if (exists $NGSJunctOnRefCommentHsh{$geneID}) {
			$junctionCriterion = "revise";
			push @commentAry, "potentialNewJunction";
		} else {
			push @{$NGSJunctOnRefCommentHsh{$geneID}}, "no";
		}

		my $newJunctionComment = join ";", @{$NGSJunctOnRefCommentHsh{$geneID}};
		
		#---summarize the intron confirmation
		my $totalIntronNum = my $confirmedIntron = my $modifiedIntron = my $unknownIntron = 0;
		my $intronScenerio = "na";
		if (exists $refIntronRngByGeneHsh{$geneID}) {
			my @intronScenerioAry;
			foreach my $intronNum (sort {$a <=> $b} keys %{$refIntronRngByGeneHsh{$geneID}}) {
				my $intronID = $geneID.":".$intronNum;
				$totalIntronNum++;
				if ($allScenarioByRefIntronIDHsh{$intronID} =~ m/confirm/) {
					$confirmedIntron++;

				} elsif ($allScenarioByRefIntronIDHsh{$intronID} =~ m/modify/) {
					$modifiedIntron++;
					$junctionCriterion = "revise";
					push @commentAry, 'modifyJunction';
					
				} elsif ($allScenarioByRefIntronIDHsh{$intronID} =~ m/unknown/) {
					$unknownIntron++;
					$junctionCriterion = "failed";
					push @commentAry, 'missedJunction';
				
				}
				push @intronScenerioAry, $allScenarioByRefIntronIDHsh{$intronID};
			}
			$intronScenerio = join "//", @intronScenerioAry;
		}
		
		${$refBasedNGSHitInfoHsh{$geneID}}{'totalIntronNum'} = $totalIntronNum;
		${$refBasedNGSHitInfoHsh{$geneID}}{'confirmedIntron'} = $confirmedIntron;
		${$refBasedNGSHitInfoHsh{$geneID}}{'modifiedIntron'} = $modifiedIntron;
		${$refBasedNGSHitInfoHsh{$geneID}}{'unknownIntron'} = $unknownIntron;
		${$refBasedNGSHitInfoHsh{$geneID}}{'intronScenerio'} = $intronScenerio;
		
		#---NRegion and CntgEnd info
		my $NRegionHit = my $cntgEndHit = "no";

		$NRegionHit = 'yes' if exists $NRegionRefHitHsh{$geneID};
		$cntgEndHit = 'yes' if exists $cntgEndRefHitHsh{$geneID};
		
		my $end5CntgEndPrxmty = ${${$cntgEndRefPrxmtyHsh{$geneID}}{"H"}}[0];
		my $end3CntgEndPrxmty = ${${$cntgEndRefPrxmtyHsh{$geneID}}{"T"}}[0];
		my $end5NRegPrxmty = ${${$NRegionRefPrxmtyHsh{$geneID}}{"H"}}[0];
		my $end3NRegPrxmty = ${${$NRegionRefPrxmtyHsh{$geneID}}{"T"}}[0];

		if (($NRegionHit eq 'yes') 
		or ($cntgEndHit eq 'yes') 
		or (($end5NRegPrxmty <= 100) and ($end5NRegPrxmty > -999)) 
		or (($end3NRegPrxmty <= 100) and ($end3NRegPrxmty > -999)) 
		or (($end5CntgEndPrxmty <= 100) and ($end5CntgEndPrxmty > -999)) 
		or (($end3CntgEndPrxmty <= 100) and ($end3CntgEndPrxmty > -999))){
			$inherentCriterion = 'failed';
			push @commentAry, 'ambiguousReg';
		}

		${$refBasedNGSHitInfoHsh{$geneID}}{'5EndCntgEndPrxmty'} = $end5CntgEndPrxmty;
		${$refBasedNGSHitInfoHsh{$geneID}}{'3EndCntgEndPrxmty'} = $end3CntgEndPrxmty;
		${$refBasedNGSHitInfoHsh{$geneID}}{'5EndNRegPrxmty'} = $end5NRegPrxmty;
		${$refBasedNGSHitInfoHsh{$geneID}}{'3EndNRegPrxmty'} = $end3NRegPrxmty;
		${$refBasedNGSHitInfoHsh{$geneID}}{'NRegionHit'} = $NRegionHit;
		${$refBasedNGSHitInfoHsh{$geneID}}{'cntgEndHit'} = $cntgEndHit;

		#----get the info and print the info
		my $NGSTransfragHit = ${$refBasedNGSHitInfoHsh{$geneID}}{'NGSTransfragHit'};
		my $NGSTransfragCovPct = ${$refBasedNGSHitInfoHsh{$geneID}}{'NGSTransfragCovPct'};
		my $end3TransfragDiff = ${$refBasedNGSHitInfoHsh{$geneID}}{'3EndTransfragDiff'};
		my $end5TransfragDiff = ${$refBasedNGSHitInfoHsh{$geneID}}{'5EndTransfragDiff'};
		my $refLength = ${$refBasedNGSHitInfoHsh{$geneID}}{'refLength'};
		my $coveredPos = ${$refBasedNGSHitInfoHsh{$geneID}}{'coveredPos'};

=pod	
		#----inactivated in v0.9 to relax the "passed" criteria
		
		if ($NGSTransfragHit > 0) {
			if ((($end3TransfragDiff <= 300) and ($end3TransfragDiff > 100))) { 
				$coverageCriterion = "revise";
				push @commentAry, 'extendedTransfragCovAt3End';
			}
			if ((($end5TransfragDiff <= 300) and ($end5TransfragDiff > 100))) { 
				$coverageCriterion = "revise";
				push @commentAry, 'extendedTransfragCovAt5End';
			}
		}
=cut			
		
		if ($NGSTransfragCovPct <= $cutoffCovPct) {
			$coverageCriterion = "failed";
			push @commentAry, "covLessThan$cutoffCovPct"."%";
		}
		
		push @commentAry, "passed" if (($inherentCriterion eq "passed") and ($coverageCriterion eq "passed") and ($junctionCriterion eq "passed"));
		my $comment = join ";", @commentAry;

		print GENEBAESDINFO join "", ((join "\t", ($geneID, $inherentCriterion, $coverageCriterion, $junctionCriterion, $coveredPos, $refLength, $NGSTransfragHit, $NGSTransfragCovPct, $end3TransfragDiff, $end5TransfragDiff, $totalIntronNum, $confirmedIntron, $modifiedIntron, $unknownIntron, $intronScenerio, $end5CntgEndPrxmty, $end3CntgEndPrxmty, $end5NRegPrxmty, $end3NRegPrxmty, $NRegionHit, $cntgEndHit, $newJunctionComment, $comment)), "\n");

	}
	close GENEBAESDINFO;
	
	return \%refBasedNGSHitInfoHsh;

}
