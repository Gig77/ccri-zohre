use warnings FATAL => qw( all );
use strict;

use Generic;
use Log::Log4perl qw(:easy);
use Tabix;
use Vcf;
use Data::Dumper;
use Getopt::Long;
use Carp;

my ($vcf_out, $header, $rejected_variants_file, $sample);
my ($rmsk_file, $simplerepeat_file, $segdup_file);
GetOptions
(
	"sample=s" => \$sample,  # sample ID
	"vcf-out=s" => \$vcf_out,  # filtered VCF output file
	"header" => \$header,  # if set, write header line to output
	"rmsk-file=s" => \$rmsk_file, # TABIX indexed UCSC table rmsk
	"simpleRepeat-file=s" => \$simplerepeat_file, # TABIX indexed UCSC table rmsk
	"segdup-file=s" => \$segdup_file # TABIX indexed UCSC table genomicSuperDups
);

# TABLE: filtered-variants
if ($header)
{
	print "sample\t";		
	print "var_type\t";
	print "somatic_status\t";
	print "status\t";
	print "rejected_because\t";
	print "chr\t";
	print "pos\t";
	print "dbSNP\t";
	print "ref\t";
	print "alt\t";
	print "gene\t";
	print "add_genes\t";
	print "impact\t";
	print "effect\t";
	print "exonIntronNum\t";
	print "dp_rem_tot\t";
	print "dp_rem_ref\t";
	print "dp_rem_var\t";
	print "freq_rem\t";
	print "dp_tum_tot\t";
	print "dp_tum_ref\t";
	print "dp_tum_var\t";
	print "freq_tum\t";
	print "both_strands_tum_var\t";
	print "somatic_pvalue\t";
	print "aa_change\t";
	print "SnpEffANN\t";
	print "repeat\t";
	print "tsegdup\n";
#	exit;
}

my $debug = 1;

my $vcf_file = $ARGV[0] or croak "ERROR: VCF file not specified\n";

croak "ERROR: --sample not specified" if (!$sample);
croak "ERROR: --rmsk-file not specified" if (!$rmsk_file);
croak "ERROR: --simpleRepeat-file not specified" if (!$simplerepeat_file);
croak "ERROR: --segdup-file not specified" if (!$segdup_file);


my $rem_sample = "NORMAL"; 
my $tum_sample = "TUMOR"; 

# get canonical UCSC transcripts
my $knownCanonicalFile = "/data_synology/max/zohre/ucsc/mm10.knownCanonical.txt";
my %canonicalIds;
open(G, $knownCanonicalFile) or die "$knownCanonicalFile";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($chrom, $chromStart, $chromEnd, $clusterId, $transcript, $protein) = split(/\t/);
	$canonicalIds{$transcript} = 1;
}
close(G);
INFO(scalar(keys(%canonicalIds))." canonical genes read from file $knownCanonicalFile");

# map ensembl transcript IDs to canonical UCSC IDs
my $knownToEnsemblFile = "/data_synology/max/zohre/ucsc/mm10.knownToEnsembl.txt";
my %isCanonical;
open(G, $knownToEnsemblFile) or die "could not open file $knownToEnsemblFile";
<G>; # skip header
while(<G>)
{
	chomp;
	my ($name, $value) = split(/\t/);
	$isCanonical{$value} = $name if (exists $canonicalIds{$name}) 
}
close(G);
INFO(scalar(keys(%isCanonical))." mappings read from file $knownToEnsemblFile");

my $rmsk = Tabix->new(-data => $rmsk_file);
my $simpleRepeat = Tabix->new(-data => $simplerepeat_file);
my $segdup = Tabix->new(-data => $segdup_file);

$| = 1; # turn on autoflush

INFO("Processing file $vcf_file...");

my $vcf = Vcf->new(file => "$vcf_file");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

if ($vcf_out) 
{
	my $cmd = "grep -P '^#' $vcf_file > $vcf_out";
	system($cmd) == 0 or die "ERROR: grep vcf header failed: $cmd\n";
	open(VCFOUT,">>$vcf_out") or die "ERROR: could not write to file $vcf_out\n";
}

# sanity check
die "ERROR: Sample name $rem_sample not found!\n" if ($rem_sample ne $samples[0] and $rem_sample ne $samples[1]);
die "ERROR: Sample name $tum_sample not found!\n" if ($tum_sample ne $samples[0] and $tum_sample ne $samples[1]);

my ($tot_var, $filtered_alt, $filtered_germ) = (0, 0, 0);
my ($numrep, $numsegdup) = (0, 0);
my %qual_num;

my %somatic_stati = 
(
	0 => 'reference',
	1 => 'germline',
	2 => 'somatic',
	3 => 'LOH',
	5 => 'unknown'
);

while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);

	my $status = $x->{FILTER}->[0];		

	$tot_var ++;
	$qual_num{$status} = $qual_num{$status} ? $qual_num{$status} + 1 : 1;
	
	##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)">
	my $somatic_status = $x->{INFO}{'SS'};
	if ($somatic_status == 0 or $somatic_status == 1) # reference or germline?
	{
		$filtered_germ ++;
		next;	
	}
	
	##INFO=<ID=SPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls">
	my $somatic_pvalue = $x->{INFO}{SPV};
	
	my $gt_rem = $x->{gtypes}{$rem_sample}{GT};
	die "ERROR: Could not determine genotype of sample $rem_sample in file $vcf_file\n" if (!defined $gt_rem or $gt_rem eq "");

	if (@{$x->{ALT}} != 1) # more than one alternative allele?
	{
		$filtered_alt ++;
		next;
	}		
	
	my ($dp_tum, $dp_rem, $freq_tum, $freq_rem, $ad_tum_ref, $ad_tum_alt, $ad_rem_ref, $ad_rem_alt, $ref_fwd, $ref_rev, $var_fwd, $var_rev);
	
	my $var_type = length($x->{REF}) == length($x->{ALT}->[0]) ? "snp" : "indel";
	
	##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
	##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
	$ad_tum_ref = $x->{gtypes}{$tum_sample}{RD};
	$ad_tum_alt = $x->{gtypes}{$tum_sample}{AD};
	$ad_rem_ref = $x->{gtypes}{$rem_sample}{RD};
	$ad_rem_alt = $x->{gtypes}{$rem_sample}{AD};

	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
	$dp_tum = $x->{gtypes}{$tum_sample}{DP};
	$dp_rem = $x->{gtypes}{$rem_sample}{DP};

	##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
	$freq_tum = $x->{gtypes}{$tum_sample}{FREQ};
	$freq_tum =~ s/\%//;
	$freq_rem = $x->{gtypes}{$rem_sample}{FREQ};		
	$freq_rem =~ s/\%//;
		
	##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">
	($ref_fwd, $ref_rev, $var_fwd, $var_rev) = split(",", $x->{gtypes}{$tum_sample}{DP4});

	my (@repeats, @dups);
	my ($chr, $pos) = ($x->{CHROM}, $x->{POS});

	# ----- annotate overlapping repeat regions
	{
		my $iter = $rmsk->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $rmsk->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[10]:$s[11]:$s[12]");
			}
		}		
	}

	{
		my $iter = $simpleRepeat->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $simpleRepeat->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@repeats, "$s[16]($s[6])");
			}
		}		
	}
	$numrep ++ if (@repeats > 0);

	# ----- annotate segmental duplications
	{
		my $iter = $segdup->query($chr, $pos-1, $pos+1);
		if ($iter and $iter->{_})
		{
			while (my $line = $segdup->read($iter)) 
			{
				my @s = split("\t", $line);
				push(@dups, "$s[4]:".sprintf("%.2f", $s[26]));
			}		
		}		
		$numsegdup ++ if (@dups > 0);
	}	
	
	my $reject = 0;
	my @reject_because;
	my $both_strands = ($var_fwd and $var_rev);

	if ($somatic_status == 2 and $ad_rem_alt > 0) { $reject = 1; push(@reject_because, "present normal"); };
	if ($somatic_status == 2 and !$both_strands) { $reject = 1; push(@reject_because, "single strand"); };
	if ($somatic_status == 5) { $reject = 1; push(@reject_because, "unknown somatic status"); };
	if ($somatic_pvalue > 0.1) { $reject = 1; push(@reject_because, "somatic p-value"); };
	if (@repeats > 0) { $reject = 1; push(@reject_because, "repetitive region"); };
	if (@dups > 0) { $reject = 1; push(@reject_because, "segmental duplication"); };
	if ($x->{ID} and $x->{ID} ne ".")  { $reject = 1; push(@reject_because, "dbSNP"); };  
	
	$line =~ s/^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t)[^\t]+/$1REJECT/ if ($reject);
	
	print VCFOUT "$line" if ($vcf_out);
	
	print "$sample\t";		
	print "$var_type\t";
	print $somatic_stati{$somatic_status}, "\t";
	print $reject ? "REJECT\t" : "$status\t";
	print join(";", @reject_because), "\t";
	print $x->{CHROM},"\t";
	print $x->{POS},"\t";
	print $x->{ID},"\t";
	print $x->{REF},"\t";
	print $x->{ALT}->[0],"\t";

	my ($gene, $add_genes, $impact, $effect, $affected_exon, $aa_change) = get_impact($x->{INFO}{ANN});
	print "$gene\t";
	print "$add_genes\t";
	print "$impact\t";
	print "$effect\t";
	print "$affected_exon\t";
	print "$dp_rem\t";
	print "$ad_rem_ref\t";
	print "$ad_rem_alt\t";
	print "$freq_rem\t";
	print "$dp_tum\t";
	print "$ad_tum_ref\t";
	print "$ad_tum_alt\t";
	print "$freq_tum\t";
	print $both_strands ? "yes" : "no", "\t";
	print defined $somatic_pvalue ? $somatic_pvalue : "", "\t";
	print "$aa_change\t";
	print "ANN=",$x->{INFO}{ANN},"\t";
	print join(',', @repeats), "\t", join(',', @dups);
	print "\n";		
}
$vcf->close();
close(VCFOUT) if ($vcf_out);
	
if ($debug)
{
	INFO("  Total number of variants: $tot_var");
	INFO("  Variants by quality:");
	foreach my $k (keys(%qual_num))
	{
		INFO("    $k: ", $qual_num{$k});
	}
	INFO("  Excluded germline variants: $filtered_germ");
	INFO("  Excluded due to missing alternative allele: $filtered_alt");
	INFO("  $numrep variants annotated with overlapping repeat.");
	INFO("  $numsegdup variants annotated with overlapping segmental duplication.");
}

# ------------------------------------------

sub get_impact
{
	my $effs = shift or die "ERROR: effect not specified";

	# determine all genes impacted by variants
	my (%genes_by_impact, %all_genes, $combined_effect, $combined_impact, %affected_exons, %aa_changes);
	foreach my $eff (split(",", $effs))
	{
		my ($allele, $annotation, $annotation_impact, $gene_name, $gene_id, $feature_type, $feature_id, 
		    $transcript_biotype, $rank, $HGVSc, $HGVSp, $cDNA_pos, $cDNA_length, $CDS_pos_CDS_length, 
		    $AA_pos_AA_length, $distance, $info) = split('\|', $eff)
				or die "ERROR: could not parse SNP effect: $eff\n";

		$feature_id =~ s/\.\d+$//; # remove version number from accession
		$feature_id =~ s/\.\d+$//; 

		if ($gene_name)
		{
			$genes_by_impact{$annotation_impact}{$gene_name}{$feature_id} = $annotation;
			$all_genes{$gene_name} = 1;
			
			if ($rank)
			{
				$affected_exons{$gene_name}{$rank}{$feature_id} = 1;
				if ($isCanonical{$feature_id})
				{
					$affected_exons{$gene_name}{'canonical'}{$rank}{$feature_id} = 1;
				}
			}
		}		 

		$aa_changes{$HGVSp} = 1 if ($HGVSp);
		$combined_impact = $annotation_impact;		
		$combined_effect = $annotation;
	}
	
	# if multiple genes are affected, preferentially chose gene with the predicted higher impact
	if ($genes_by_impact{HIGH})
	{
		$combined_impact = "HIGH";
	}
	elsif ($genes_by_impact{MODERATE})
	{
		$combined_impact = "MODERATE";
	}
	elsif ($genes_by_impact{LOW})
	{
		$combined_impact = "LOW";
	}
	elsif ($genes_by_impact{MODIFIER})
	{
		$combined_impact = "MODIFIER";
	}
	
	# if multiple genes in highest category, arbitrarily choose first of alphabetical order
	my ($gene, $add_genes) = ("", "");
	if (keys(%all_genes) > 0)
	{
		my @sorted_genes = sort keys(%{$genes_by_impact{$combined_impact}});
		$gene = $sorted_genes[0]; # first choice is first in alphabetically sorted list
		delete $all_genes{$gene};
		$add_genes = join(",", keys(%all_genes));
	
		foreach my $t (keys(%{$genes_by_impact{$combined_impact}{$gene}}))
		{
			if ($isCanonical{$t}) {
				$combined_effect = $genes_by_impact{$combined_impact}{$gene}{$t};
				last;		
			}
		}
	}

	my @aff_exons;
	foreach my $g (keys(%affected_exons))
	{
		if (exists $affected_exons{$g}{'canonical'}) # is there a known canonical transcript for this gene? --> only use these
		{
			foreach my $e (keys(%{$affected_exons{$g}{'canonical'}}))
			{
				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{'canonical'}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
		}
		else
		{
			foreach my $e (keys(%{$affected_exons{$g}}))
			{
				next if ($e eq 'canonical');

				my @transcripts;
				foreach my $t (keys(%{$affected_exons{$g}{$e}}))
				{
					push(@transcripts, "$g:$t");
				}
				push(@aff_exons, "$e (".join(";", @transcripts).")");
			}
			
		}
	}

	return ($gene, $add_genes, $combined_impact, $combined_effect, 
			@aff_exons > 0 ? join(",", @aff_exons) : "", join(";", keys(%aa_changes)));
}
