use warnings FATAL => qw( all );
use strict;

use Vcf;
use Getopt::Long;

my $min_dp = 20;
my $max_dp = 300;
my $max_pl_primary_gt = 0;      # lower is more stringent (0-255)
my $min_pl_secondary_gt = 100;  # higher is more stringent (0-255)
my $only_dbsnp = 0;
my $max_vdb = 0.001;              # Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)
my $max_rpb = 0.001;              # Mann-Whitney U test of Read Position Bias (bigger is better)
my $max_mqb = 0.001;              # Mann-Whitney U test of Mapping Quality Bias (bigger is better)
my $max_bqb = 0.001;              # Mann-Whitney U test of Base Quality Bias (bigger is better)

my ($vcf_file, $sample, $not_homozygous_in_sample);
GetOptions
(
	"vcf-file=s" => \$vcf_file,                                 # VCF input file
	"sample=s" => \$sample,                                     # sample name
	"not-homozygous-in-sample=s" => \$not_homozygous_in_sample  # do not output BAFs for variants that are homozygous in this sample (optional)
);

die "ERROR: --vcf-file not specified" if (!$vcf_file);

my $vcf = Vcf->new(file => "$vcf_file");
$vcf->parse_header();
my (@samples) = $vcf->get_samples();

die "ERROR: Sample $sample not present in VCF file $vcf_file\n" if (!grep( /^$sample$/, @samples ));
die "ERROR: Sample $not_homozygous_in_sample not present in VCF file $vcf_file\n" if ($not_homozygous_in_sample and !grep( /^$not_homozygous_in_sample$/, @samples ));

#$| = 1; # turn on autoflush

my $total_variants = 0;
my $kept_variants = 0;
my $filtered_dbsnp = 0;
my $filtered_mas = 0;
my $filtered_vdb = 0;
my $filtered_indel = 0;
my $filtered_rpb = 0;
my $filtered_mqb = 0;
my $filtered_bqb = 0;
my $filtered_dp_low = 0;
my $filtered_dp_high = 0;
my $filtered_pl_prim = 0;
my $filtered_pl_sec = 0;
my $filtered_hom = 0;

while (my $line = $vcf->next_line())
{
	my $x = $vcf->next_data_hash($line);
	
	$total_variants ++;

	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter) =  ($x->{CHROM}, $x->{POS}, $x->{ID}, $x->{REF}, $x->{ALT}, $x->{QUAL}, $x->{FILTER});
	
	# skip non-dbSNP sites
	if ($only_dbsnp and $id eq ".") {
		$filtered_dbsnp ++;
		next;
	}
	
	# skip multi-allelic sites
	if (@{$alt} != 1) {
		$filtered_mas ++;
		next;
	}
	$alt = $alt->[0];

	# skip indels
	if (length($ref) != length($alt)) {
		$filtered_indel ++;
		next;
	}

	# variant distance bias?
	my $vdb = $x->{INFO}{VDB};
	if (defined $vdb and $vdb < $max_vdb) {
		$filtered_vdb ++;
		next;
	}

	# read position bias?
	my $rpb = $x->{INFO}{RPB};
	if (defined $rpb and $rpb < $max_rpb) {
		$filtered_rpb ++;
		next;
	}
	
	# mapping quality bias?
	my $mqb = $x->{INFO}{MQB};
	if (defined $mqb and $mqb < $max_mqb) {
		$filtered_mqb ++;
		next;
	}

	# base quality bias?
	my $bqb = $x->{INFO}{BQB};
	if (defined $bqb and $bqb < $max_bqb) {
		$filtered_bqb ++;
		next;
	}
	
	# read depth too low/high?
	my $dp = $x->{gtypes}{$sample}{DP};
	die "ERROR: Missing genotype field 'DP' for sample $sample:\n$line" if (!defined $dp);
	if ($dp < $min_dp) {
		$filtered_dp_low ++;
		next;
	}
	if ($dp > $max_dp) {
		$filtered_dp_high ++;
		next;
	}

	# sufficient genotype quality?
	my $pl = $x->{gtypes}{$sample}{PL};
	die "ERROR: Missing genotype field 'PL' for sample $sample:\n$line" if (!defined $pl);
	my @pl_sorted = sort { $a <=> $b } split(",", $pl); 
	if ($pl_sorted[0] > $max_pl_primary_gt) {
		$filtered_pl_prim ++;
		next;
	}
	if ($pl_sorted[1] < $min_pl_secondary_gt) {
		$filtered_pl_sec ++;
		next;
	}
	
	if ($not_homozygous_in_sample) {
		my $gt = $x->{gtypes}{$not_homozygous_in_sample}{GT};
		die "ERROR: Missing genotype field 'GT' for sample $not_homozygous_in_sample:\n$line" if (!defined $gt);
		if ($gt eq "0/0" or $gt eq "1/1") {
			$filtered_hom ++;
			next;
		}
	}
	
	my $ad_str = $x->{gtypes}{$sample}{AD};
	die "ERROR: Missing genotype field 'AD' for sample $sample:\n$line" if (!defined $ad_str);

	my @ad = split(",", $ad_str);
	my $baf = $ad[1] / ($ad[0] + $ad[1]);
	
	print "$chrom\t".($pos-1)."\t$pos\t$ref:$alt\t$baf\n";
	
	$kept_variants ++;
}

print STDERR "$total_variants total variants for sample $sample.\n";
print STDERR "Filtered $filtered_dbsnp non-dbSNP variants".($only_dbsnp ? "" : " (filter disabled)")."\n";
print STDERR "Filtered $filtered_mas multi-allelic variants.\n";
print STDERR "Filtered $filtered_indel indels.\n";
print STDERR "Filtered $filtered_vdb variants with variant distance bias < $max_vdb.\n";
print STDERR "Filtered $filtered_rpb variants with read position bias < $max_rpb.\n";
print STDERR "Filtered $filtered_mqb variants with mapping quality bias < $max_mqb.\n";
print STDERR "Filtered $filtered_bqb variants with base quality bias < $max_bqb.\n";
print STDERR "Filtered $filtered_dp_low variants with read depth < $min_dp.\n";
print STDERR "Filtered $filtered_dp_high variants with read depth > $max_dp.\n";
print STDERR "Filtered $filtered_pl_prim variants with primary genotype quality > $max_pl_primary_gt.\n";
print STDERR "Filtered $filtered_pl_sec variants with secondary genotype quality < $min_pl_secondary_gt.\n";
print STDERR "Filtered $filtered_hom because homozygous in sample $not_homozygous_in_sample.\n" if ($not_homozygous_in_sample);
print STDERR "Kept $kept_variants variants for sample $sample.\n";