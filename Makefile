export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete intermediate files
.SECONDEXPANSION: # allow functions in dependency list
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=/mnt/projects/zohre
FASTQ=/data_synology_nfs/max/zohre/exome_maus/fastq
REFGENOME=/data_synology/max/zohre/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
CHRSIZES=/data_synology/max/zohre/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.chrsizes
COVERBED=/data_synology/max/zohre/agilent/S0276129_Covered.mm10.bed
DOCKER=docker run -i --rm --net=host -e DOCKER_UID=$$(id -u) -e DOCKER_UNAME=$$(id -un) -e DOCKER_GID=$$(id -g) -e DOCKER_GNAME=$$(id -gn) -e DOCKER_HOME=$$HOME \
       -v /home:/home \
       -v /data_synology:/data_synology \
       -v /mnt/biohazard/home/cf/zohre/results:$(PROJECT_HOME)/results \
       -v /data_synology/max/zohre:$(PROJECT_HOME)/data:ro \
       -v /data_synology/christian/generic/data/current:/mnt/projects/generic/data:ro \
       -w $$(pwd)
BWA=flock -x .lock /data_synology/software/bwa-0.7.12/bwa
SAMTOOLS=/data_synology/software/samtools-0.1.19/samtools
SAMTOOLS131=/data_synology/software/samtools-1.3.1/samtools
BEDTOOLS=/data_synology/software/bedtools-2.17.0/bin/bedtools
PICARD=$(DOCKER) biowaste:5000/ccri/picard-2.2.2 java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /usr/picard/picard.jar
GATK=java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar /data_synology/software/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
VARSCAN=java -jar /data_synology/software/varscan-2.4.2/VarScan.v2.4.2.jar
SNPEFF=java -Xmx4g -jar /data_synology/software/snpeff-4.2/snpEff.jar -v GRCm38.82
SNPSIFT=java -Xmx4g -jar /data_synology/software/snpeff-4.2/SnpSift.jar
DBSNP=/data_synology/max/zohre/snps/ftp.ncbi.nih.gov/snp/organisms/mouse_10090/VCF/mouse.dbSnp146.vcf.gz
VCFCONCAT=/data_synology/software/vcftools_0.1.10/bin/vcf-concat
VCFSORT=/data_synology/software/vcftools_0.1.10/bin/vcf-sort
BCFTOOLS=/data_synology/software/bcftools-1.3.1/bcftools
VCFTOOLS=/data_synology/software/vcftools-0.1.14/bin/vcftools
VCFUTILS=/data_synology/software/bcftools-1.3.1/vcfutils.pl
SEGMENTCOV=/mnt/projects/hdall/scripts/cnv/get-segment-coverage.pl
PLOTSEGCOV=/mnt/projects/zohre/scripts/plot-segment-coverage.R

PAIRS=test 11291 11682 11689 11689_2 11232 11221Up 11746_1 9193_1 11746_2 9193_2
SAMPLES_TUMOR=$(addsuffix T,$(PAIRS)) 
SAMPLES_NORMAL=$(addsuffix K,$(subst Up,,$(subst _2,,$(subst _1,,$(PAIRS))))) 
SAMPLES=$(SAMPLES_TUMOR) $(SAMPLES_NORMAL)

all: bwa picard gatk varscan snpeff cna filtered-variants sex/allsamples.readcount.sexchr.txt

#include /mnt/projects/generic/scripts/rna-seq/fastqc.mk

clean:
	rm -rf bwa picard gatk varscan snpeff filtered-variants tmp
	
#-----------	
# ALIGNMENT, SORTING, MARK DUPLICATES, INDEXING
#-----------

.PHONY: bwa
bwa: $(foreach S, $(SAMPLES), bwa/$S.bwa.sorted.dupmarked.bam.bai)
	 
bwa/%.bwa.bam: $(REFGENOME) $(FASTQ)/%_R1.txt.gz $(FASTQ)/%_R2.txt.gz
	mkdir -p bwa
	$(BWA) mem \
		-t 50 \
		-R '@RG\tID:RG1\tPL:Illumina\tSM:$*' \
		$(word 1,$^) $(word 2,$^) $(word 3,$^) | \
			$(SAMTOOLS) view -bhS - \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

bwa/%.bwa.sorted.bam: bwa/%.bwa.bam
	$(SAMTOOLS) sort -@ 5 -o $< $* 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	rm $<

bwa/%.bwa.sorted.dupmarked.bam: bwa/%.bwa.sorted.bam
	mkdir -p picard
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /data_synology/software/picard-tools-1.114/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	rm $<

bwa/%.bwa.sorted.dupmarked.bam.bai: bwa/%.bwa.sorted.dupmarked.bam
	rm -f $@
	$(SAMTOOLS) index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------
# GATK INDEL REALIGNMENT
#-----------

.PHONY: gatk
gatk: $(foreach S, $(SAMPLES), gatk/$S.bwa.sorted.dupmarked.realigned.bam)

gatk/%.intervalList.intervals: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai
	mkdir -p gatk
	$(GATK) \
		-T RealignerTargetCreator \
		-R $(REFGENOME) \
		-I $< \
		-o $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
   
gatk/%.bwa.sorted.dupmarked.realigned.bam: bwa/%.bwa.sorted.dupmarked.bam gatk/%.intervalList.intervals
	$(GATK) \
		-T IndelRealigner \
		-R $(REFGENOME) \
		-I $(word 1,$^) \
		-targetIntervals $(word 2,$^) \
		-o $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	rm -f gatk/%.bwa.sorted.dupmarked.realigned.bam.part.bai

gatk/all_normals.bwa.sorted.dupmarked.realigned.bam: $(foreach S, $(SAMPLES_NORMAL), gatk/$S.bwa.sorted.dupmarked.realigned.bam)
	$(SAMTOOLS) merge $@.part $^
	mv $@.part $@
	
gatk/%.bwa.sorted.dupmarked.realigned.bam.bai: gatk/%.bwa.sorted.dupmarked.realigned.bam
	rm -f $@
	$(SAMTOOLS) index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------	
# PICARD METRICS
#-----------	

.PHONY: picard
picard: picard/combined.picard.insert_size_metrics.tsv \
	    picard/combined.picard.alignment_summary_metrics.tsv \
	    picard/combined.picard.insert_size_histogram.pdf \
	    picard/combined.picard.base_distribution_by_cycle.pdf \
	    picard/combined.picard.hs_metrics.tsv

picard/%.multiplemetrics: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai
	mkdir -p picard && rm -f $@
	$(PICARD) CollectMultipleMetrics \
		INPUT=$< \
		OUTPUT=picard/$* \
		VALIDATION_STRINGENCY=LENIENT \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=MeanQualityByCycle \
		2>&1 | $(LOG)
	touch $@

picard/%.hs_metrics: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai $(COVERBED)
	mkdir -p picard && rm -f $@
	$(SAMTOOLS) view -H $< 2>&1 1> picard/$*.bed | $(LOG)
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,"+",$$4 }' $(word 3,$^) >> picard/$*.bed
	$(PICARD) CalculateHsMetrics \
		BAIT_INTERVALS=picard/$*.bed \
		TARGET_INTERVALS=picard/$*.bed \
		INPUT=$< \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=$(REFGENOME) \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.bed
	mv $@.part $@

picard/combined.picard.insert_size_metrics.tsv: $(foreach S, $(SAMPLES), picard/$S.multiplemetrics)
	awk 'NR==7 { print "SAMPLE\t" $$_}' $(subst .multiplemetrics,.insert_size_metrics,$<) > $@.part
	for FILE in $(subst .multiplemetrics,.insert_size_metrics,$^) ; do awk 'NR==8 { print FILENAME "\t" $$_}' "$$FILE" >> $@.part ; done
	mv $@.part $@

picard/combined.picard.alignment_summary_metrics.tsv: $(foreach S, $(SAMPLES), picard/$S.multiplemetrics)
	awk 'NR==7 { print "SAMPLE\t" $$_}' $(subst .multiplemetrics,.alignment_summary_metrics,$<) > $@.part
	for FILE in $(subst .multiplemetrics,.alignment_summary_metrics,$^) ; do awk 'NR==8||NR==9 {print FILENAME "\t" $$_}' "$$FILE" >> $@.part ; done
	mv $@.part $@

picard/combined.picard.insert_size_histogram.pdf: $(foreach S, $(SAMPLES), picard/$S.multiplemetrics)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $(subst .multiplemetrics,.insert_size_histogram.pdf,$^)
	mv $@.part $@

picard/combined.picard.base_distribution_by_cycle.pdf: $(foreach S, $(SAMPLES), picard/$S.multiplemetrics)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $(subst .multiplemetrics,.base_distribution_by_cycle.pdf,$^)
	mv $@.part $@

picard/combined.picard.hs_metrics.tsv: $(foreach S, $(SAMPLES), picard/$S.hs_metrics)
	awk 'NR==7' $< > $@.part
	for FILE in $^ ; do echo $$FILE && awk 'NR==8' "$$FILE" >> $@.part ; done
	mv $@.part $@
	
#-----------	
# VARSCAN
#-----------	

.PHONY: varscan
varscan: $(foreach P, $(PAIRS), varscan/$P.varscan.vcf) \
		 varscan/all_normals.varscan.vcf.gz

# tumor/normal comparison
varscan/%.varscan.vcf: gatk/%T.bwa.sorted.dupmarked.realigned.bam \
                       gatk/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.bwa.sorted.dupmarked.realigned.bam \
                       gatk/%T.bwa.sorted.dupmarked.realigned.bam.bai \
                       gatk/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.bwa.sorted.dupmarked.realigned.bam.bai \
                       $(VCFCONCAT) $(VCFSORT)
	mkdir -p varscan
	$(VARSCAN) somatic \
		<($(SAMTOOLS) mpileup -f $(REFGENOME) -q 10 -B --rf 2 --ff 3852  $(word 2,$^) $(word 1,$^)) \
		--mpileup 1 \
		--min-coverage 2 \
		--min-strands2 2 \
		--min-var-freq 0.1 \
		--normal-purity 1 \
		--tumor-purity 0.95 \
		--p-value 1 \
		--somatic-p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		--output-snp $@.snp \
		--output-indel $@.indel \
		2>&1 | $(LOG)
	$(VCFCONCAT) $@.snp $@.indel | $(VCFSORT) -c > $@.part
	rm $@.snp $@.indel
	mv $@.part $@

# all variants in combined normals for filtering of somatic mutations
varscan/all_normals.varscan.vcf.gz: gatk/all_normals.bwa.sorted.dupmarked.realigned.bam gatk/all_normals.bwa.sorted.dupmarked.realigned.bam.bai
	mkdir -p varscan
	$(VARSCAN) mpileup2cns \
		<($(SAMTOOLS) mpileup -f $(REFGENOME) -q 10 -B --rf 2 --ff 3852 $<) \
		--min-coverage 3 \
		--min-reads2 3 \
		--min-avg-qual 20 \
		--min-var-freq 0 \
		--p-value 1 \
		--strand-filter 0 \
		--output-vcf 1 \
		--variants 1 \
		2>&1 1>$@.part | $(LOG)
	bgzip -c $@.part > $@
	tabix -p vcf $@ 

#-----------	
# SNP calling
#-----------

snps/allsamples.vcf.gz:  $(foreach S, $(SAMPLES), gatk/$S.bwa.sorted.dupmarked.realigned.bam)
	mkdir -p snps
	$(SAMTOOLS131) mpileup -t DP,AD,ADF,ADR -ugf $(REFGENOME) $^ | \
		$(BCFTOOLS) call -vmO u | \
		$(BCFTOOLS) filter -i'%QUAL>10' -O z -o $@.part -
	mv $@.part $@
	tabix -p vcf $@
	
snps/allsamples.noMissing.vcf.gz: snps/allsamples.vcf.gz
	$(VCFTOOLS) --gzvcf $< --remove-indv testT --remove-indv testK --max-missing 1 --recode --stdout | bgzip -c > $@.part
	mv $@.part $@
	tabix -p vcf $@
  
#-----------	
# SNPEFF
#-----------

.PHONY: snpeff
snpeff: $(foreach P, $(PAIRS), snpeff/$P.varscan.dbsnp.snpeff.vcf.bgz.tbi)

snpeff/%.varscan.dbsnp.vcf: varscan/%.varscan.vcf $(DBSNP)
	mkdir -p snpeff
	$(SNPSIFT) annotate -v -a $(DBSNP) $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

snpeff/%.varscan.dbsnp.snpeff.vcf: snpeff/%.varscan.dbsnp.vcf
	mkdir -p snpeff
	$(SNPEFF) -stats snpeff/$*.snpeff.summary.html $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

snpeff/%.vcf.bgz.tbi: snpeff/%.vcf
	bgzip -c $< > snpeff/$*.vcf.bgz.part
	mv snpeff/$*.vcf.bgz.part snpeff/$*.vcf.bgz
	/data_synology/software/tabix-0.2.6/tabix -p vcf snpeff/$*.vcf.bgz

#-----------	
# FILTERED VARIANTS LIST
#-----------	

.PHONY: filtered-variants
filtered-variants: $(foreach P, $(PAIRS), filtered-variants/$P.tsv)

filtered-variants/%.tsv: snpeff/%.varscan.dbsnp.snpeff.vcf /mnt/projects/zohre/scripts/filter-variants.pl varscan/all_normals.varscan.vcf.gz
	mkdir -p filtered-variants
	perl /mnt/projects/zohre/scripts/filter-variants.pl \
		$< \
		--sample $* \
		--header \
		--rmsk-file /data_synology/max/zohre/ucsc/mm10.rmsk.nochr.sorted.txt.gz \
		--simpleRepeat-file /data_synology/max/zohre/ucsc/mm10.simpleRepeat.nochr.sorted.txt.gz \
		--segdup-file /data_synology/max/zohre/ucsc/mm10.genomicSuperDups.nochr.sorted.txt.gz \
		--normal-variants varscan/all_normals.varscan.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

#-----------	
# CNA
#-----------	
	
.PHONY: cna
cna: cna/allpatients.genome-coverage.pdf

cna/allpatients.genome-coverage.pdf: $(foreach P, $(filter-out test, $(PAIRS)), cna/$P.genome-coverage.pdf)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@
	rm $^

cna/%.snp-profile.pdf: cna/%T.genome-coverage.tsv \
               		   cna/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.genome-coverage.tsv \
               		   snpeff/%T.varscan.dbsnp.vcf \
               		   snpeff/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.varscan.dbsnp.vcf \
               		   /mnt/projects/generic/scripts/snp-profile.R
	mkdir -p snp-profile
	Rscript /mnt/projects/generic/scripts/snp-profile.R \
		--sample-id $* \
		--tumor-coverage $(word 1, $^) \
		--normal-coverage $(word 2, $^) \
		--tumor-vcf $(word 3, $^) \
		--normal-vcf $(word 4, $^) \
		--plot-output-file $@.part
	mv $@.part $@

cna/%.genome-coverage.pdf: cna/%T.genome-coverage.tsv \
               		       cna/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.genome-coverage.tsv \
						   cna/gcPerBin.tsv \
						   $(PLOTSEGCOV)
	mkdir -p cna
	Rscript $(PLOTSEGCOV) \
		--patient $* \
		--tumor $(word 1,$^) \
		--normal $(word 2,$^) \
		--gccontent $(word 3,$^) \
		--output $@.part
	mv $@.part $@

cna/%.genome-coverage.tsv: gatk/%.bwa.sorted.dupmarked.realigned.bam $(SEGMENTCOV) $(CHRSIZES)
	mkdir -p cna
	$(SAMTOOLS) depth -Q 10 $< | \
		perl $(SEGMENTCOV) \
			--sample $* \
			--bin-size 250000 \
			--chr-sizes $(CHRSIZES) \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	

cna/gcPerBin.tsv: $(REFGENOME) $(CHRSIZES)
	$(BEDTOOLS) nuc -fi $(REFGENOME) -bed <($(BEDTOOLS) makewindows -g $(CHRSIZES) -w 250000) | cut -f 1-2,5 | grep -v "^#" > $@.part
	mv $@.part $@

#cna/%.coverage.bedtools.txt: bwa/%.bwa.sorted.dupmarked.bam $(COVERBED)
#	mkdir -p cna
#	$(SAMTOOLS) view -bq 1 -F 3852 $< | \
#		$(BEDTOOLS) coverage -hist -abam - -b $(word 2, $^) | grep ^all > $@.part
#	mv $@.part $@
	
#----------
# sex
#----------

sex/allsamples.readcount.sexchr.txt: $(foreach S, $(SAMPLES), sex/$S.readcount.sexchr.txt)
	echo sample readsY readsTotal > $@.part
	cat $^ >> $@.part
	mv $@.part $@
	
sex/%.readcount.sexchr.txt: bwa/%.bwa.sorted.dupmarked.bam
	mkdir -p sex
	echo $* $$($(SAMTOOLS) view -c -f 2 -F 3840 -q 40 $< Y) $$($(SAMTOOLS) view -c -f 2 -F 3840 -q 40 $<) > $@.part
	mv $@.part $@