export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
.SECONDEXPANSION: # allow functions in dependency list
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=/mnt/projects/zohre
FASTQ=/data_synology_nfs/max/zohre/exome_maus/fastq
REFGENOME=/data_synology/max/zohre/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
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
PICARD=$(DOCKER) biowaste:5000/ccri/picard-2.2.2 java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /usr/picard/picard.jar
GATK=java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar /data_synology/software/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar
VARSCAN=java -jar /data_synology/software/varscan-2.3.6/VarScan.v2.3.6.jar
SNPEFF=java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar /data_synology/software/snpeff-4.2/snpEff.jar

#TUMOR_0=11291 11682 11689 11232 11221Up    
#TUMOR_1=11746 9193
#TUMOR_2=11746 9193
#SAMPLES=$(foreach S, $(TUMOR_0), $(S)T $(S)K) \
#        $(foreach S, $(TUMOR_1), $(S)_1T $(S)K) \
#        $(foreach S, $(TUMOR_2), $(S)_2T $(S)K)

SAMPLES=test 11291 11682 11689 11232 11221Up 11746_1 9193_1 11746_2 9193_2

all: bwa picard gatk varscan snpeff filtered-variants

#include /mnt/projects/generic/scripts/rna-seq/fastqc.mk

clean:
	rm -rf bwa picard gatk varscan snpeff filtered-variants tmp
	
#-----------	
# ALIGNMENT, SORTING, MARK DUPLICATES, INDEXING
#-----------

.PHONY: bwa
bwa: $(foreach S, $(SAMPLES), bwa/$(S)T.bwa.sorted.dupmarked.bam.bai bwa/$$(subst Up,,$$(subst _2,,$$(subst _1,,$S)))K.bwa.sorted.dupmarked.bam.bai)
	 
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
gatk: $(foreach S, $(SAMPLES), gatk/$(S)T.bwa.sorted.dupmarked.realigned.bam gatk/$$(subst Up,,$$(subst _2,,$$(subst _1,,$S)))K.bwa.sorted.dupmarked.realigned.bam)

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
	
gatk/%.bwa.sorted.dupmarked.realigned.bam.bai: gatk/%.bwa.sorted.dupmarked.realigned.bam
	rm -f $@
	$(SAMTOOLS) index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------	
# PICARD METRICS
#-----------	

.PHONY: picard
picard: $(foreach S, $(SAMPLES), picard/$(S)T.multiplemetrics   picard/$$(subst Up,,$$(subst _2,,$$(subst _1,,$S)))K.multiplemetrics \
                                 picard/$(S)T.hs_metrics        picard/$$(subst Up,,$$(subst _2,,$$(subst _1,,$S)))K.hs_metrics)

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

#-----------	
# VARSCAN
#-----------	

.PHONY: varscan
varscan: $(foreach S, $(SAMPLES), varscan/$S.varscan.vcf)

varscan/%.varscan.vcf: gatk/%T.bwa.sorted.dupmarked.realigned.bam \
                       gatk/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.bwa.sorted.dupmarked.realigned.bam \
                       gatk/%T.bwa.sorted.dupmarked.realigned.bam.bai \
                       gatk/$$(subst Up,,$$(subst _2,,$$(subst _1,,%)))K.bwa.sorted.dupmarked.realigned.bam.bai
	mkdir -p varscan
	$(VARSCAN) somatic \
		<($(SAMTOOLS) view -u -q 10 -F 1024 -f 2 $(word 1,$^) | $(SAMTOOLS) mpileup -f $(REFGENOME) -) \
		<($(SAMTOOLS) view -u -q 10 -F 1024 -f 2 $(word 2,$^) | $(SAMTOOLS) mpileup -f $(REFGENOME) -) \
		varscan/$*.varscan.part \
		--min-coverage 2 \
		--min-strands2 2 \
		--min-var-freq 0.2 \
		--normal-purity 1 \
		--tumor-purity 0.95 \
		--p-value 1 \
		--somatic-p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		2>&1 | $(LOG)
	mv varscan/$*.varscan.part.snp.vcf varscan/$*.varscan.snp.vcf
	mv varscan/$*.varscan.part.indel.vcf varscan/$*.varscan.indel.vcf

#-----------	
# SNPEFF
#-----------

.PHONY: snpeff
snpeff: $(foreach S, $(SAMPLES), snpeff/$S.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi)

snpeff/%.varscan.dbsnp.vcf: varscan/%.varscan.vcf /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	mkdir -p snpeff
	PWD=$(pwd)
	(cd /data_synology/software/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar annotate \
		-v /mnt/projects/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.varscan.dbsnp.snpeff.vcf: snpeff/%.varscan.dbsnp.vcf
	PWD=$(pwd)
	(cd /data_synology/software/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar snpEff.jar -v -stats $(PWD)/snpeff/$*.snpeff.summary.html GRCm38.82 $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.varscan.dbsnp.snpeff.dbNSFP.vcf: snpeff/%.varscan.dbsnp.snpeff.vcf /mnt/projects/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz
	PWD=$(pwd)
	(cd /data_synology/software/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar dbnsfp \
		-v /mnt/projects/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz \
		-collapse \
		-f SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,SiPhy_29way_logOdds,LRT_pred,LRT_score,MutationTaster_pred,MutationTaster_score,MutationAssessor_pred,MutationAssessor_score,FATHMM_pred,FATHMM_score,RadialSVM_pred,RadialSVM_score,GERP++_RS,1000Gp1_AF,1000Gp1_AFR_AF,1000Gp1_EUR_AF,1000Gp1_AMR_AF,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,Uniprot_acc,Interpro_domain, \
		$(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.vcf.bgz.tbi: snpeff/%.vcf
	bgzip -c $< > snpeff/$*.vcf.bgz.part
	mv snpeff/$*.vcf.bgz.part snpeff/$*.vcf.bgz
	/data_synology/software/tabix-0.2.6/tabix -p vcf snpeff/$*.vcf.bgz

#---------------
# COVERAGE PLOT
#---------------

coverage/coverage.png: $(foreach S, $(SAMPLES_SURESELECTXT), coverage/$S.coverage.bedtools.txt) /mnt/projects/oskar/scripts/coverage-plot.R
	Rscript /mnt/projects/oskar/scripts/coverage-plot.R
	mv $@.part $@
	
coverage/%.coverage.bedtools.txt: bwa/%.bwa.sorted.dupmarked.bam /mnt/projects/oskar/data/S06588914_Covered.nochr.bed
	mkdir -p coverage
	$(SAMTOOLS) view -bq 1 -F 3852 $< | /data_synology/software/bedtools-2.17.0/bin/bedtools coverage -hist -abam - -b $(word 2, $^) | grep ^all > $@.part
	mv $@.part $@


#-----------	
# FILTERED VARIANTS LIST
#-----------	

.PHONY: filtered-variants
filtered-variants: $(foreach S, $(SAMPLES), filtered-variants/$S.tsv)

filtered-variants/%.tsv: snpeff/%.varscan.dbsnp.snpeff.dbNSFP.vcf /mnt/projects/oskar/scripts/filter-variants.pl
	mkdir -p filtered-variants
	perl /mnt/projects/oskar/scripts/filter-variants.pl \
		$< \
		--patient $* \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file /mnt/projects/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--exac-file /mnt/projects/generic/data/ExAC/ExAC.r0.2.sites.vep.vcf.gz \
		--clinvar-file /mnt/projects/generic/data/clinvar/clinvar_20140929.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
	
filtered-variants/StA.tsv: snpeff/StA_BS.varscan.dbsnp.snpeff.dbNSFP.vcf snpeff/StA_ES.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi snpeff/StA_SS.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi /mnt/projects/oskar/scripts/filter-variants.pl
	mkdir -p filtered-variants
	perl /mnt/projects/oskar/scripts/filter-variants.pl \
		$< \
		--patient StA_BS \
		--mother snpeff/StA_ES.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz \
		--father snpeff/StA_SS.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz \
		--rmsk-file /mnt/projects/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file /mnt/projects/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file /mnt/projects/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file /mnt/projects/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible /mnt/projects/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro /mnt/projects/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file /mnt/projects/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file /mnt/projects/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file /mnt/projects/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--exac-file /mnt/projects/generic/data/ExAC/ExAC.r0.2.sites.vep.vcf.gz \
		--clinvar-file /mnt/projects/generic/data/clinvar/clinvar_20140929.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@

#------------
# PUBLISH ON WEB SERVER FOR MEDGEN TO DOWNLOAD
#------------

published/%: filtered-variants/%.tsv picard/%.agilent-sureselect.hs_metrics gatk/%.bwa.sorted.dupmarked.realigned.bam
	mkdir -p published
	scp $^ cf@biotrash:/var/www/html/medgen
	ssh cf@biotrash 'mv /var/www/html/medgen/$*.tsv /var/www/html/medgen/$*.xls'
	touch $@