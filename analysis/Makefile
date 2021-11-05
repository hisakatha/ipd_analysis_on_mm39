SHELL := bash
.DELETE_ON_ERROR:
.SECONDEXPANSION:
all: $$(targets)
.PHONY: all

DSLOP := 20
FA := /glusterfs/hisakatha/ensembl_mm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
FAI := $(FA:=.fai)

fst_wt := load_and_save.wu2016_wt.R.done
fst_ko := load_and_save.wu2016_ko.R.done

$(fst_wt):
	Rscript load_and_save.wu2016_wt.R && touch $@

$(fst_ko):
	Rscript load_and_save.wu2016_ko.R && touch $@
targets += $(fst_wt) $(fst_ko)

# Avoid parallel execution of heavy jobs
$(fst_wt): | $(fst_ko)

deep_kinetics_wt := deep_kinetics_region.wu2016_wt.bed
deep_kinetics_ko := deep_kinetics_region.wu2016_ko.bed
DEEP_KINETICS := $(deep_kinetics_wt) $(deep_kinetics_ko)

$(deep_kinetics_wt): $(fst_wt)
	Rscript convert_kinetics_deep_region_into_bed.wu2016_wt.R
$(deep_kinetics_ko): $(fst_ko)
	Rscript convert_kinetics_deep_region_into_bed.wu2016_ko.R

DEEP_KINETICS_MERGED_SORTED := $(DEEP_KINETICS:=.merged.sorted)
$(DEEP_KINETICS_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_SLOP := $(DEEP_KINETICS_MERGED_SORTED:=.slop$(DSLOP).bed)
$(DEEP_KINETICS_SLOP): %.slop$(DSLOP).bed: %
	bedtools slop -i $< -g $(FAI) -b $(DSLOP) > $@ && touch $@
DEEP_KINETICS_SLOP_FA := $(DEEP_KINETICS_SLOP:=.fa)
$(DEEP_KINETICS_SLOP_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) > $@ && touch $@
# For counting the size of the regions
DEEP_KINETICS_SLOP_MERGED := $(DEEP_KINETICS_SLOP:=.merged)
$(DEEP_KINETICS_SLOP_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $@ && touch $@
targets += $(DEEP_KINETICS_SLOP_FA) $(DEEP_KINETICS_SLOP_MERGED)

m6A_gff_wt := ipd_summary.m6A_cov25.wu2016_wt.gff
$(m6A_gff_wt): ../wu2016_wt_no_chunk/ipd_summary.m6A_cov25.gff
	ln -s $< $@
m6A_gff_ko := ipd_summary.m6A_cov25.wu2016_ko.gff
$(m6A_gff_ko): ../wu2016_ko_no_chunk/ipd_summary.m6A_cov25.gff
	ln -s $< $@
m6A_context_wt := ipd_summary.m6A_cov25.wu2016_wt.gff.fa
$(m6A_context_wt): ../wu2016_wt_no_chunk/ipd_summary.m6A_cov25.gff.fa
	ln -s $< $@
m6A_context_ko := ipd_summary.m6A_cov25.wu2016_ko.gff.fa
$(m6A_context_ko): ../wu2016_ko_no_chunk/ipd_summary.m6A_cov25.gff.fa
	ln -s $< $@
targets += $(m6A_gff_wt) $(m6A_gff_ko) $(m6A_context_wt) $(m6A_context_ko)

MOD := $(m6A_gff_wt) $(m6A_gff_ko)
MSLOP := 20

MOD_SLOP := $(MOD:.gff=.slop$(MSLOP).gff)
$(MOD_SLOP): %.slop$(MSLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(MSLOP) > $@
MOD_FULL := $(MOD_SLOP:.gff=.fullLength.gff)
$(MOD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(MSLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@
MOD_FULL_FA := $(MOD_FULL:=.fa)
$(MOD_FULL_FA): %.fullLength.gff.fa: %.fullLength.gff
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
MOD_FULL_MERGED := $(MOD_FULL:=.merged)
$(MOD_FULL_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 3,6,7 -o distinct,mean,distinct > $@

targets += $(MOD_FULL_FA) $(MOD_FULL_MERGED)
