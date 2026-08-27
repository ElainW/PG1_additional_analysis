# Regenerate every analysis output from scratch, in dependency order.
#   make all     run scripts 00..16 (each reads the previous ones' results/)
#   make clean   remove generated results/ and figures/ (keeps data/ cache)
# Needs: python3 with pandas numpy scipy matplotlib seaborn (biopython for 17),
# and the tracked data/ cache (UniProt for 08/15; 1V4S.pdb for 17). No network.

PY = python3
SD = scripts

SCRIPTS = \
	00_check_directionality \
	01_compute_spearman \
	02_spearman_by_taxa \
	03_persite_ranking \
	04_subset_figs \
	05_fig6_familybest \
	06_hxk4_site_disagreement \
	07_w99_deepdive \
	08_gck_structure_map \
	09_w99_subset_figs \
	10_assay_overlap \
	11_hxk4_lowagreement_batch \
	12_abundance_physchem \
	13_agreement_vs_accuracy \
	14_quartile_by_agreement \
	15_functional_disagreement_generalization \
	16_collapse_modeltype_withinsite \
	17_gck_3d_distances \
	21_generalization_candidates \
	22_disease_gof_lof

.PHONY: all clean

all:
	@echo "=== run started: $$(date) ===" | tee RUN.log
	@for s in $(SCRIPTS); do \
		echo ">>> $$s.py  ($$(date +%H:%M:%S))" | tee -a RUN.log; \
		( cd $(SD) && $(PY) $$s.py ) >> ../RUN.log 2>&1 \
			|| { echo "!!! FAILED: $$s.py (see RUN.log)"; exit 1; }; \
	done
	@echo "=== run finished: $$(date) ===" | tee -a RUN.log

clean:
	rm -f results/*.csv
	rm -rf figures
