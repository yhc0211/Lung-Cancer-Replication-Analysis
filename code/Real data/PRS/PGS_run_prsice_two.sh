#!/bin/bash
# Run PRSice per chromosome for each of the three SNP sets:
#   gwas (ILCCO GWAS-significant)
#   meta (ILCCO + MVP meta-analysis)
#   rep  (Replication / csmGmm)
#
# - Continues past per-chromosome failures (no `set -e`)
# - Auto-retries with --extract <out>.valid when PRSice flags duplicate SNP IDs
# - Skips (tag, chr) combos whose .bed doesn't exist

set -uo pipefail   # NOTE: no -e, so one bad chr doesn't kill the whole loop

# Load R if needed (uncomment and adjust the version that `module avail R` shows)
# module load R/4.3.1

DATA=/home/ychang11/csmGmm_reproduce/Lung/PGS/data
PRSICE_R=/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung/PGS/PRSice/PRSice.R
PRSICE_BIN=/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung/PGS/PRSice/bin/PRSice

PHENO=${DATA}/pheno.txt
COV=${DATA}/cov.txt
COV_COLS=age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

run_one () {
  local tag=$1 chr=$2
  local target=${DATA}/${tag}_ukb_chr${chr}_filtered
  local out=${DATA}/${tag}_prs_chr${chr}_r2

  if [[ ! -f ${target}.bed ]]; then
    echo "skip ${tag} chr${chr} (no .bed)"
    return
  fi

  local base_args=(
    --prsice ${PRSICE_BIN}
    --base   ${DATA}/base_${tag}.txt
    --target ${target}
    --out    ${out}
    --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2
    --stat BETA --beta --pvalue P
    --pheno  ${PHENO} --pheno-col Phenotype
    --cov    ${COV}   --cov-col ${COV_COLS}
    --binary-target T
    --clump-kb 250 --clump-r2 0.5
    --bar-levels 0.05,0.1,1
    --all-score --fastscore --thread 8 --ignore-fid
  )

  echo "==> ${tag} chr${chr}"
  Rscript ${PRSICE_R} "${base_args[@]}" || true

  # PRSice writes <out>.valid when it detects duplicate SNP IDs, then halts.
  # If .all_score wasn't produced but .valid was, retry with --extract.
  if [[ -f ${out}.valid && ! -f ${out}.all_score ]]; then
    echo ">> ${tag} chr${chr}: duplicates found, retrying with --extract ${out}.valid"
    Rscript ${PRSICE_R} "${base_args[@]}" --extract ${out}.valid || true
  fi

  if [[ -f ${out}.all_score ]]; then
    echo "ok ${tag} chr${chr}"
  else
    echo "!! ${tag} chr${chr} FAILED — see ${out}.log"
  fi
}

for tag in gwas meta rep; do
  for chr in $(seq 1 22); do
    run_one ${tag} ${chr}
  done
done

echo "All done."