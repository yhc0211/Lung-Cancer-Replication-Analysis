# Introduction

Code to produce a model-based replication analysis for lung cancer GWAS summary statistics.


# Code Documentation

This folder contains the analysis code for the lung cancer replication analysis project.

The code is organized into two main components:

1. `Simulation/`: simulation studies comparing the proposed replication method with competing approaches.
2. `Real data/`: real lung cancer GWAS analyses, including two-cohort and three-cohort replication analyses, `repfdr` and kernel method comparisons, incongruity checks, PRS analyses, and FAVOR annotation analyses.

---

## Folder Structure

```text
code/
├── Simulation/
│   ├── plot_simulation.R
│   ├── sim_three.R
│   ├── sim_three.lsf
│   ├── sim_two.R
│   └── sim_two.lsf
│
└── Real data/
    ├── FAVOR/
    │   ├── Plot_FAVOR_LUSC_revised.R
    │   └── Plot_FAVOR_Overall_revised.R
    │
    ├── PRS/
    │   ├── PGS_run_prsice_three.sh
    │   ├── PGS_run_prsice_two.sh
    │   ├── PRS_revised_three.R
    │   ├── PRS_revised_three.lsf
    │   ├── PRS_revised_two.R
    │   └── PRS_revised_two.lsf
    │
    ├── three cohort/
    │   ├── plot_three_rep.R
    │   ├── three_incon_adeno_scc.R
    │   ├── three_incon_adeno_scc.lsf
    │   ├── three_incon_overall.R
    │   ├── three_incon_overall.lsf
    │   ├── three_rep_adeno.R
    │   ├── three_rep_adeno.lsf
    │   ├── three_rep_overall.R
    │   ├── three_rep_overall.lsf
    │   ├── three_rep_scc.R
    │   ├── three_rep_scc.lsf
    │   ├── three_repfdr_adeno.R
    │   ├── three_repfdr_adeno.lsf
    │   ├── three_repfdr_overall.R
    │   ├── three_repfdr_overall.lsf
    │   ├── three_repfdr_scc.R
    │   └── three_repfdr_scc.lsf
    │
    └── two cohort/
        ├── plot_two_rep.R
        ├── plot_two_rep.lsf
        ├── two_incon_overall.R
        ├── two_incon_overall.lsf
        ├── two_rep_adeno.R
        ├── two_rep_adeno.lsf
        ├── two_rep_overall.R
        ├── two_rep_overall.lsf
        ├── two_rep_scc.R
        ├── two_rep_scc.lsf
        ├── two_repfdr_adeno.R
        ├── two_repfdr_adeno.lsf
        ├── two_repfdr_overall.R
        ├── two_repfdr_overall.lsf
        ├── two_repfdr_scc.R
        └── two_repfdr_scc.lsf
```
## Analysis Workflow

- **Run simulation studies**
  - `Simulation/sim_two.R`: runs the two-cohort simulation study.
  - `Simulation/sim_three.R`: runs the three-cohort simulation study.

- **Summarize simulation results**
  - `Simulation/plot_simulation.R`: aggregates simulation outputs and generates simulation summary figures/tables.

- **Run two-cohort real-data replication analyses**
  - `Real data/two cohort/two_rep_overall.R`: overall lung cancer analysis.
  - `Real data/two cohort/two_rep_adeno.R`: lung adenocarcinoma analysis.
  - `Real data/two cohort/two_rep_scc.R`: lung squamous cell carcinoma analysis.

- **Run three-cohort real-data replication analyses**
  - `Real data/three cohort/three_rep_overall.R`: overall lung cancer analysis.
  - `Real data/three cohort/three_rep_adeno.R`: lung adenocarcinoma analysis.
  - `Real data/three cohort/three_rep_scc.R`: lung squamous cell carcinoma analysis.

- **Run comparison analyses using `repfdr`**
  - `two_repfdr_*.R`: applies `repfdr` & kernel method in the two-cohort setting.
  - `three_repfdr_*.R`: applies `repfdr` & kernel method in the three-cohort setting.

- **Check incongruous SNP ranking patterns**
  - `two_incon_overall.R`: checks incongruity in the two-cohort overall lung cancer analysis.
  - `three_incon_overall.R`: checks incongruity in the three-cohort overall lung cancer analysis.
  - `three_incon_adeno_scc.R`: checks incongruity in subtype-specific three-cohort analyses.

- **Generate final real-data figures and tables**
  - `plot_two_rep.R`: generates two-cohort real-data result summaries.
  - `plot_three_rep.R`: generates three-cohort real-data result summaries.

- **Run downstream PRS analysis**
  - Scripts in `Real data/PRS/` run PRSice and summarize polygenic risk score results.

- **Run FAVOR annotation visualization**
  - Scripts in `Real data/FAVOR/` generate functional annotation plots for selected SNPs.
