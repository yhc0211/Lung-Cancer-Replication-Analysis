# Introduction

Code to produce a model-based replication analysis for lung cancer GWAS summary statistics.


# Code Documentation

This folder contains the analysis code for the lung cancer replication analysis project.

The code is organized into two main components:

1. `Simulation/`: simulation studies comparing the proposed replication method with competing approaches.
2. `Real data/`: real lung cancer GWAS analyses, including two-cohort and three-cohort replication analyses, `repfdr` comparisons, incongruity checks, PRS analyses, and FAVOR annotation analyses.

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
