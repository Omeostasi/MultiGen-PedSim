# Project Overview and Quick Start Guide
This simulation framework models multi-generational demographic histories and complex genetic architectures. Although developed to study maternal and offspring traits, the underlying code is fully generalized and can be adapted to any pair of correlated traits or alternative populations.

## 📂 Repository Structure

```text
├── simulation.R                  # Main script
├── blupf90_preprocessing.R       # Prepropreccis for BLUPF90 analysis
├── snakefile                     # Snakemake workflow manager
├── src/
│   ├── founderHaplotypesGenerator.R  # Coalescence simulation via AlphaSimR/runMacs2
│   ├── pop_generation.R              # Generate different generations
│   ├── subsetting.R                  # Greedy unrelated sampling & pedigree pruning algorithms
│   ├── helpers.R                     # Mathematical helper functions
│   └── my_writePlink.R               # Modified PLINK exporter
|-- test/
|   |-- simulation_test.R             # quickSimulation
└── README.md                     

```

---

##  Prerequisites & Setup

### R Environment

Ensure you have R installed along with the following required packages:

```R
install.packages(c("AlphaSimR", "MASS", "Matrix", "data.table", "pedigree", "MCMCglmm"))

```

### Workflow Dependencies

* **Snakemake**: For workflow orchestration.
* **BLUPF90 Suite**: (`renumf90`, `preGSf90`, `gibbsf90+`) must be accessible in your system `PATH` if executing quantitative genetic variance estimation steps.

---

##  How to Use

### 1. Running via Snakemake 

The snakefile can be arranged to work in another HPC or bash envoirnment.
The only modification required is the addition of the conda envoirnments. 
Each rule can be modified based on the resources avaialble and the size of the simulation. 
Through the snakemake file is possible to change the parameters governing the phenotype simulation and the genetic components of the population. 
CHanging the params structure inside the snakefile allow parallelizaiton of different genetic component simulations and different runs of the BLUPF90 analysis. 
One run of the workflow consinst of the generaiton of the popultion, the phenotype simulation, the blupf90anlysis for 4 different models: pedigree only, univariate ssgblup for trait 1, univariate ssgblup for trait 2 and bivariate model ssgblup. 

### 2. Running the R Simulation Directly

`simulation.R` accepts 12 positional command-line arguments to easily scale simulations on high-performance computing (HPC) clusters:

```bash
Rscript simulation.R [run_id] [prevCMC] [prevPregCMC] [prevASC] [male_ratio] [betaPreg] [var_d] [var_m] [var_c] [rho_CMC_d] [rho_CMC_m] [rho_d_m]

```

| Parameter | Description | Default Value |
| --- | --- | --- |
| `run_id` | Unique identifier used also as the random seed | `1` |
| `prevCMC_mother` | Lifetime prevalence of the maternal condition | `0.10` |
| `prevPregCMC` | Prevalence of pregnancy-related episodes per mother | `0.16` |
| `prevASC_child` | Population prevalence of child trait (both sexes combined) | `0.02` |
| `ASC_male_ratio` | Male-to-female ratio for diagnostic asymmetry | `5.0` (5:1) |
| `betaPreg` | Direct impact of pregnancy exposure on child liability | `0.35` |
| `var_d` | Proportion of child variance from direct additive genetic effects | `0.20` |
| `var_m` | Proportion of variance from maternal additive effects (nurture) | `0.10` |
| `var_c` | Proportion of variance from shared maternal environment | `0.10` |
| `rho_CMC_d` | Correlation between maternal trait and child direct SNPs | `0.30` |
| `rho_CMC_m` | Correlation between maternal trait and maternal nurture SNPs | `0.20` |
| `rho_d_m` | Correlation between child direct and maternal indirect genetic components | `0.10` |


The R script can be run directly if there is no need for parallelization of multiple different population/phenotype simulations.
Another purpuose of using the R script directly is to change the demographic charateristics of the simulated population.
The idea behind the workflow is to generate 1 population under different genetic components to test the heritbility estiamtes of  BLUPF90,
but if there is the possibility of changing important characteristic of the population and create different stratified generations through direct modifications of the parametrs.
The simulation.R script explains what to modify.

Inside the the test folder is available a script which is a small scale simulation that can be run locally without snakemake and blupf90 to test and get used to the simulation features. 
It's just a smaller scale of the simulation with smaller values. 
---

## Extending the Code for Your Own Population

This software is built explicitly to be a customizable template for your own research projects:

### Changing Demographics

Open `simulation.R` and look for the vectors `mothers_fraction`, `lambdaKids`, `p_new_partner`, and `mean_unions_dad`. These control the behavior of each successive generation. By default, they contain 5 entries corresponding to the **Greatest, Silent, Boomers, GenX, and Millennial** cohorts estimated from Statistics Denmark. You can expand or contract these arrays to model alternative global populations.

### Modifying the Traits

* Change the variable parameters in the `snakefile` or command line script to model any maternal trait (e.g., maternal autoimmune status, educational attainment, or socioeconomic context) impacting an offspring outcome (e.g., ADHD, birthweight, or cardiovascular risk).
* Adjust the genetic parameters ($\text{var}_d$, $\text{var}_m$, correlations) to shift from a strictly direct genetic model to one dominated by environment or genetic nurture.

---

## Feel fre to use it