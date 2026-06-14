# MultiGen-PedSim — Simulation Framework for Stratified Multi-Generation Pedigrees
 
This simulation framework models multi-generational demographic histories and complex genetic architectures. Although developed to study maternal and offspring traits, the underlying code is fully generalizable and can be adapted to any pair of correlated traits or alternative populations.
 
---
 
## 📂 Repository Structure
 
```text
├── simulation.R                      # Main simulation script
├── blupf90_preprocessing.R           # Preprocessing for BLUPF90 analysis
├── snakefile                         # Snakemake workflow manager
├── src/
│   ├── founderHaplotypesGenerator.R  # Coalescent simulation via AlphaSimR/runMacs2
│   ├── pop_generation.R              # Multi-generation population generator
│   ├── subsetting.R                  # Greedy unrelated sampling and pedigree pruning
│   ├── helpers.R                     # Mathematical helper functions
│   └── my_writePlink.R               # Modified PLINK exporter for large genotype matrices
├── test/
│   └── simulation_test.R             # Small-scale local simulation for testing
└── README.md
```
 
---
 
## Prerequisites & Setup
 
### R Environment
 
Ensure you have R installed along with the following required packages:
 
```r
install.packages(c("AlphaSimR", "MASS", "Matrix", "data.table", "pedigree", "MCMCglmm"))
```
 
### Workflow Dependencies
 
- **Snakemake**: For workflow orchestration.
- **BLUPF90 Suite**: `renumf90`, `preGSf90`, and `gibbsf90+` must be accessible in your system `PATH` if running variance component estimation steps.
---
 
## How to Use
 
### 1. Running via Snakemake
 
The Snakefile can be adapted to work in any HPC or bash environment. The only modification required is the addition of the conda environments. Each rule can be modified based on the available resources and the size of the simulation.
 
The genetic component parameters governing the phenotype simulation can be changed directly through the Snakefile. Modifying the `params` structure enables parallelization of different genetic component scenarios and multiple BLUPF90 analysis runs.
 
One run of the workflow consists of: population generation, phenotype simulation, and BLUPF90 analysis across four model configurations — pedigree-only, univariate ssGBLUP for trait 1, univariate ssGBLUP for trait 2, and a bivariate ssGBLUP model.
 
### 2. Running the R Simulation Directly
 
`simulation.R` accepts 12 positional command-line arguments, making it easy to scale simulations on HPC clusters:
 
```bash
Rscript simulation.R [run_id] [prevCMC] [prevPregCMC] [prevASC] [male_ratio] [betaPreg] [var_d] [var_m] [var_c] [rho_CMC_d] [rho_CMC_m] [rho_d_m]
```
 
| Parameter | Description | Default |
| --- | --- | --- |
| `run_id` | Unique identifier, also used as the random seed | `1` |
| `prevCMC_mother` | Lifetime prevalence of the maternal condition | `0.10` |
| `prevPregCMC` | Prevalence of pregnancy-related episodes per mother | `0.16` |
| `prevASC_child` | Population prevalence of the child trait (both sexes combined) | `0.02` |
| `ASC_male_ratio` | Male-to-female ratio for diagnostic asymmetry | `5.0` (5:1) |
| `betaPreg` | Direct effect of pregnancy exposure on child liability | `0.35` |
| `var_d` | Proportion of child variance from direct additive genetic effects | `0.20` |
| `var_m` | Proportion of variance from maternal additive effects (genetic nurture) | `0.10` |
| `var_c` | Proportion of variance from shared maternal environment | `0.10` |
| `rho_CMC_d` | Correlation between maternal trait and child direct SNP effects | `0.30` |
| `rho_CMC_m` | Correlation between maternal trait and maternal nurture SNP effects | `0.20` |
| `rho_d_m` | Correlation between child direct and maternal indirect genetic components | `0.10` |
 
Running `simulation.R` directly is useful when parallelization across multiple scenarios is not needed, or when modifying the demographic characteristics of the simulated population. The core idea of the workflow is to generate one population and test BLUPF90 heritability estimates under different genetic component configurations. The script is commented to guide any modifications to the population parameters.
 
The `test/` folder contains a small-scale version of the simulation that can be run locally without Snakemake or BLUPF90, useful for getting familiar with the simulation features before running the full pipeline.
 
---
 
## Extending the Code for Your Own Population
 
This framework is designed to be a customizable template for your own research:
 
### Changing Demographics
 
Open `simulation.R` and look for the vectors `mothers_fraction`, `lambdaKids`, `p_new_partner`, and `mean_unions_dad`. These control the demographic behaviour of each successive generation. By default, they contain five entries (5 generations) with parameters estimated from Statistics Denmark register data. These arrays can be expanded or contracted to model alternative populations.
 
### Modifying the Traits
 
- Change the variable parameters in the Snakefile or via command-line arguments to model any maternal trait — such as maternal autoimmune status, educational attainment, or socioeconomic context — and any offspring outcome, such as ADHD, birthweight, or cardiovascular risk.
- Adjust the genetic parameters (`var_d`, `var_m`, and the correlation values) to shift from a purely direct genetic model to one dominated by environmental effects or genetic nurture.
---
 
## Citation
 
If you use this framework in your research, please cite the associated thesis and any relevant dependencies (AlphaSimR, BLUPF90, Snakemake).
 
---
 
*This project is open for use and adaptation for research purposes.*