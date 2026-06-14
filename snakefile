# PARAMETERS
PARAMS = {
    1: {
        "prevCMC_mother": 0.10,
        "prevPregCMC": 0.16,
        "prevASC_child": 0.02,
        "ASC_male_ratio": 5,
        "betaPreg": 0.35,
        "var_d": 0.0.20,
        "var_m": 0.10,
        "var_c": 0.10,
        "rho_CMC_d": 0.30,
        "rho_CMC_m": 0.20,
        "rho_d_m": 0.10,
    },
    # Add more runs here if needed, create another dictionary with the same keys and different values,
    # and add the run_id to the expand function in the all rule.
    2: {
        "prevCMC_mother": 0.10,
        "prevPregCMC": 0.16,
        "prevASC_child": 0.02,
        "ASC_male_ratio": 5,
        "betaPreg": 0.35,
        "var_d": 0.35,
        "var_m": 0.10,
        "var_c": 0.10,
        "rho_CMC_d": 0.30,
        "rho_CMC_m": 0.20,
        "rho_d_m": 0.10,
    }
}

run_id = list(PARAMS.keys()) # key must be integer for the run ids

# This checks that everything has been done and it gives the final output the workflow will produce. 
# Leads to the the rules that produce the final output. If the final output is missing, it will check which rule is missing and run it.
rule all:
    input:
        expand("results/job_done_run_{n}.txt", n=run_id)



# Rules called simulation_run_id are the ones that will produce the population and/or the phenotypes files.
# Order args to the R script are: run_id, prevCMC_mother, prevPregCMC, prevASC_child, ASC_male_ratio, betaPreg, var_d, var_m, var_c, rho_CMC_d, rho_CMC_m, rho_d_m.
# If not given they will go to default
rule simulation_run:
    output:
        "results/simulation_run_{n}_done.txt"
    resources:
        mem_mb=1200000,
        disk_mb=200000,
        runtime=7000,
        cpus_per_task=10,
    conda:
        # Conda envoirnment with AlphaSimR, Tidyverse, 
    log:
        "logs/sim_run_{n}.log"
    params:
        prevCMC_mother=lambda wc: PARAMS[int(wc.n)]["prevCMC_mother"],
        prevPregCMC=lambda wc: PARAMS[int(wc.n)]["prevPregCMC"],
        prevASC_child=lambda wc: PARAMS[int(wc.n)]["prevASC_child"],
        ASC_male_ratio=lambda wc: PARAMS[int(wc.n)]["ASC_male_ratio"],
        betaPreg=lambda wc: PARAMS[int(wc.n)]["betaPreg"],
        var_d=lambda wc: PARAMS[int(wc.n)]["var_d"],
        var_m=lambda wc: PARAMS[int(wc.n)]["var_m"],
        var_c=lambda wc: PARAMS[int(wc.n)]["var_c"],
        rho_CMC_d=lambda wc: PARAMS[int(wc.n)]["rho_CMC_d"],
        rho_CMC_m=lambda wc: PARAMS[int(wc.n)]["rho_CMC_m"],
        rho_d_m=lambda wc: PARAMS[int(wc.n)]["rho_d_m"]
    threads: 10
    shell:
        """
        date > {log} 2>&1

        Rscript simulation.R \
            {wildcards.n} \
            {params.prevCMC_mother} \
            {params.prevPregCMC} \
            {params.prevASC_child} \
            {params.ASC_male_ratio} \
            {params.betaPreg} \
            {params.var_d} \
            {params.var_m} \
            {params.var_c} \
            {params.rho_CMC_d} \
            {params.rho_CMC_m} \
            {params.rho_d_m} \
            >> {log} 2>&1 && touch {output}
        """

# This rule produces all input files for blupf90:
# pedigree, phenotype files, and all renumf90 par files for all BLUP folders.
rule blupf90_pre_run:
    output:
        "results/blupf90_pre_processing_run_{n}_done.txt"
    input:
        "results/simulation_run_{n}_done.txt"
    resources:
        mem_mb=20000,
        disk_mb=50000,
        runtime=30,
        cpus_per_task=10,
    conda:
        # Envoirnment with tidyverse and data.table
    log:
        "logs/blupf90_pre_processing_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1
        Rscript blupf90_preprocessing.R {wildcards.n} > {log} 2>&1 && touch {output}
        """

# Runs renumf90 in BLUP_PED and BLUP_GENO_ASC, then preGSf90 once in BLUP_GENO_ASC.
# Copies GimA22i and genotype files to BLUP_GENO_CMC so it can reuse the precomputed G inverse.
# Reruns renumf90 in BLUP_GENO_ASC with the expanded model par file (maternal effects).
# Runs renumf90 in BLUP_GENO_CMC with the CMC + PE model par file.
# All downstream Gibbs rules depend on this rule being done.
rule renumf90_run:
    output:
        "results/renumf90_run_{n}_done.txt"
    input:
        "results/blupf90_pre_processing_run_{n}_done.txt"
    resources:
        mem_mb=100000,
        disk_mb=50000,
        runtime=45,
        cpus_per_task=10,
    conda:
        # envoirnment with blupf90
    log:
        "logs/renumf90_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1

        # --- BLUP_PED ---
        echo "[INFO] Running renumf90 in BLUP_PED for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] BLUP_PED renumf90 start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_ped && renumf90 renumf90.par) >> {log} 2>&1
        echo "[TIME] BLUP_PED renumf90 end: $(date)" >> {log} 2>&1

        # --- BLUP_GENO_ASC: copy and subsample genotype files ---
        echo "[INFO] Copying .map file to BLUP_GENO_ASC for run {wildcards.n}" >> {log} 2>&1
        cp results/run{wildcards.n}/*.map results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.map >> {log} 2>&1

        echo "[INFO] Shuffling and subsampling 30k animals from .ped for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] shuf start: $(date)" >> {log} 2>&1
        shuf -n 30000 results/run{wildcards.n}/*.ped -o results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.ped >> {log} 2>&1
        echo "[TIME] shuf end: $(date)" >> {log} 2>&1

        # --- BLUP_GENO_ASC: renumf90 with simple model (needed for preGSf90) ---
        echo "[INFO] Running renumf90 (pregibbs) in BLUP_GENO_ASC for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] BLUP_GENO_ASC renumf90 pregibbs start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_asc && renumf90 renumf90_pregibbs.par) >> {log} 2>&1
        echo "[TIME] BLUP_GENO_ASC renumf90 pregibbs end: $(date)" >> {log} 2>&1

        # --- BLUP_GENO_ASC: run preGSf90 once to build GimA22i ---
        echo "[INFO] Running preGSf90 in BLUP_GENO_ASC for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] preGSf90 start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_asc && preGSf90 renf90.par) >> {log} 2>&1
        echo "[TIME] preGSf90 end: $(date)" >> {log} 2>&1

        # --- Copy preGSf90 outputs to BLUP_GENO_CMC so it can reuse GimA22i ---
        echo "[INFO] Copying preGSf90 outputs to BLUP_GENO_CMC for run {wildcards.n}" >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/GimA22i results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.ped results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.map results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.ped_XrefID results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/freqdata.count.after.clean results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/freqdata.count results/run{wildcards.n}/blupf90_geno_cmc/ >> {log} 2>&1

        # --- Copy preGSf90 outputs to BLUP_COR_DIR so it can reuse GimA22i ---
        echo "[INFO] Copying preGSf90 outputs to BLUP_COR_DIR for run {wildcards.n}" >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/GimA22i results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.ped results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.map results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/sim_genotyped_30k.ped_XrefID results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/freqdata.count.after.clean results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1
        cp results/run{wildcards.n}/blupf90_geno_asc/freqdata.count results/run{wildcards.n}/blupf90_cor_asc_cmc/ >> {log} 2>&1

        # --- BLUP_GENO_ASC: rerun renumf90 with expanded model (maternal effects) ---
        # This produces the final renf90.par that gibbsf90+ will use for ASC
        echo "[INFO] Running renumf90 (expanded ASC model) in BLUP_GENO_ASC for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] BLUP_GENO_ASC renumf90 expanded start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_asc && renumf90 renumf90_gibbs_asc.par) >> {log} 2>&1
        echo "[TIME] BLUP_GENO_ASC renumf90 expanded end: $(date)" >> {log} 2>&1

        # --- BLUP_GENO_CMC: run renumf90 with CMC + PE model ---
        # This produces the renf90.par that gibbsf90+ will use for CMC
        echo "[INFO] Running renumf90 (CMC + PE model) in BLUP_GENO_CMC for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] BLUP_GENO_CMC renumf90 start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_cmc && renumf90 renumf90_gibbs_pregcmc.par) >> {log} 2>&1
        echo "[TIME] BLUP_GENO_CMC renumf90 end: $(date)" >> {log} 2>&1

        date >> {log} 2>&1
        touch {output}
        """

# Gibbs sampling for pedigree-based ASC model
rule gibbs_ped_run:
    output:
        "results/gibbs_ped_run_{n}_done.txt"
    input:
        "results/renumf90_run_{n}_done.txt"
    resources:
        mem_mb=10000,
        disk_mb=20000,
        runtime=6000, # time based on the number of samples. ca. 493 rounds/hour
        cpus_per_task=10,
    conda:
        # Envoirnment with BLUPF90.
    log:
        "logs/gibbs_ped_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_PED start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_ped && gibbsf90+ renf90.par --samples 5000 --burnin 2500 --interval 10) >> {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_PED end: $(date)" >> {log} 2>&1
        date >> {log} 2>&1
        touch {output}
        """

# Gibbs sampling for genomic expanded ASC model (direct + maternal genetic + maternal PE)
rule gibbs_geno_asc_run:
    output:
        "results/gibbs_geno_asc_run_{n}_done.txt"
    input:
        "results/renumf90_run_{n}_done.txt"
    resources:
        mem_mb=180000,
        disk_mb=20000,
        runtime=6000, # time based on the number of samples. ca. 160 rounds/hour
        cpus_per_task=10,
    conda:
        # Envoirnment with BLUPF90
    log:
        "logs/gibbs_geno_asc_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_GENO_ASC start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_asc && gibbsf90+ renf90.par --samples 2000 --burnin 1000 --interval 10) >> {log} 2>&1 
        echo "[TIME] gibbsf90+ BLUP_GENO_ASC end: $(date)" >> {log} 2>&1
        date >> {log} 2>&1
        touch {output}
        """

# Gibbs sampling for genomic CMC model (animal + permanent environment)
rule gibbs_geno_cmc_run:
    output:
        "results/gibbs_geno_cmc_run_{n}_done.txt"
    input:
        "results/renumf90_run_{n}_done.txt"
    resources:
        mem_mb=180000,
        disk_mb=20000,
        runtime=6000, # time based on number of sample. ca. 326 rounds/hour
        cpus_per_task=10,
    conda:
        # Conda envoirnment with BLUPF90
    log:
        "logs/gibbs_geno_cmc_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_GENO_CMC start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_geno_cmc && gibbsf90+ renf90.par --samples 4000 --burnin 2000 --interval 10) >> {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_GENO_CMC end: $(date)" >> {log} 2>&1
        date >> {log} 2>&1
        touch {output}
        """



# Runs renumf90 for the bivariate ASC-CMC correlation model, then runs gibbsf90+.
# Depends on renumf90_run being done so that GimA22i is already in BLUP_COR_DIR.
rule renumf90_cor_run:
    output:
        "results/renumf90_cor_run_{n}_done.txt"
    input:
        "results/renumf90_run_{n}_done.txt"
    resources:
        mem_mb=50000,
        disk_mb=50000,
        runtime=45,
        cpus_per_task=10,
    conda:
        # Conda envoirnment with BLUPF90
    log:
        "logs/renumf90_cor_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1

        echo "[INFO] Running renumf90 (correlation model) in BLUP_COR_DIR for run {wildcards.n}" >> {log} 2>&1
        echo "[TIME] renumf90 correlation start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_cor_asc_cmc && renumf90 renumf90_gibbs_correlation.par) >> {log} 2>&1
        echo "[TIME] renumf90 correlation end: $(date)" >> {log} 2>&1

        date >> {log} 2>&1
        touch {output}
        """

# Gibbs sampling for bivariate ASC-CMC correlation model
# This is the most computationally intensive run - requires many more samples for convergence
rule gibbs_cor_run:
    output:
        "results/gibbs_cor_run_{n}_done.txt"
    input:
        "results/renumf90_cor_run_{n}_done.txt"
    resources:
        mem_mb=200000,
        disk_mb=50000,
        runtime=7000, # time based on number of samples. ca. 78 rounds/hour
        cpus_per_task=10,
    conda:
        # Conda envoirnment with BLUPF90.
    log:
        "logs/gibbs_cor_run_{n}.log"
    threads: 10
    shell:
        """
        date > {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_COR start: $(date)" >> {log} 2>&1
        (cd results/run{wildcards.n}/blupf90_cor_asc_cmc && gibbsf90+ renf90.par --samples 3900 --burnin 1000 --interval 10) >> {log} 2>&1
        echo "[TIME] gibbsf90+ BLUP_COR end: $(date)" >> {log} 2>&1
        date >> {log} 2>&1
        touch {output}
        """

# Final sentinel rule: all three Gibbs runs must finish before the run is marked done
rule blupf90_run:
    output:
        "results/job_done_run_{n}.txt"
    input:
        "results/gibbs_ped_run_{n}_done.txt",
        "results/gibbs_geno_asc_run_{n}_done.txt",
        "results/gibbs_geno_cmc_run_{n}_done.txt",
        "results/gibbs_cor_run_{n}_done.txt"
    resources:
        mem_mb=1200,
        disk_mb=200,
        runtime=4,
        cpus_per_task=1,
    conda:
        "/home/giovincavallo/miniforge3/envs/alphasim"
    log:
        "logs/blupf90_run_{n}.log"
    threads: 1
    shell:
        """
        date > {log} 2>&1
        touch {output}
        """