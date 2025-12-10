# Expect: --config cfg=cfg_1 (or cfg_2, ...)
cfg_id = config.get("cfg", "cfg_1")  # default to cfg_1 if not given

rule all:
    input:
        f"results/{cfg_id}/submodels/MODEL_ANALYSIS_DONE.txt"

rule train_submodels:
    threads: 100
    input:
        cfg = f"configs/{cfg_id}.yaml"
    output:
        f"results/{cfg_id}/submodels/SUBMODELS_DONE.txt"
    shell:
        """
        Rscript src/R/sub_model_train.R {input.cfg} {threads}
        """

rule train_sg:
    threads: 10
    input:
        cfg      = f"configs/{cfg_id}.yaml",
        sub_done = f"results/{cfg_id}/submodels/SUBMODELS_DONE.txt"
    output:
        f"results/{cfg_id}/submodels/SG_DONE.txt"
    shell:
        """
        Rscript src/R/sg_model_train.R {input.cfg} {threads}
        """

rule model_analysis:
    threads: 1
    input:
        cfg    = f"configs/{cfg_id}.yaml",
        sg_done = f"results/{cfg_id}/submodels/SG_DONE.txt"
    output:
        f"results/{cfg_id}/submodels/MODEL_ANALYSIS_DONE.txt"
    shell:
        """
        Rscript src/R/model_analysis.R {input.cfg} {threads}
        """
