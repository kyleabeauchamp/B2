import pandas as pd
from density_simulation_parameters import DATA_PATH, MOLECULES_PER_BOX

# You should have already run count_classes.py to generate this file, which lists the experiments to run.
EXPERIMENTS = pd.read_csv("./inputs/data_dielectric.csv")
EXPERIMENTS["IDENTIFIER"] = EXPERIMENTS["cas"] + "_" + EXPERIMENTS["Temperature, K"].astype('str')
EXPERIMENTS = EXPERIMENTS.set_index("IDENTIFIER")
EXPERIMENTS["MOLECULES"] = MOLECULES_PER_BOX
EXPERIMENTS = EXPERIMENTS.T.to_dict()

IDENTIFIERS = list(EXPERIMENTS.keys())
ALL_CAS = {v["cas"] for k, v in EXPERIMENTS.items()}

rule all:
    input:
        monomer_mol2 = expand("results/monomers/{cas}.mol2", cas=ALL_CAS),
        equil_dcd = expand("results/equil/{identifier}.dcd", identifier=IDENTIFIERS),
        production_dcd = expand("results/production/{identifier}.dcd", identifier=IDENTIFIERS),

        #si_csv_in = "results/tables/full_filtered_data.csv",
        #si_csv = "results/tables/data_with_metadata.csv",

        #latex = "results/tables/functional_groups_latex.txt",
        predictions_csv = expand("results/tables/predictions/{identifier}.csv", identifier=IDENTIFIERS),
        all_predictions_csv = "results/tables/predictions.csv",


rule build_monomer:
    output:
        pdb = "results/monomers/{cas}.pdb",
        mol2 = "results/monomers/{cas}.mol2",
        frcmod = "results/monomers/{cas}.frcmod",
    shell:
        "code/lbrun.py build_monomer "
        "--cas=\"'{wildcards.cas}'\" "  # NB: need extra quotes in fire!
        "--pdb={output.pdb} --mol2={output.mol2} --frcmod={output.frcmod};"


rule build_box:
    params:
        cas = lambda wildcards: EXPERIMENTS[wildcards.identifier]["cas"],
        molecules = lambda wildcards: EXPERIMENTS[wildcards.identifier]["MOLECULES"],
    input:
        mol2 = lambda wildcards: "results/monomers/{cas}.mol2".format(cas=EXPERIMENTS[wildcards.identifier]["cas"]),
        frcmod = lambda wildcards: "results/monomers/{cas}.frcmod".format(cas=EXPERIMENTS[wildcards.identifier]["cas"]),
    output:
        pdb = "results/packmol_boxes/{identifier}.pdb",
        prmtop = "results/tleap/{identifier}.prmtop",
        inpcrd = "results/tleap/{identifier}.inpcrd",
    shell:
        "code/lbrun.py build_box "
        "--in_mol2={input.mol2} "
        "--in_frcmod={input.frcmod} "
        "--out_pdb={output.pdb} "
        "--out_inpcrd={output.inpcrd} "
        "--out_prmtop={output.prmtop} "
        "--n_monomers={params.molecules};"


rule equilibrate:
    params:
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        prmtop = "results/tleap/{identifier}.prmtop",
        inpcrd = "results/tleap/{identifier}.inpcrd",
    output:
        dcd = "results/equil/{identifier}.dcd",
        pdb = "results/equil/{identifier}.pdb",
    shell:
        "code/lbrun.py equilibrate "
        "--in_prmtop={input.prmtop} "
        "--in_inpcrd={input.inpcrd} "
        "--out_dcd={output.dcd} "
        "--out_pdb={output.pdb} "
        "--temperature={params.temperature};"


rule production:
    params:
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        prmtop = "results/tleap/{identifier}.prmtop",
        pdb = "results/equil/{identifier}.pdb",
    output:
        dcd = "results/production/{identifier}.dcd",
        pdb = "results/production/{identifier}.pdb",
        csv = "results/production/{identifier}.csv",
    shell:
        "code/lbrun.py production "
        "--in_prmtop {input.prmtop} "
        "--in_pdb {input.pdb} "
        "--out_dcd {output.dcd} "
        "--out_csv {output.csv} "
        "--temperature {params.temperature};"

rule build_si_table:
    input:
        incsv = "results/tables/full_filtered_data.csv",
    output:
        outcsv = "results/tables/data_with_metadata.csv",
    shell:
        "code/create_data_table_for_si.py {input.csv} {output.outcsv};"


rule list_functional_groups:
    input:
        incsv = "results/tables/data_dielectric.csv",
    output:
        latex = "results/tables/functional_groups_latex.txt",
    shell:
        "code/list_functional_groups.py {input.incsv} > {output.latex};"


rule generate_prediction:
    params:
        cas = lambda wildcards: EXPERIMENTS[wildcards.identifier]["cas"],
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        prmtop = "results/tleap/{identifier}.prmtop",
        dcd = "results/production/{identifier}.dcd",
        pdb = "results/production/{identifier}.pdb",
        csv = "results/production/{identifier}.csv",
    output:
        csv = "results/tables/predictions/{identifier}.csv",
    shell:
        "code/munge_output_amber.py predict"
        "--in_prmtop={input.prmtop} "
        "--in_csv={input.csv} "
        "--in_dcd={input.dcd} "
        "--out_csv={output.csv} "
        "--cas=\"{params.cas}\" "
        "--temperature {params.temperature};"


rule merge_predictions:
    input:
        csv = expand("results/tables/predictions/{identifier}.csv", identifier=IDENTIFIERS),
    output:
        csv = "results/tables/predictions.csv",
    shell:
        "code/munge_output_amber.py merge"
        "--incsv={input.csv} --outcsv={output.csv};"