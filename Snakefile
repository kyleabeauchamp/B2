import pandas as pd
from density_simulation_parameters import MOLECULES_PER_BOX

# I'm assuming you have already generated these two CSV files, which lists the experiments to run.
FULLDATA_CSV = "inputs/full_filtered_data.csv",
INPUT_DIELECTRIC = "inputs/data_dielectric.csv"

EXPERIMENTS = pd.read_csv(INPUT_DIELECTRIC)
EXPERIMENTS["IDENTIFIER"] = EXPERIMENTS["cas"] + "_" + EXPERIMENTS["Temperature, K"].astype('str')
EXPERIMENTS = EXPERIMENTS.set_index("IDENTIFIER")
EXPERIMENTS["MOLECULES"] = MOLECULES_PER_BOX
EXPERIMENTS = EXPERIMENTS[0:5]  # Just do a couple for testing
EXPERIMENTS = EXPERIMENTS.T.to_dict()

BLACKLIST = {"126068-67-5_313.2", "126068-67-5_303.2", "110-97-4_313.2", "118662-30-9 298.2", "118662-30-9_298.2", "118662-30-9_293.2", "105-59-9_318.2", "105-59-9_293.2",
"105-59-9_308.2"}  # Some molecules don't converge enough to run bootstraps with the limited steps trial run , so exclude them.
EXPERIMENTS = {key: val for key, val in EXPERIMENTS.items() if key not in BLACKLIST}

IDENTIFIERS = list(EXPERIMENTS.keys())
ALL_CAS = {v["cas"] for k, v in EXPERIMENTS.items()}

#FFXML = None  # Use this for GAFF
FFXML = "forcefield/smirnoff99Frosst.ffxml"  # Use this line for smirff

rule all:
    input:
        monomer_mol2 = expand("results/monomers/{cas}.mol2", cas=ALL_CAS),
        equil_dcd = expand("results/equil/{identifier}.dcd", identifier=IDENTIFIERS),
        production_dcd = expand("results/production/{identifier}.dcd", identifier=IDENTIFIERS),

        si_csv = "results/tables/data_with_metadata.csv",

        # latex = "results/tables/functional_groups_latex.txt",  # NB: I haven't tested this yet, no OE installation right now...
        predictions_csv = expand("results/tables/predictions/{identifier}.csv", identifier=IDENTIFIERS),
        all_predictions_csv = "results/tables/predictions.csv",

        dens_pdf = "results/figures/densities_thermoml.pdf",

rule build_monomer:
    output:
        mol2 = "results/monomers/{cas}.mol2",
        frcmod = "results/monomers/{cas}.frcmod",
    shell:
        "code/lbrun.py build_monomer "
        "\"'{wildcards.cas}'\" "  # NB: need extra quotes in fire to not try to cast CAS string as integers
        "{output.mol2} "
        "{output.frcmod} ;"


rule build_box:
    params:
        cas = lambda wildcards: EXPERIMENTS[wildcards.identifier]["cas"],
        molecules = lambda wildcards: EXPERIMENTS[wildcards.identifier]["MOLECULES"],
        ff_flag = lambda wildcards: ("--ffxml=%s" % FFXML if FFXML is not None else ""),
    input:
        mol2 = lambda wildcards: "results/monomers/{cas}.mol2".format(cas=EXPERIMENTS[wildcards.identifier]["cas"]),
        frcmod = lambda wildcards: "results/monomers/{cas}.frcmod".format(cas=EXPERIMENTS[wildcards.identifier]["cas"]),
    output:
        pdb = "results/packmol_boxes/{identifier}.pdb",
        top = "results/tleap/{identifier}.top",
        gro = "results/tleap/{identifier}.gro",
    shell:
        "code/lbrun.py build_box "
        "--in_mol2={input.mol2} "
        "--in_frcmod={input.frcmod} "
        "--out_pdb={output.pdb} "
        "--out_gro={output.gro} "
        "--out_top={output.top} "
        "--n_monomers={params.molecules} "
        "{params.ff_flag};"


rule equilibrate:
    params:
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        top = rules.build_box.output.top,
        gro = rules.build_box.output.gro,
    output:
        dcd = "results/equil/{identifier}.dcd",
        pdb = "results/equil/{identifier}.pdb",
    shell:
        "code/lbrun.py equilibrate "
        "--in_top={input.top} "
        "--in_gro={input.gro} "
        "--out_dcd={output.dcd} "
        "--out_pdb={output.pdb} "
        "--temperature={params.temperature};"


rule production:
    params:
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        top = rules.build_box.output.top,
        pdb = rules.equilibrate.output.pdb,
    output:
        dcd = "results/production/{identifier}.dcd",
        csv = "results/production/{identifier}.csv",
    shell:
        "code/lbrun.py production "
        "--in_top {input.top} "
        "--in_pdb {input.pdb} "
        "--out_dcd {output.dcd} "
        "--out_csv {output.csv} "
        "--temperature {params.temperature};"

rule build_si_table:
    output:
        outcsv = "results/tables/data_with_metadata.csv",
    shell:
        "code/create_data_table_for_si.py {FULLDATA_CSV} {output.outcsv};"


rule list_functional_groups:
    output:
        latex = "results/tables/functional_groups_latex.txt",
    shell:
        "code/list_functional_groups.py {INPUT_DIELECTRIC} > {output.latex};"


rule generate_prediction:
    params:
        cas = lambda wildcards: EXPERIMENTS[wildcards.identifier]["cas"],
        temperature = lambda wildcards: EXPERIMENTS[wildcards.identifier]["Temperature, K"],
    input:
        top = rules.build_box.output.top,
        dcd = rules.production.output.dcd,
        pdb = rules.equilibrate.output.pdb,
        csv = rules.production.output.csv,
    output:
        csv = "results/tables/predictions/{identifier}.csv",
    shell:
        "code/munge_output_gromacs.py predict "
        "--in_top={input.top} "
        "--in_pdb {input.pdb} "
        "--in_csv={input.csv} "
        "--in_dcd={input.dcd} "
        "--out_csv={output.csv} "
        "--cas=\"'{params.cas}'\" "  # NB: need extra quotes in fire!
        "--temperature {params.temperature};"


rule merge_predictions:
    input:
        csv = expand("results/tables/predictions/{identifier}.csv", identifier=IDENTIFIERS),
    output:
        csv = "results/tables/predictions.csv",
    params:
        csv = ",".join(expand("results/tables/predictions/{identifier}.csv", identifier=IDENTIFIERS)),
    shell:
        "code/munge_output_gromacs.py merge "
        "--incsv='{params.csv}' --outcsv={output.csv};"


rule plot_tbv:
    input:
        pred_csv = rules.merge_predictions.output.csv,
        data_with_metadata = rules.build_si_table.output.outcsv,
    output:
        dens_pdf = "results/figures/densities_thermoml.pdf",
        diff_pdf = "results/figures/densities_differences_thermoml.pdf",
        diel_pdf = "results/figures/dielectrics_thermoml.pdf",
        nocorr_pdf = "results/figures/dielectrics_thermoml_nocorr.pdf",
    shell:
        "code/plot_tbv.py "
        "--expt_csv={input.data_with_metadata} "
        "--pred_csv={input.pred_csv} "
        "--dens_pdf={output.dens_pdf} "
        "--diff_pdf={output.diff_pdf} "
        "--diel_pdf={output.diel_pdf} "
        "--nocorr_pdf={output.nocorr_pdf} "
