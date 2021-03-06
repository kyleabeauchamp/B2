from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from density_simulation_parameters import CUTOFF
import parmed


def build_mixture_prmtop(gaff_mol2_filenames, box_filename, out_top, out_gro, ffxml):
    """Analog of openmoltools.amber.build_mixture_prmtop which uses SMIRNOFF forcefield (from github.com/open-forcefield-group/smarty) to parameterize small molecules, rather than GAFF.

Parameters
----------
mol2_filenames : list(str)
    Filenames of GAFF flavored mol2 files.  Each must contain exactly
    ONE solute.
box_filename : str
    Filename of PDB containing an arbitrary box of the mol2 molecules.
out_top : str
    output GROMACS topology filename.  Should have suffix .top
out_gro : str
    output gro filename.  Should have suffix .gro
ffxml : str
    filename containing input SMIRNOFF FFXML file for use in parameterizing the system


Returns
-------
success : bool
    True or False as to success of operation

Notes
-----
This can be easily broken if there are missing, duplicated, or
inconsistent ligand residue names in your box, mol2, and frcmod files.
You can use mdtraj to edit the residue names with something like
this: trj.top.residue(0).name = "L1"
"""
    from openeye import oechem
    from openforcefield.typing.engines.smirnoff import ForceField, PME
    from openforcefield.utils import read_molecules, get_data_filename, generateTopologyFromOEMol

    # Read in molecules
    oemols = []
    for mol2file in set(gaff_mol2_filenames):
        mol = oechem.OEGraphMol()
        ifs = oechem.oemolistream(mol2file)
        flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
        ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
        oechem.OEReadMolecule(ifs, mol)
        oechem.OETriposAtomNames(mol)
        oemols.append(mol)

    # Read in PDB file to get topology
    pdb = app.PDBFile(box_filename)

    # Load forcefield
    ff = ForceField(ffxml)

    # Construct system; charging not needed as mol2 files already have charges here
    system = ff.createSystem(pdb.topology, oemols, nonbondedMethod=PME, nonbondedCutoff=CUTOFF)

    # Dump to AMBER format
    structure = parmed.openmm.topsystem.load_topology(pdb.topology, system, pdb.positions)
    structure.save(out_top, overwrite=True)
    structure.save(out_gro, overwrite=True)

    return True
