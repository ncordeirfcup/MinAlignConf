import os, time
import sys
import gzip
import pickle
from itertools import combinations
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.PropertyMol import PropertyMol


def prep_input(fname, id_field_name, nconf, energy, rms, seed):
    input_format = 'smi' 
    for mol, mol_name, act, mol_id in read_input(fname, input_format, None,True, True):
        yield mol, mol_name, nconf, energy, rms, seed, act, mol_id

def sorted_confids(mol):
    sorted_list = []

    for conf in mol.GetConformers():
        #ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())
        if ff is None:
            print(Chem.MolToSmiles(mol))
        else:
            sorted_list.append((conf.GetId(), ff.CalcEnergy()))

    if sorted_list:
        sorted_list = sorted(sorted_list, key=lambda x: x[1])

    return sorted_list


def remove_confs(mol, energy, rms):
    e = []
    for conf in mol.GetConformers():
        #ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())

        if ff is None:
            print(Chem.MolToSmiles(mol))
            return
        e.append((conf.GetId(), ff.CalcEnergy()))
    e = sorted(e, key=lambda x: x[1])

    if not e:
        return

    kept_ids = [e[0][0]]
    remove_ids = []

    for item in e[1:]:
        if item[1] - e[0][1] <= energy:
            kept_ids.append(item[0])
        else:
            remove_ids.append(item[0])

    if rms is not None:
        rms_list = [(i1, i2, AllChem.GetConformerRMS(mol, i1, i2)) for i1, i2 in combinations(kept_ids, 2)]
        while any(item[2] < rms for item in rms_list):
            for item in rms_list:
                if item[2] < rms:
                    i = item[1]
                    remove_ids.append(i)
                    break
            rms_list = [item for item in rms_list if item[0] != i and item[1] != i]

    for cid in set(remove_ids):
        mol.RemoveConformer(cid)


def gen_confs(mol, mol_name, nconf, energy, rms, seed, act, mol_id):
    print("Hello")
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=nconf, maxAttempts=700, randomSeed=seed)

    if len(cids) == 0:
        confs = gen_confs_obabel(mol, nconf=nconf)
        for conf in confs:
            mol.AddConformer(conf.GetConformer())
        cids = list(range(0, len(confs)))

    for cid in cids:
        try:
            #AllChem.MMFFOptimizeMolecule(mol, confId=cid)
            AllChem.UFFOptimizeMolecule(mol, confId=cid)
        except:
            continue
    remove_confs(mol, energy, rms)
    return mol_name, mol, act, mol_id


def read_input(fname, input_format=None, id_field_name=None, sanitize=True, removeHs=True):
    """
    fname - is a file name, None if STDIN
    input_format - is a format of input data, cannot be None for STDIN
    id_field_name - name of the field containing molecule name, if None molecule title will be taken
    """
    if input_format is None:
        tmp = os.path.basename(fname).split('.')
        if tmp == 'gz':
            input_format = '.'.join(tmp[-2:])
        else:
            input_format = tmp[-1]
    input_format = input_format.lower()
    if fname is None:    # handle STDIN
        if input_format == 'sdf':
            suppl = __read_stdin_sdf(sanitize=sanitize, removeHs=removeHs)
        elif input_format == 'smi':
            suppl = __read_stdin_smiles(sanitize=sanitize)
        else:
            raise Exception("Input STDIN format '%s' is not supported. It can be only sdf, smi." % input_format)
    elif input_format in ("sdf", "sdf.gz"):
        print("ok")
        suppl = __read_sdf(os.path.abspath(fname), input_format, id_field_name, sanitize, removeHs)
    elif input_format == 'pkl':
        suppl = __read_pkl(os.path.abspath(fname))
    elif input_format in ('smi'):
        suppl = __read_smiles(os.path.abspath(fname), sanitize)

    else:
        raise Exception("Input file format '%s' is not supported. It can be only sdf, sdf.gz, smi, pkl." % input_format)
    for mol_tuple in suppl:
        yield mol_tuple


def __read_smiles(fname, sanitize=True):
    with open(fname) as f:
        for line in f:
            tmp = line.strip().split(',')

            if not tmp:
                continue
            mol = Chem.MolFromSmiles(tmp[0], sanitize=sanitize)
            if mol is not None:
                if len(tmp) > 1:
                    mol_title = tmp[1]
                else:
                    mol_title = Chem.MolToSmiles(mol, isomericSmiles=True)
                if len(tmp) > 2:
                    act = tmp[2]
                    if act.lower() == 'active':
                        act = 1
                    elif act.lower() == 'inactive':
                        act = 0
                else:
                    act = None

                if len(tmp) > 3:
                    mol_id = tmp[3]
                else:
                    mol_id = tmp[1]

                yield mol, mol_title, act, mol_id

            else:
                print('Error mol', line)


def process(fname,outname,nconf,energy,rms,seed, writer):
    #def prep_input(fname, id_field_name, nconf, energy, rms, seed):
    for mol, mol_name, nconf, energy, rms, seed, act, mol_id in prep_input(fname, None, nconf,energy,rms,seed):
        mol_name, mol, act, mol_id=gen_confs(mol, mol_name, nconf, energy, rms, seed, act, mol_id)
        ids_sorted = sorted_confids(mol)
        if not ids_sorted:
           pass
           print(Chem.MolToSmiles(mol), mol_name)

        mol.SetProp("_Name", mol_name)
        mol.SetProp("Act", str(act))
        mol.SetProp("Mol", mol_id)
    

        for confId, energ in ids_sorted:
            name = '{name}_{confId}'.format(name=mol_name, confId=confId)
            mol.SetProp("_Name", name)
            writer.write(mol, confId=confId)

