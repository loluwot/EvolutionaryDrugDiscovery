from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from random import randint
import os
import numpy as np
import deepchem as dc
from deepchem.molnet import load_pdbbind_grid
protein_path = '3lpt.pdb'
possible_atoms = ['C', 'N', 'O', 'S', 'Cl', 'Br', 'F', 'I']
possible_bonds = [[1,2,3], [1,2,3], [1,2], [1,2], [1], [1], [1], [1]]
bond_type = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]
pdbbind_tasks, pdbbind_datasets, transformers = load_pdbbind_grid(split='random', subset='full')
train_data, valid_data, test_data = pdbbind_datasets
print(train_data.X[0].shape)
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
reload = dc.models.MultitaskRegressor(len(pdbbind_tasks), train_data.X.shape[1],dropouts=[.25],learning_rate=0.0003,weight_init_stddevs=[.1],batch_size=64, model_dir='pdbbind_test')
print(reload.get_checkpoints())
reload.load_from_pretrained(reload, model_dir='pdbbind_test')
reload.evaluate(train_data, [metric], transformers)
feature = dc.feat.RdkitGridFeaturizer(voxel_width=16.0,feature_types=["ecfp", "splif", "hbond", "salt_bridge"],ecfp_power=9,splif_power=9,flatten=True)
# scores = []
# for i in range(0,100):
	# grid = feature.featurize_complexes(['gen1/ligand' + str(i) + '.sdf'], ['3lpt.pdb'])
	# #print(grid[0].shape)
	# scores.extend(reload.predict_on_batch(grid[0]))
atom_dict = {
	'H':[1],
	'C':[1,2,3],
	'N':[1,2,3],
	'O':[1,2],
	'S':[1,2],
	'Cl':[1],
	'Br':[1],
	'F':[1],
	'I':[1]
}
atom_valences = {
	'H':[1],
	'C':[4],
	'N':[3],
	'O':[2],
	'S':[2,6],
	'Cl':[1],
	'Br':[1],
	'F':[1],
	'I':[1]
}
def init_gen (generation, final_size, n_children, hashed, n):
	while (len(generation) < final_size):
		rand_mol = generation[randint(0, len(generation)-1)]
		child = mutate(rand_mol)
		if (Chem.RDKFingerprint(child) not in hashed):
			hashed.append(Chem.RDKFingerprint(child))
			generation.append(child)
	count = 0
	os.mkdir('gen' + str(n))
	for m in generation:
		w = Chem.SDWriter('gen' + str(n)+ '/ligand' + str(count) + '.sdf')
		w.write(m)
		count = count + 1
	img=Draw.MolsToGridImage(generation,molsPerRow=4,subImgSize=(300,300))
	img.save('generation' + str(n) +'.png')
	return generation
			
def mutate(mol):
	rand_op = randint(0,3)
	if (rand_op == 0):
		return Chem.RemoveHs(addAtom(mol))
	elif (rand_op == 1):
		return Chem.RemoveHs(addBond(mol))
	# elif (rand_op == 2):
		# return Chem.RemoveHs(removeAtom(mol))
	elif (rand_op == 2):
		return Chem.RemoveHs(removeBond(mol))
	elif (rand_op == 3):
		return Chem.RemoveHs(replaceAtom(mol))
	
def score(n, gen_size, gen):
	new_gen = []
	scores = []
	for i in range(0,gen_size):
		grid = feature.featurize_complexes(['gen'+str(n)+'/ligand' + str(i) + '.sdf'], [protein_path])
		#print(grid[0].shape)
		scores.extend(reload.predict_on_batch(grid[0]))
	sorted = scores.copy()
	sorted.sort(reverse=True)
	#top_10_score = []
	for i in range(0, 10):
		new_gen.append(scores.index(sorted[i]))
		
	print(str(sorted[0]) + " maximum score------------")
	print(new_gen)
	next_gen = [gen[x] for x in new_gen]
	return next_gen
def addAtom(mol):	
	bondType = randint(1,3)
	possible_id = []
	for atom in mol.GetAtoms():
		#print(atom.GetSymbol() + "  atom")
		if (bondType in atom_dict[atom.GetSymbol()] and atom.GetImplicitValence() >= bondType):
			possible_id.append(atom.GetIdx())
	if (len(possible_id) > 0):
		#print(str(possible_id) + "possible")
		#print(bondType)
		possible = []
		for atom in possible_atoms:
			if (bondType in atom_dict[atom]):
				possible.append(atom)
		#print(possible)
		if (len(possible)>0):
			id1 = possible_id[randint(0, len(possible_id)-1)]
			a = Chem.rdchem.Atom(possible[randint(0, len(possible)-1)])
			#a = Chem.rdchem.Atom(12)
			editable = Chem.rdchem.EditableMol(mol)
			id2 = editable.AddAtom(a)
			editable.AddBond(id1, id2, bond_type[bondType-1])
			newMol = editable.GetMol()
			#Chem.SanitizeMol(newMol)
			newMol.UpdatePropertyCache()
			#for atom in newMol.GetAtoms():
				#print(str(atom.GetImplicitValence()) + " " + atom.GetSymbol())
			newMol = Chem.AddHs(newMol)
			AllChem.Compute2DCoords(newMol)
			#Draw.MolToFile(newMol,'test.png')
			AllChem.EmbedMolecule(newMol)
			#print(newMol.GetNumAtoms())
			return newMol
		return mol
	return mol
def addBond(mol):
	bondType = randint(1,3)
	possible_id = []
	for atom in mol.GetAtoms():
		if (bondType in atom_dict[atom.GetSymbol()] and atom.GetImplicitValence() >= bondType):
			possible_id.append(atom.GetIdx())
	#print(possible_id)
	if len(possible_id) > 0:
		index = randint(0, len(possible_id)-1)
		#print(bondType)
		if len(possible_id) > 1:	
			id1 = possible_id[index]
			del possible_id[index]
			id2 = possible_id[randint(0, len(possible_id)-1)]
			editable = Chem.rdchem.EditableMol(mol)
			if mol.GetBondBetweenAtoms(id1, id2) is not None:
				editable.RemoveBond(id1, id2)
			editable.AddBond(id1, id2, bond_type[bondType-1])
			newMol = editable.GetMol()
			newMol.UpdatePropertyCache()
			newMol = Chem.AddHs(newMol)
			AllChem.Compute2DCoords(newMol)
			AllChem.EmbedMolecule(newMol)
			Draw.MolToFile(newMol, 'test.png')
			return newMol
		return mol
	return mol
def removeAtom(mol):
	atoms = mol.GetAtoms()
	if (len(atoms) > 1):
		r_atom = atoms[randint(0, len(atoms)-1)]
		editable = Chem.rdchem.EditableMol(mol)
		editable.RemoveAtom(r_atom.GetIdx())
		newMol = editable.GetMol()
		newMol.UpdatePropertyCache()
		newMol = Chem.AddHs(newMol)
		AllChem.Compute2DCoords(newMol)
		frags = Chem.rdmolops.GetMolFrags(newMol, asMols=True)
		#print(str(frags) + "frags")
		largest_frag = frags[0]
		for mol_frag in frags:
			if (mol_frag.GetNumAtoms() >= largest_frag.GetNumAtoms()):
				largest_frag = mol_frag
		largest_frag.UpdatePropertyCache()
		AllChem.Compute2DCoords(largest_frag)		
		AllChem.EmbedMolecule(largest_frag)
		Draw.MolToFile(largest_frag, 'test.png')
		return largest_frag
	return mol
def removeBond(mol):
	Chem.Kekulize(mol, clearAromaticFlags=True)
	bonds = mol.GetBonds()
	if (len(bonds)>0):
		r_bond = bonds[randint(0, len(bonds)-1)]
		editable = Chem.rdchem.EditableMol(mol)
		editable.RemoveBond(r_bond.GetBeginAtomIdx(), r_bond.GetEndAtomIdx())
		newMol = editable.GetMol()
		newMol.UpdatePropertyCache()
		newMol = Chem.AddHs(newMol)
		AllChem.Compute2DCoords(newMol)
		Chem.Kekulize(newMol)
		frags = Chem.rdmolops.GetMolFrags(newMol,sanitizeFrags=True, asMols=True)
		largest_frag = frags[0]
		for mol_frag in frags:
			if (mol_frag.GetNumAtoms() >= largest_frag.GetNumAtoms()):
				largest_frag = mol_frag
		largest_frag.UpdatePropertyCache()
		AllChem.Compute2DCoords(largest_frag)	
		Chem.Kekulize(largest_frag)
		AllChem.EmbedMolecule(largest_frag)
		Draw.MolToFile(largest_frag, 'test.png')
		return largest_frag
	return mol
def replaceAtom(mol):
	Chem.Kekulize(mol, clearAromaticFlags=True)
	atoms = mol.GetAtoms()
	r_atom = atoms[randint(0, len(atoms)-1)]
	possible = []
	for s in possible_atoms:
		if (r_atom.GetTotalValence() in atom_valences[s]):
			possible.append(s)
	#print(possible)
	a = Chem.rdchem.Atom(possible[randint(0, len(possible)-1)])
	editable = Chem.rdchem.EditableMol(mol)
	editable.ReplaceAtom(r_atom.GetIdx(), a)
	newMol = editable.GetMol()
	newMol.UpdatePropertyCache()
	newMol = Chem.AddHs(newMol)
	AllChem.Compute2DCoords(newMol)
	Chem.Kekulize(newMol)
	AllChem.EmbedMolecule(newMol)
	Draw.MolToFile(newMol, 'test.png')
	return newMol
init_mol = Chem.MolFromSmiles("[nH]1c2c(c(c(CC(=O)O)c1=O)c1ccccc1)cc(Cl)cc2")
AllChem.EmbedMolecule(init_mol)
init_generation = [init_mol]
num_generations = 2;
final_generation_size = 40
hashed = [Chem.RDKFingerprint(init_mol)]
for i in range(1,30):
	init_generation = init_gen(init_generation, final_generation_size, 10, hashed,i)
	init_generation = score(i, final_generation_size, init_generation)
	
