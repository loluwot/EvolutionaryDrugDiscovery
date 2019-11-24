from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from random import randint
possible_atoms = ['C', 'N', 'O', 'S', 'Cl', 'Br', 'F', 'I', 'P']
possible_bonds = [[1,2,3], [1,2,3], [1,2], [1,2], [1], [1], [1], [1], [1,2]]
bond_type = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]
atom_dict = {
	'C':[1,2,3],
	'N':[1,2,3],
	'O':[1,2],
	'S':[1,2],
	'Cl':[1],
	'Br':[1],
	'F':[1],
	'I':[1],
	'P':[1,2]
}
atom_valences = {
	'C':[4],
	'N':[3],
	'O':[2],
	'S':[2,6],
	'Cl':[1],
	'Br':[1],
	'F':[1],
	'I':[1],
	'P':[3,5]
}
# def mutate(mol, generation, n):
	# for _ in range(n):
		
def addAtom(mol):	
	bondType = randint(1,3)
	possible_id = []
	for atom in mol.GetAtoms():
		print(atom.GetSymbol() + "  atom")
		if (bondType in atom_dict[atom.GetSymbol()] and atom.GetImplicitValence() >= bondType):
			possible_id.append(atom.GetIdx())
	print(possible_id)
	print(bondType)
	possible = []
	for atom in possible_atoms:
		if (bondType in atom_dict[atom]):
			possible.append(atom)
	print(possible)
	id1 = possible_id[randint(0, len(possible_id)-1)]
	a = Chem.rdchem.Atom(possible[randint(0, len(possible)-1)])
	#a = Chem.rdchem.Atom(12)
	editable = Chem.rdchem.EditableMol(mol)
	id2 = editable.AddAtom(a)
	editable.AddBond(id1, id2, bond_type[bondType-1])
	newMol = editable.GetMol()
	#Chem.SanitizeMol(newMol)
	newMol.UpdatePropertyCache()
	for atom in newMol.GetAtoms():
		print(str(atom.GetImplicitValence()) + " " + atom.GetSymbol())
	newMol = Chem.AddHs(newMol)
	AllChem.Compute2DCoords(newMol)
	#Draw.MolToFile(newMol,'test.png')
	AllChem.EmbedMolecule(newMol)
	print(newMol.GetNumAtoms())
	return newMol
	
def addBond(mol):
	bondType = randint(1,3)
	possible_id = []
	for atom in mol.GetAtoms():
		if (bondType in atom_dict[atom.GetSymbol()] and atom.GetImplicitValence() >= bondType):
			possible_id.append(atom.GetIdx())
	print(possible_id)
	index = randint(0, len(possible_id)-1)
	print(bondType)
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
def removeAtom(mol):
	atoms = mol.GetAtoms()
	r_atom = atoms[randint(0, len(atoms)-1)]
	editable = Chem.rdchem.EditableMol(mol)
	editable.RemoveAtom(r_atom.GetIdx())
	newMol = editable.GetMol()
	newMol.UpdatePropertyCache()
	newMol = Chem.AddHs(newMol)
	AllChem.Compute2DCoords(newMol)
	frags = Chem.rdmolops.GetMolFrags(newMol,sanitizeFrags=True, asMols=True)
	largest_frag = frags[0]
	for mol_frag in frags:
		if (mol_frag.GetNumAtoms() >= largest_frag.GetNumAtoms()):
			largest_frag = mol_frag
	largest_frag.UpdatePropertyCache()
	AllChem.Compute2DCoords(largest_frag)		
	AllChem.EmbedMolecule(largest_frag)
	Draw.MolToFile(largest_frag, 'test.png')
def removeBond(mol):
	bonds = mol.GetBonds()
	r_bond = bonds[randint(0, len(bonds)-1)]
	editable = Chem.rdchem.EditableMol(mol)
	editable.RemoveBond(r_bond.GetBeginAtomIdx(), r_bond.GetEndAtomIdx())
	newMol = editable.GetMol()
	newMol.UpdatePropertyCache()
	newMol = Chem.AddHs(newMol)
	AllChem.Compute2DCoords(newMol)
	frags = Chem.rdmolops.GetMolFrags(newMol,sanitizeFrags=True, asMols=True)
	largest_frag = frags[0]
	for mol_frag in frags:
		if (mol_frag.GetNumAtoms() >= largest_frag.GetNumAtoms()):
			largest_frag = mol_frag
	largest_frag.UpdatePropertyCache()
	AllChem.Compute2DCoords(largest_frag)		
	AllChem.EmbedMolecule(largest_frag)
	Draw.MolToFile(largest_frag, 'test.png')
def replaceAtom(mol):
	atoms = mol.GetAtoms()
	r_atom = atoms[randint(0, len(atoms)-1)]
	possible = []
	for s in possible_atoms:
		if (r_atom.GetTotalValence() in atom_valences[s]):
			possible.append(s)
	print(possible)
	a = Chem.rdchem.Atom(possible[randint(0, len(possible)-1)])
	editable = Chem.rdchem.EditableMol(mol)
	editable.ReplaceAtom(r_atom.GetIdx(), a)
	newMol = editable.GetMol()
	newMol.UpdatePropertyCache()
	newMol = Chem.AddHs(newMol)
	AllChem.Compute2DCoords(newMol)
	AllChem.EmbedMolecule(newMol)
	Draw.MolToFile(newMol, 'test.png')
init_mol = Chem.MolFromSmiles("COCOCCOC")
init_generation = [init_mol]
num_generations = 2;
final_generation_size = 100
#addAtom(init_mol)
replaceAtom(init_mol)
