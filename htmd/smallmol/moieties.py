from rdkit import Chem
from htmd.smallmol.smallmol import SmallMol

AROMATIC = Chem.rdchem.BondType.AROMATIC
SINGLE = Chem.rdchem.BondType.SINGLE
DOUBLE = Chem.rdchem.BondType.DOUBLE
TRIPLE = Chem.rdchem.BondType.TRIPLE
SP3 = Chem.rdchem.HybridizationType.SP3
SP2 = Chem.rdchem.HybridizationType.SP2
SP =  Chem.rdchem.HybridizationType.SP

class MoietyFragmenter:
    _heteroatoms = ['O', 'S', 'N', 'P'] 
    _halogens = ['Cl', 'F', 'Br' ]
    
    _bondtypes = [AROMATIC, SINGLE, DOUBLE, TRIPLE]
    _colors = [(1,0,0), (0,1,0), (0, 0,1), (1,1,0), (1,0,1), (0,1,1)]

    def __init__(self, mol):
        class WrongInputType(Exception):
            pass

        # Process as Rdkit molecule
        if isinstance(mol, SmallMol):
            self._mol = mol
 
        else:
            raise WrongInputType("Not a valid SmallMol object. Provide a valid one")

        self.moieties = self.getMoieties()
        

    def _isHeteroAtom(self, atom):
        symbol = atom.GetSymbol()
        if symbol not in self._heteroatoms:
            return False
    
        return True

    def _isHalogenAtom(self, atom):
        symbol = atom.GetSymbol()
        if symbol not in self._halogens:
            return False
    
        return True

    def _isConnectedTo(self, atom,  atomtype, bondtype=None, returnatoms=False):
        # bontype: AROMATIC, SINGLE, DOUBLE, TRIPLE
        atomIdx = atom.GetIdx()
        bonds = atom.GetBonds()
    
        atomMatched = []
        isconnect = False
    
        for b in bonds:
            atomOther = b.GetOtherAtom(atom)
            if b.GetBondType() != bondtype and bondtype != None: 
                continue
        
            if atomOther.GetSymbol() == atomtype:
                atomMatched.append(atomOther)
                isconnect = True
    
        if returnatoms:
            return isconnect, atomMatched
        return isconnect
    

    def getCarbonConnected(self, atom, hybridization=None,):
        '''hybridization: None--> all, sp2, sp1, sp3
        '''
        bondtype = SINGLE if hybridization == 'sp3' else DOUBLE if hybridization == 'sp2' else TRIPLE if hybridization == 'sp' else None

        isconnect, carbonatoms =  self._isConnectedTo(atom, 'C', bondtype, returnatoms=True)
        return carbonatoms
            
    def _getAcetalFg(self, heteroatoms):
        from collections import Counter
        carbon_sp3 = [ c for h_a in heteroatoms for c in  self.getCarbonConnected(h_a, hybridization='sp3')]
        carbonIdx_sp3 = [c.GetIdx() for c in carbon_sp3]

        carbon_acetals = [ self._mol._mol.GetAtomWithIdx(k) for k, v in Counter(carbonIdx_sp3).items() if v >= 2 ]
        hetero_acetals = [ [atom for atom in c.GetNeighbors() if self._isHeteroAtom(atom) ] for c in carbon_acetals  ]
        #print(carbon_acetals, [ c.GetIdx() for c in carbon_acetals])
        #print(hetero_acetals, [ [h_a.GetSymbol() for h_a in atoms] for atoms in hetero_acetals  ])
        
        return [ Moiety([c] + h_as) for c, h_as in zip(carbon_acetals, hetero_acetals) ]
        # return [ Moiety(atoms) for atoms in fgs ]
        
    def _getThreeTermRing(self, heteroatoms):
        hetero_inThreeRings = [h_a for h_a in heteroatoms if h_a.IsInRingSize(3) ]
        fgs = [ [h_a] +  [ a for a in h_a.GetNeighbors() if a.IsInRingSize(3)] for h_a in hetero_inThreeRings  ]
        return [ Moiety(atoms) for atoms in fgs ]

    def _getAlchenes(self, carbonatoms):
        import itertools
        carbonsSp2 = [ c for c in carbonatoms if c.GetHybridization() == SP2 and not c.GetIsAromatic() ]
        carbonsIdxs_Sp2 = [ c.GetIdx() for c in carbonsSp2 ]
        putative_bonded = list(itertools.combinations(carbonsIdxs_Sp2, 2))
                
        
        bonded_atoms = [ tuple(sorted(c_p)) for c_p in self._mol.get_bonds()[0].tolist() ]
   
        fgs = [[self._mol._mol.GetAtomWithIdx(c) for c in c_p] for c_p in putative_bonded if c_p in bonded_atoms ]
        return [ Moiety(atoms) for atoms in fgs ]

    def _getAlchines(self, carbonatoms):
        import itertools# first FG
        carbonsSp = [ c for c in carbonatoms if c.GetHybridization() == SP and not c.GetIsAromatic() ]
        carbonsIdxs_Sp = [ c.GetIdx() for c in carbonsSp ]
        putative_bonded = list(itertools.combinations(carbonsIdxs_Sp, 2))
        
        bonded_atoms = [ tuple(sorted(c_p)) for c_p in self._mol.get_bonds()[0].tolist() ]
   
        fgs = [[self._mol._mol.GetAtomWithIdx(c) for c in c_p] for c_p in putative_bonded if c_p in bonded_atoms ]
        return [ Moiety(atoms) for atoms in fgs ]

        
    def getMoieties(self):
        smallmol = self._mol 

        atoms = smallmol._mol.GetAtoms()

        fgs = []

        # heteroatoms as rdkit atom object
        heteroatoms = [ a for  a in  atoms if self._isHeteroAtom(a) ]

        # carbonatoms as rdkit atom object
        carbonatoms = [a for a in atoms if a.GetSymbol() == 'C']
        # Init Moieties
        # the moieties can be of three types         
        #   1) containing hetero atoms
        #   2) containing halogen atoms
        #   3) without hetero atoms

        # 1) get moieties with hetero atoms 
        # growing Fg with carbons sp2 sp3
        fg_hetero_carbonsSp2_carbonsSp = [[h_a] + self.getCarbonConnected(h_a, hybridization='sp2') + self.getCarbonConnected(h_a, hybridization='sp') for h_a in heteroatoms ]
        fg_hetero_carbonsSp2_carbonsSp = [Moiety(atoms) for atoms in fg_hetero_carbonsSp2_carbonsSp if len(atoms) != 1] 
        
        if len(fg_hetero_carbonsSp2_carbonsSp) != 0:
            # first FG
            fgs = fgs + fg_hetero_carbonsSp2_carbonsSp

        # merge Fg if heteroatoms part of acetal
        hetero_others = []
        if len(heteroatoms) > len(fg_hetero_carbonsSp2_carbonsSp):
            if len(fg_hetero_carbonsSp2_carbonsSp) == 0:
                hetero_others = [h_a for h_a in heteroatoms ]
            else:
                heteroIdxs_all = [(h_a.GetIdx()) for h_a in heteroatoms]
                heteroIdxs_matched = [a.GetIdx() for fg in fgs for a in fg.atoms if self._isHeteroAtom(a) ] 
                hetero_others = [self._mol._mol.GetAtomWithIdx(h_idx) for h_idx in list(set(heteroIdxs_all) - set(heteroIdxs_matched)) ]
                
        if len(hetero_others) >= 2:
            # second FG
            fg_hetero_acetal = self._getAcetalFg(hetero_others) 
            if len(fg_hetero_acetal) != 0:
                fgs = fgs + fg_hetero_acetal

        # check for 3 term rings with N, O, S
        if len(fgs) == 0:
            hetero_others = heteroatoms
        else:
            hetero_others = [h_a for h_a in heteroatoms for fg in fgs if h_a not in fg.atoms]
        if len(hetero_others) != 0:
            # third FG
            fg_hetero_threeTermRings = self._getThreeTermRing(hetero_others) 
            if len(fg_hetero_threeTermRings) != 0:
                fgs = fgs + fg_hetero_threeTermRings
        
        # all others heteroatoms not in previous criteria
        heteroIdxs_all = [(h_a.GetIdx()) for h_a in heteroatoms]
        heteroIdxs_matched = [a.GetIdx() for fg in fgs for a in fg.atoms if self._isHeteroAtom(a) ] 
        atoms_hetero_others = [[self._mol._mol.GetAtomWithIdx(h_idx)] for h_idx in list(set(heteroIdxs_all) - set(heteroIdxs_matched)) ]
        if len(atoms_hetero_others) != 0:
            # fourth FG
            fgs = fgs + [ Moiety(atom) for atom in atoms_hetero_others ]

        # 2) halogens
        halogenatoms = [ a for  a in  atoms if self._isHalogenAtom(a) ]
        if len(halogenatoms) != 0:
            # fifth FG
            fgs = fgs + [ Moiety(atom) for atom in halogenatoms ]

        # 3) without heteroatoms
        fg_alchenes = self._getAlchenes(carbonatoms)
        if len(fg_alchenes ) != 0:
            # sixth FG
            fgs = fgs + fg_alchenes

        fg_alchines = self._getAlchines(carbonatoms)
        if len(fg_alchines ) != 0:
            # seventh FG
            fgs = fgs + fg_alchines

        fgs = self._merge(fgs)

        for fg in fgs:
            fg.completeFg(self._mol._mol)

        return fgs

    def _merge(self, fgs):
        fgs_merged = []
        
        completed = True
        for fg in fgs:
            if fg in fgs_merged:
                continue
            ismerged = False
            for fg_merged in fgs_merged:
                if fg_merged.isMerged(fg):
                    ismerged = True
            if ismerged:
                continue
    
            isConnected = False
            for fg_merged in fgs_merged:
                isConnected = fg_merged.isConnected(fg)
                if isConnected:
                    completed = False
                    fg_merged.merge(fg)
                    break

            if not isConnected: 
                fgs_merged.append(fg)
                
        if not completed:
            self._merge(fgs_merged)

        return fgs_merged

    def depictFGs(self, fgs=None, filename=None, ipython=False, optimize=False, atomnames=False):
        from rdkit.Chem.Draw import IPythonConsole
        from rdkit.Chem.Draw import MolToImage
        from rdkit.Chem.Draw import rdMolDraw2D
        from IPython.display import SVG
        from rdkit.Chem.AllChem import EmbedMolecule

        if fgs == None:
            fgs = self.moieties

        drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        opts = drawer.drawOptions()
        if atomnames:
            for i in range(self._mol._mol.GetNumAtoms()):
                opts.atomLabels[i] = self._mol._mol.GetAtomWithIdx(i).GetSymbol()+str(i)
        
        highlightAtoms = [a for fg in fgs for a in fg.AtomsIdx ] + [a for fg in fgs for a in fg.EnviromentsIdx ]
        highlightColors = { a : self._colors[i%len(self._colors)] for i in range(len(fgs)) for a in fgs[i].AtomsIdx }

        if optimize:
            EmbedMolecule(self._mol._mol)
        

        drawer.DrawMolecule(self._mol._mol, highlightAtoms=highlightAtoms, highlightBonds=[], highlightAtomColors=highlightColors)
        
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        if filename != None:
            f = open(filename, 'w')
            f.write(svg)
            f.close()
        
        if ipython:
            svg = svg.replace('svg:', '')
            return SVG(svg)
        else:
            return None

class Moiety:
    _heteroatoms = ['O', 'S', 'N', 'P'] 
    _halogens = ['Cl', 'F', 'Br' ]
    
    _bondtypes = [AROMATIC, SINGLE, DOUBLE, TRIPLE]
    _colors = [(1,0,0), (0,1,0), (0, 0,1), (1,1,0), (1,0,1), (0,1,1)]
    
    def __init__(self, atoms):
        
        atoms = [atoms] if not isinstance(atoms, list) else atoms
        self.atoms = atoms
        self.bonds = self._getBonds(atoms)
        self.enviroments = []
        self._mol = None

        self.bondsenvironments = []

        self.name = None
        self.subcategory = None
        
    def _getBonds(self, atoms, exclude=None):
        atoms_idx = [ a.GetIdx() for a in atoms]
        bonds_found = [] if exclude == None else [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in exclude ]
        Bonds = []
        for a in atoms:
            bonds = a.GetBonds()
            for b in bonds:
                a1 = b.GetBeginAtomIdx()
                a2 = b.GetEndAtomIdx()
                if (a1 in atoms_idx and a2 in atoms_idx) and (a1, a2) not in bonds_found:
                    bonds_found.append((a1, a2))
                    Bonds.append(b)

        return Bonds

    @property
    def AtomsIdx(self):
        return [ atom.GetIdx() for atom in self.atoms ]
        
    @property
    def EnviromentsIdx(self):
        return [ atom.GetIdx() for atom in self.enviroments ]
        
    @property
    def Enviroments(self):
        return self.enviroments

    @property
    def Atoms(self):
        return  self.atoms 

    @property
    def Bonds(self):
        return self.bonds

    @property
    def BondsIdx(self):
        return [bond.GetIdx() for bond in self.bonds]

    @property
    def BondsEnvironments(self):
        return self.bondsenvironments

    @property
    def BondsEnvironmentsIdx(self):
        return [b.GetIdx() for bond in self.bondsenvironments]


    def _hasAtomIdx(self, atomidx):
        if atomidx in self.AtomsIdx:
            return True
        return False

    def isMerged(self, moiety):
        moiety_atoms_idx = moiety.AtomsIdx
        for atom in moiety_atoms_idx:
            if not self._hasAtomIdx(atom):
                return False
        return True

    def isConnected(self, moiety):
        moiety_atoms_idx = [ a.GetIdx() for atom in moiety.Atoms for a in atom.GetNeighbors() ]
        isconnected = False
        for atom in moiety_atoms_idx:
            if atom in self.AtomsIdx:
                isconnected = True
        return isconnected

    def merge(self, moiety):
        self.atoms = self.atoms + moiety.atoms
        self.bonds = self._getBonds(self.atoms) 

    def completeFg(self, mol):
        self._getEnvironments()
        self._keepHydrogens()

        self._mol = self._getMoietyMol(mol)

        #self.name = self._getName()
        
    def _getMoietyMol(self, mol):
        envAtoms = self.Enviroments
        envAtoms_ids = self.EnviromentsIdx

        atoms = self.Atoms
        atoms_ids = self.AtomsIdx

        other_atoms_del = []

        for env_atomId in envAtoms_ids:
            env_atom = mol.GetAtomWithIdx(env_atomId)
            atoms_bonded = [ b.GetOtherAtomIdx(env_atomId) for b in env_atom.GetBonds() ]
            atoms_bonded_elements = [ mol.GetAtomWithIdx(a).GetSymbol() for a in atoms_bonded ]
            start_atoms_del = [ aId for aId, aEl in zip(atoms_bonded, atoms_bonded_elements) if aId not in atoms_ids and aEl != 'H' ] 
            for aId_del in start_atoms_del:
                other_atoms_del.extend(self._findAtomToDel(mol, aId_del, env_atomId))

        other_atoms_del = list(set(other_atoms_del))
        new_mol = self.trim_atoms(mol, envAtoms_ids, other_atoms_del)

        return new_mol

    def trim_atoms(self,  mol, atomIds, atomIds_delete):
        rwmol = Chem.RWMol(mol)
        for atomId in atomIds:
            atom = rwmol.GetAtomWithIdx(atomId)
            #atom.SetAtomicNum(0)
            atom.SetIsotope(0)
            atom.SetFormalCharge(0)
            atom.SetIsAromatic(False)
            atom.SetNumExplicitHs(0)
        
        for aId in sorted(atomIds_delete, reverse=True):
            rwmol.RemoveAtom(aId)
            
        return rwmol.GetMol()

    def _get_atoms_to_visit(self, atom, seen_ids):
        neighbor_ids = []
        neighbor_atoms = []
        #print(atom.GetIdx())
        for bond in atom.GetBonds():
            neighbor_atom = bond.GetOtherAtom(atom)
            #print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            neighbor_id = neighbor_atom.GetIdx()
            #print('Neighbors: ',  neighbor_id)
            if neighbor_id not in seen_ids:
                neighbor_ids.append(neighbor_id) 
                neighbor_atoms.append(neighbor_atom)
            #print(neighbor_ids)
            #print()
        
        return neighbor_ids, neighbor_atoms

    def _findAtomToDel(self, mol, start_atomId, ignore_atomId):
        seen_ids = {start_atomId, ignore_atomId}
        atom_ids = [start_atomId]
        stack = [mol.GetAtomWithIdx(start_atomId)]
    
        while stack:
            atom = stack.pop()
            atom_ids_to_visit, atoms_to_visit = self._get_atoms_to_visit(atom, seen_ids)
            atom_ids.extend(atom_ids_to_visit)
            stack.extend(atoms_to_visit)
            seen_ids.update(atom_ids_to_visit)
        return atom_ids


    def _keepHydrogens(self):

        atoms = [a for a in self.Atoms if (a.GetSymbol() == 'C' and a.GetHybridization() == SP2) or a.GetSymbol() in self._heteroatoms ]
        bond_double = [ b for b in self.Bonds if b.GetBondType() == DOUBLE ]

        hydrogens_alkene = []
        if len(bond_double) != 0:
            carbon_atoms = [[b.GetBeginAtom(), b.GetEndAtom() ] for b in bond_double if [b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()] == ['C', 'C'] ]
            hydrogens_alkene = [ a.GetIdx() for b in carbon_atoms for c in b for a in c.GetNeighbors() if a.GetSymbol() == 'H' ]
        
        # Doubt: hydrogens on alkene. Useful for cys/trans. I prefer clean it
        hydrogens = [ atom for a in  atoms for atom in a.GetNeighbors() if atom.GetSymbol() == 'H' if atom.GetIdx() not in hydrogens_alkene]
        self.atoms = self.atoms + [ h for h in hydrogens if not self._hasAtomIdx(h.GetIdx()) ]
        self.bonds = self._getBonds(self.atoms)

        

    def _getName(self):
        from collections import Counter

        elements = self.get_elements()
        atoms = self.Atoms
        bonds = self.Bonds
        counts_atoms = Counter(elements) 

        # print("Naming")
        # print(">>> ", counts_atoms)

        atoms =  self.get_elements()
        if set(['O', 'N', 'S']) < set(atoms):
            # print("Moiety with S, N, O")
            name, subcategory = self._getNameSNO()

    def _getNameSNO(self):

        # sulfunamide
        name, subcategory = ('sulfunamide', self._getSubcategory('N')) if self._isSulfunamide() else (None, None)

        return name, subcategory

    def _isSulfunamide(self):
        from collections import Counter 
        _ref = [ ('N', SINGLE),  ('O', DOUBLE), ('O', DOUBLE),  ]
        
        elements = self.get_elements()
        atoms = self.Atoms
        bonds = self.Bonds
        enviroments = self.Enviroments
        bondsenvironments = self.BondsEnvironments
        counts_atoms = Counter(elements)

        if counts_atoms['O'] != 2 or counts_atoms['N'] != 1 or counts_atoms['S'] != 1:
            return False
        
        atom_sulfur = atoms[elements.index('S')]
        atoms_n, bonds_n = self._getNeighbors(atom_sulfur)
        
        atoms_bonds_pair = list(zip(atoms_n, bonds_n))

        return Counter(atoms_bonds_pair) == Counter(_ref)
        

    def _getSubcategory(self, atom):
        print("DAJE")

    def _getNeighbors(self, atom):
        import itertools
        atom_idx = atom.GetIdx()
        atoms_bonded = [ a for a in atom.GetNeighbors() ]
        atomIdx_bonded = [ a.GetIdx() for a in atoms_bonded ]
        atomsElement_bonded = [ a.GetSymbol() for a in atoms_bonded]

        atoms_pair = list(itertools.product([atom_idx], atomIdx_bonded))
        atoms_pair = [ sorted( list(p) ) for p in atoms_pair ] 
        
        _atoms = []
        _bonds = []
        for b in self.Bonds:
            atomsIdx_inBond = sorted([b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx()])
            if atomsIdx_inBond in atoms_pair:
                _other_atomIdx = [aIdx for aIdx in atomsIdx_inBond if aIdx != atom_idx][0]
                _idx = atomIdx_bonded.index(_other_atomIdx)
                _other_atom = atomsElement_bonded[_idx]
                _atoms.append(_other_atom)
                _bonds.append(b.GetBondType())

        return _atoms, _bonds        

    def _getEnvironments(self):
        carbons = [atom for a in self.Atoms for atom in a.GetNeighbors() if atom.GetSymbol() == 'C']
        carbons_env = [c for c in carbons if c.GetIdx() not in self.AtomsIdx]
        
        self.enviroments = carbons_env
        self.bondsenvironments = self._getBonds(self.enviroments + self.atoms, exclude=self.Bonds)

    def get_elements(self):
        elements = [ atom.GetSymbol() for atom in self.Atoms ]
        return elements          

    def depict(self, filename=None, ipython=False):
        from rdkit.Chem.Draw import IPythonConsole
        from rdkit.Chem.Draw import MolToImage
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem.AllChem import EmbedMolecule
        from IPython.display import SVG
        from rdkit.Chem import RWMol, MolFromSmiles, Atom, BondType, ChiralType

        _ = MolFromSmiles('C')
        rmol = RWMol(_)

        dict_old_new_idx = {}
        n = 1
        for  a in self.atoms:
            old_idx = a.GetIdx()
            rmol.AddAtom(a)
            dict_old_new_idx[old_idx] = n 
            n+=1
        
        for a in self.enviroments:
            old_idx = a.GetIdx()
            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
            a.SetIsAromatic(0)
            rmol.AddAtom(a)
            dict_old_new_idx[old_idx] = n 
            n+=1

        for b in self.Bonds:
             rmol.AddBond(dict_old_new_idx[b.GetBeginAtomIdx()], dict_old_new_idx[b.GetEndAtomIdx()], b.GetBondType())
        for b in self.bondsenvironments:
            rmol.AddBond(dict_old_new_idx[b.GetBeginAtomIdx()], dict_old_new_idx[b.GetEndAtomIdx()], b.GetBondType())
        
        rmol.RemoveAtom(0)

        rmol.UpdatePropertyCache(strict=False)
        EmbedMolecule(rmol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        
        drawer.DrawMolecule(rmol)
        
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        if filename != None:
            f = open(filename, 'w')
            f.write(svg)
            f.close()
        
        if ipython:
            svg = svg.replace('svg:', '')
            return SVG(svg)
        else:
            return None
