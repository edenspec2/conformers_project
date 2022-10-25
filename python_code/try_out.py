# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 23:17:15 2022

@author: edens
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PyMol
from rdkit.Chem.Subshape import SubshapeAligner, SubshapeBuilder, SubshapeObjects
mols = [m for m in Chem.SDMolSupplier('cdk2.sdf')]
for m in mols:
    molid = m.GetProp('id')
    m.SetProp('_Name', molid) #_Name prop is required for align with shape-it

ref = Chem.Mol(mols[0].ToBinary())
probe = Chem.Mol(mols[1].ToBinary())
AllChem.CanonicalizeConformer(ref.GetConformer())
builder = SubshapeBuilder.SubshapeBuilder()
builder.gridDims = (20.,20.,10)
builder.gridSpacing=0.5
builder.winRad = 4.
refShape = builder.GenerateSubshapeShape(ref)
probeShape = builder.GenerateSubshapeShape(probe)
 
aligner = SubshapeAligner.SubshapeAligner()
algs = aligner.GetSubshapeAlignments(ref, refShape, probe, probeShape, builder)
 
alg = algs[0]
AllChem.TransformMol(probe, alg.transform)
newprobeShape = builder(probe)
v = PyMol.MolViewer()
v.ShowMol(ref, name='ref')
SubshapeObjects.DisplaySubshape(v, refShape, 'ref_Shape')
v.server.do('set transparency=0.5')
v.ShowMol(probe, name='probe', showOnly=False)
SubshapeObjects.DisplaySubshape(v, newprobeShape, 'prob_Shape')
v.server.do('set transparency=0.5')
v.GetPNG()