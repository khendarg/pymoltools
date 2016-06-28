#!/usr/bin/env python2
"""
Copyleft 2016 by Kevin Hendargo, all rights reversed
"""
from __future__ import print_function
import pymol, Bio.PDB

def seqload(files):
	if type(files) is str: files = [files]
	strucs = []
	for fn in files:
		parseme = None
		if "pdb" in fn:
			parseme = Bio.PDB.PDBParser()
		else:
			parseme = Bio.PDB.MMCIFParser()
			parsemore = Bio.PDB.MMCIF2Dict.MMCIF2Dict(fn)
			for i in sorted(parsemore.keys()): print(i, "::=", parsemore[i])
			#printme = zip(parsemore["_pdbx_unobs_or_zero_occ_residues.auth_asym_id"], parsemore["_pdbx_unobs_or_zero_occ_residues.auth_comp_id"])
			#_pdbx_poly_seq_scheme.ndb_seq_num
			#_pdbx_poly_seq_scheme.mon_id
			#_pdbx_poly_seq_scheme.pdb_mod_id
			printme = zip(parsemore["_pdbx_poly_seq_scheme.ndb_seq_num", parsemore["_pdbx_poly_seq_scheme.mon_id"], parsemore["_pdbx_poly_seq_scheme.pdb_mod_id"])

			for i in printme: print(i)
		strucs.append(parseme.get_structure(fn[0:fn.index(".")], fn))

pymol.cmd.extend("seqload", seqload)
