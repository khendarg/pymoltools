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
			#for i in sorted(parsemore.keys()): print(i, "::=", parsemore[i])
			printme = zip(\
#absolute
#parsemore["_pdbx_poly_seq_scheme.asym_id"], \
#relative/"caption" chain ID
parsemore["_pdbx_poly_seq_scheme.pdb_strand_id"], \
#parsemore["_pdbx_poly_seq_scheme.entity_id"], \
parsemore["_pdbx_poly_seq_scheme.mon_id"], \
parsemore["_pdbx_poly_seq_scheme.pdb_mon_id"], \
#parsemore["_pdbx_poly_seq_scheme.seq_id"], \
#parsemore["_pdbx_poly_seq_scheme.pdb_seq_num"], \
parsemore["_pdbx_poly_seq_scheme.ndb_seq_num"], \
parsemore["_pdbx_poly_seq_scheme.auth_seq_num"] \
)
			#s1 = ""
			s2 = ""
			for i in printme: print(i[1], i[2])
			#for i in printme: #print(i[1], i[2])
				#s1 += Bio.PDB.protein_letters_3to1[i[1]]
				if i[2] == "?": s2 += "X"
				else: s2 += Bio.PDB.protein_letters_3to1[i[2]]
			#print(""+s1)
			print(""+s2)
			#print("----"+"MSVAVETFGFFMSALGLLMLGLTLSNSYWRVSTVHGNVITTNTIFENLWYSCATDSLGVSNCWDFPSMLALSGYVQGCRALMITAILLGFLGLFLGMVGLRCTNVGNMDLSKKAKLLAIAGTLHILAGACGMVAISWYAVNITTDFFNPLYAGTKYELGPALYLGWSASLLSILGGICVFSTCCCSSKEEPATRAGLPYKPSTVVIPRATSDESDISFGKYGKNAYV")

		#strucs.append(parseme.get_structure(fn[0:fn.index(".")], fn))

pymol.cmd.extend("seqload", seqload)
