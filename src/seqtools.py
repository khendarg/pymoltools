#!/usr/bin/env python2
"""
Copyleft 2016 by Kevin Hendargo, all rights reversed
"""
from __future__ import print_function
import pymol, Bio.PDB
import os

def xclip(selection="sele"):
	if selection not in pymol.cmd.get_names("all", enabled_only=1): selection = "all"

	pymol.cmd.save(".asdftmp.fa~", selection, format="fasta")

	os.system("xclip .asdftmp.fa~")
	os.remove(".asdftmp.fa~")

#def seqload(*files):
#    """
#DESCRIPTION
#
#    "seqload" loads sequence data from a list of files
#
#    """
#	strucs = []
#	for fn in files:
#		parseme = None
#		if "pdb" in fn:
#			parseme = Bio.PDB.PDBParser()
#		else:
#			parsemore = Bio.PDB.MMCIF2Dict.MMCIF2Dict(fn)
#
#pymol.cmd.extend("seqload", seqload)
pymol.cmd.extend("xclip", xclip)
