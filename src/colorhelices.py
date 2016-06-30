#!/usr/bin/env python2
"""
Copyleft 2016 by Kevin Hendargo, all rights reversed
"""
from __future__ import print_function

import pymol
import os 
import re
import Bio.PDB
import subprocess 
import math

def dictcat(x, indentcount=0):
    dlm = "    "
    if type(x) is dict:
        for i in sorted(x.keys()):
            if type(x[i]) is dict: 
                print(dlm*(indentcount+1) + str(i))
                dictcat(x[i], indentcount + 1)
            else:
                print(dlm*indentcount + str(i))
                print(dlm*(indentcount+1) + str(x[i]))
    else: 
        print(dlm*indentcount + str(i))
        print(dlm*(indentcount+1) + str(x[i]))

def hsv2rgb(h, s, v):
    """
DESCRIPTION

    Converts HSV tuples to RGB string representations
    """
    hsv = (float(h), float(s), float(v))
    chroma = hsv[2] * hsv[1] / 100
    hueimage = hsv[0]/60
    x = chroma * (1 - abs(hueimage % 2 - 1))
    if 0 <= hueimage and hueimage < 1: rgbtuple = (chroma, x, 0)
    elif 1 <= hueimage and hueimage < 2: rgbtuple = (x, chroma, 0)
    elif 2 <= hueimage and hueimage < 3: rgbtuple = (0, chroma, x)
    elif 3 <= hueimage and hueimage < 4: rgbtuple = (0, x, chroma)
    elif 4 <= hueimage and hueimage < 5: rgbtuple = (x, 0, chroma)
    elif 5 <= hueimage and hueimage < 6: rgbtuple = (chroma, 0, x)
    else: rgbtuple = (0,0,0)
    m = hsv[2] - chroma
    rgb = (int((rgbtuple[0] + m)*2.555), int((rgbtuple[1] + m)*2.555), int((rgbtuple[2] + m)*2.555))

    return "0x%02.x%02.x%02.x" % rgb 

def gradient(n, hue_start=0, hue_end=240, s=90, v=90):
    """
DESCRIPTION

    Returns a list of RGB hexes for a continuous hue gradient from a starting hue to an ending hue with constant saturation and value
    """
    n = int(n) - 1
    colors = []
    if n:
        for h in range(hue_start, hue_end + 1, int(math.ceil((hue_end - hue_start)/n))):
            colors.append(hsv2rgb(h, s, v))
        return colors
    else: 
        return [hsv2rgb((hue_start + hue_end)/2.0, s, v)]

def str2bool(s):
    try:
        if int(s): return True
        else: return False
    except ValueError:
        if "n" and "f" in s: return False
        else: return True

class Helix:
    def __init__(self, start, end, chain="A"):
        self.start = start
        self.end = end
        self.color = "0xff0000"
        self.chain = chain
    def __str__(self):
        return self.start + "-" + self.end

def stride(selection="sele"):
    """
DESCRIPTION

    Run Stride on a selection

USAGE

    stride [selection]

ARGUMENTS

    selection = string: selection of interest {default: sele}
    """
    #boilerplate for detecting whether a selection (usually sele) is enabled
    #runs stride on everything if no appropriate selection is available

    if selection not in pymol.cmd.get_names("all", True): selection = "(all)"

    #save state to avoid cluttering the user's workspace

    starting_names = set(pymol.cmd.get_names("all"))
    enabled_names = pymol.cmd.get_names("all", True)

    for o in enabled_names: pymol.cmd.split_chains(o)

    chains = {}

    #create a nice, hierarchical dict with all the chains...

    for c in pymol.cmd.get_names("all", True):
        try: chains[c[0:-2]][c[-1]] = []
        except KeyError: chains[c[0:-2]] = {c[-1]:[]}

    #...and iterate over it in an orderly fashion

    for o in sorted(chains.keys()):
        for c in sorted(chains[o].keys()): 

            #write a temporary pdb file with the coordinates of the current chain

            pymol.cmd.save("tmp.pdb", o + " and c. " + c)

            #run stride on the resulting pdb file and remove it

            struck = subprocess.check_output(["stride", "tmp.pdb"])
            os.remove("tmp.pdb")

            #assign more helices based on stride's assignments
            #TODO: check if other ss (even no ss) should be reassigned as well

            for l in re.split("\n", struck): 
                if "LOC" and "Helix" in l:
                    h = Helix(l[22:27].strip(), l[39:45].strip(), l[28])
                    try: chains[o][c].append(h)
                    except KeyError: chains[o][c] = [h]
                    #pymol.cmd.color("red", o + " and c. " + c + " and i. " + l[22:27].strip() + "-" + l[39:45].strip())

    #restore state

    ending_names = set(pymol.cmd.get_names("all"))

    for i in ending_names^starting_names: pymol.cmd.delete(i)
    for i in enabled_names: pymol.cmd.enable(i)

    #return the hierarchical dict containing chains and now assigned helices

    return chains

def hmmtop(selection="sele"):
    """
DESCRIPTION

    hmmtop runs hmmtop on loaded sequences

NOTES

    Many if not most loaded sequences will be incomplete relative to gene product sequences. Use with caution.
    """
    if selection not in pymol.cmd.get_names("all", True): selection = "(all)"

    starting_names = pymol.cmd.get_names("all")
    enabled_names = pymol.cmd.get_names("all", True)

    pymol.cmd.split_chains(selection)

    helices = {}

    for i in pymol.cmd.get_names("all", True):

        p = subprocess.Popen(["hmmtop", "-if=--", "-of=--"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        tmss, junk = p.communicate(pymol.cmd.get_fastastr(i))

        tmss = tmss.split()

        o = tmss[2][:-2]
        c = tmss[2][-1]
        try: helices[o][c] = []
        except KeyError: helices[o] = {c:[]}

        for i in range(5, len(tmss) - 1, 2):
            helices[o][c].append(Helix(tmss[i], tmss[i+1], c))

    ending_names = pymol.cmd.get_names("all")
    for i in list(set(ending_names)^set(starting_names)):
        pymol.cmd.delete(i)
    for i in enabled_names:
        pymol.cmd.enable(i)
    return helices
    
pymol.cmd.extend("hmmtop", hmmtop)

def get_helices(coordfile):

    helices = {}
    if "pdb" in coordfile:
        fh = open(coordfile)
        for l in fh:
            if "HELIX" in l[0:6]:
                #helices.append(Helix(l[21:25].strip(), l[33:37].strip(), l[19]))
                try: helices[l[19]].append(Helix(l[21:25].strip(), l[33:37].strip(), l[19]))
                except KeyError: helices[l[19]] = [Helix(l[21:25].strip(), l[33:37].strip(), l[19])]
    else:
        #CIF
        ssdict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(coordfile)
        #print(ssdict["_struct_conf"])
        sstype = ssdict["_struct_conf.conf_type_id"]
        sschain = ssdict["_struct_conf.beg_auth_asym_id"]
        ssstart = ssdict["_struct_conf.beg_auth_seq_id"]
        ssend = ssdict["_struct_conf.end_auth_seq_id"]
        all_ss = zip(sstype, sschain, ssstart, ssend)
        for ss in all_ss:
            if "HELX" in ss[0]:
                try: helices[ss[1]].append(Helix(ss[2], ss[3], ss[1]))
                except KeyError: helices[ss[1]] = [Helix(ss[2], ss[3], ss[1])]
    return helices

def get_fasta_mapping(files, index=1):
    if type(files) is str: files = [files]
    mapping = {}
    for fn in files:
        parseme = None
        if "pdb" in fn:
            parseme = Bio.PDB.PDBParser()
            #print("Resolution:",Bio.PDB.parse_pdb_header(fn)["resolution"])
            
        else:
            parseme = Bio.PDB.MMCIFParser()
            #print(can_seq(fn, "cif"))

        stucco = parseme.get_structure(fn[0:fn.index(".")], fn)
        #try: stucco = parseme.get_structure(fn[0:fn.index(".")], fn)
        #except Bio.PDB.PDBExceptions.PDBConstructionException: 
            
        mapping[stucco.get_id()] = {}
        #TODO: Fix broken PDBs

        fastas = ""
        for c in stucco.get_chains():
            #renumbering stuff
            n = index
            mapping[c.get_full_id()[0]][c.get_id()] = {}

            #residue-level stuff
            for i in c.get_residues():
                if Bio.PDB.is_aa(i):
                    #renumbering stuff
                    mapping[c.get_full_id()[0]][c.get_id()][n] = i.get_id()[1]
                    n += 1
    return mapping

def get_fasta(files, index=1):
    if type(files) is str: files = [files]
    fastas = ""
    for fn in files:
        parseme = None
        if "pdb" in fn:
            parseme = Bio.PDB.PDBParser()
            #print("Resolution:",Bio.PDB.parse_pdb_header(fn)["resolution"])
            
        else:
            parseme = Bio.PDB.MMCIFParser()
            #print(can_seq(fn, "cif"))

        stucco = parseme.get_structure(fn[0:fn.index(".")], fn)

        for c in stucco.get_chains():
            l = 0
            chainseq = ""

            #residue-level stuff
            for i in c.get_residues():
                if Bio.PDB.is_aa(i):
                    chainseq += Bio.PDB.protein_letters_3to1[i.get_resname()]
                    l += 1
                    if not l % 80: chainseq += "\n"
            if chainseq:
                chainhead = ">" + c.get_full_id()[0].upper() + ":" + c.get_id() + "\n"
                fastas += chainhead + chainseq.strip() + "\n"
            
    return fastas
def get_chains():
    chains = []

    starting_objs = pymol.cmd.get_names("all", False)
    enabled_objs = pymol.cmd.get_names("all", True)

    for i in enabled_objs:
        pymol.cmd.split_chains(i)

    ending_objs = pymol.cmd.get_names("all", False)


    chains = list(set(starting_objs)^set(ending_objs))
    for i in chains: pymol.cmd.delete(i)

    for i in enabled_objs: pymol.cmd.enable(i)

    actual_chains = {}
    for c in chains:
        try: actual_chains[c[0:-2]].append(c[-1])
        except KeyError: actual_chains[c[0:-2]] = [c[-1]]

    return actual_chains

def ph(helxdict, color="red", fix=False):
    num = 0
    for o in helxdict.keys():
        if fix:
            pymol.cmd.save(o + ".tmp.pdb", o)
            m = get_fasta_mapping(o + ".tmp.pdb")
            os.remove(o + ".tmp.pdb")
            for c in helxdict[o].keys():
                for h in helxdict[o][c]:
                    a = str(m[o][c][int(h.start)])
                    b = str(m[o][c][int(h.end)])
                    pymol.cmd.color(color, o + " and c. " + c + " and i. " + a + "-" + b)
                    num += 1
                num = 0
        else:
            for c in helxdict[o].keys():
                for h in helxdict[o][c]:
                    #pymol.cmd.color("red", o + " and c. " + c + " and i. " + str(h))
                    a = str(h.start)
                    b = str(h.end)
                    pymol.cmd.color(color, o + " and c. " + c + " and i. " + a + "-" + b)
                    num += 1
                num = 0

def tms_test(base="gray", helix="yellow", tms_helix="green"):
    pymol.cmd.color(base, "enabled")
    ph(stride(), helix, fix=False)
    ph(hmmtop(), tms_helix, fix=True)

def paint_test(selection="sele", gray=0):
    if selection not in pymol.cmd.get_names("all", True): selection = "(all)"

    if str2bool(gray): pymol.cmd.color("gray", selection)

    pymol.cmd.color("white", selection + "and ss H and r. LEU+ALA+ILE+VAL+PHE")
    pymol.cmd.color("blue", selection + " and r. ARG+LYS")
    pymol.cmd.color("purple", selection + "and ss H and r. TRP+TYR")

def paint_tmss(start_hue=0, end_hue=240, expand=0, shade=0.8, termini=False, gray=False, offset=0):
    """
DESCRIPTION

    "paint_tmss" calculates probable TMSs using HMMTOP and paints them sequentially

USAGE

    paint_tmss[ start_hue[, end_hue[, expand[, shade[, termini[, gray]]]]]]

ARGUMENTS

    start_hue = int: First hue in gradient {default: 0}

    end_hue = int: Ending hue in gradient {default: 240}

    expand = int: Number of residues to expand predicted TMSs in each direction {default: 0}

    shade = float: How much to shade each additional chain {default: 1.0}

    termini = bool: Whether to paint first/last TMS according to helix orientation

    gray = bool: Whether to gray out the rest of the structure(s)

    offset = int: How much to shift starting hue  for each additional chain {default: 0}

SEE ALSO

    paint_tmss_orig
    """
    stuff = hmmtop()
    start_hue = int(start_hue)
    end_hue = int(end_hue)
    expand = int(expand)
    shade = float(shade)
    offset = int(offset)
    shade_now = 0.9 * 100
    objid = 0

    if str2bool(gray): pymol.cmd.color("gray")

    for o in stuff.keys():
        pymol.cmd.save(o + ".tmp.pdb", o)
        m = get_fasta_mapping(o + ".tmp.pdb")
        os.remove(o + ".tmp.pdb")


        for c in stuff[o].keys():
            for i in range(len(stuff[o][c])):
                hc = (gradient(len(stuff[o][c]), start_hue + objid*offset, end_hue + objid*offset, v=int(shade_now)), stuff[o][c][i])
            #for hc in zip(gradient(len(stuff[o][c]), start_hue, end_hue, v=shade_now), stuff[o][c]):
                a = str(int(m[o][c][int(hc[1].start)]) - expand)
                b = str(int(m[o][c][int(hc[1].end)]) + expand)
                pymol.cmd.color(hc[0][i], o + " and c. " + c + " and i. " + a + "-" + b)
                #print([hc[0][i], o + " and c. " + c + " and i. " + a + "-" + b])
                if termini:
                    pymol.cmd.color("nitrogen", o + " and c. " + c + " and i. " + str(int(a)-1+1))
                    pymol.cmd.color("oxygen", o + " and c. " + c + " and i. " + str(int(b)+1-1))
            shade_now *= shade
            objid += 1

#next colorer:
#accepts (id1, id2, (helixin1, helixin2), (helixin1, helixin2), (helixin1, helixin2))

def paint_tmss_orig(filename, start_hue=0, end_hue=240, expand=0, shade=0.8, termini=False, gray=False, offset=0):
    """
DESCRIPTION

    paint_tmss_orig colors TMSs based on associated UniProt sequences

USAGE

    paint_tmss_orig filename[, start_hue[, end_hue[, expand[, shade[, termini[, gray[, offset]]]]]]]

ARGUMENTS

    filename = str: File to check for associated UniProt accessions

    start_hue = int: First hue in gradient {default:0}

    end_hue = int: Last hue in gradient {default:240}

    expand = int: Number of residues to expand predicted TMSs in each direction {default:0}

    shade = float: How much to shade each successive chain colored {default:0.8}

    termini = bool: Whether to paint first/last TMS residues to display helix orientation

    gray = bool: Whether to gray out the rest of the structure(s)

    offset = int: How much to shift the starting hue for each successive chain colored {default:0}

SEE ALSO

    paint_tmss
    """

    start_hue = int(start_hue)
    end_hue = int(end_hue)
    expand = int(expand)
    shade = float(shade)
    offset = int(offset)
    shade_now = 0.9 * 100
    objid = 0

    if str2bool(gray): pymol.cmd.color("gray")

    seqs = {}

    if "pdb" in filename: 
        f = open(filename)
        ids = {}
        for line in f:
            if "DBREF" not in line: continue
            if line[26:29] != "UNP": continue
            ids[line[12]] = line[33:39]
        f.close()
        for i in sorted(ids):
            s = subprocess.check_output(["curl", "http://www.uniprot.org/uniprot/"+ids[i]+".fasta"])
            s = s.split("\n")
            seqs[i] = ">" + filename + ":" + i + "\n"
            for ss in s[1:]:
                seqs[i] += ss

    else: 
        parseme = Bio.PDB.MMCIF2Dict.MMCIF2Dict(filename)
        raw_seqs = zip(\
parseme["_pdbx_poly_seq_scheme.pdb_strand_id"],\
parseme["_pdbx_poly_seq_scheme.mon_id"],\
parseme["_pdbx_poly_seq_scheme.auth_seq_num"])#,\

        for l in raw_seqs:
            try:
                seqs[l[0]] += Bio.PDB.protein_letters_3to1[l[1]]
            except KeyError:
                seqs[l[0]] = ">" + filename + ":" + l[0] + "\n" + Bio.PDB.protein_letters_3to1[l[1]]
    if len(pymol.cmd.get_names("objects", True)) == 1: target = pymol.cmd.get_names("objects", True)[0]

    helices = {}

    try: o = target
    except NameError: o = filename[:-4]

    for k in sorted(seqs.keys()):

        #prepare to run hmmtop

        p = subprocess.Popen(["hmmtop", "-if=--", "-of=--"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #pipe in the sequence, pipe out the output

        tmss, junk = p.communicate(seqs[k])

        tmss = tmss.split()

        #o = tmss[2][:-2]
        c = tmss[2][-1]
        try: helices[o][c] = []
        except KeyError: helices[o] = {c:[]}

        for i in range(5, len(tmss) - 1, 2):
            helices[o][c].append(Helix(tmss[i], tmss[i+1], c))

    for c in helices[o].keys():
        for i in range(len(helices[o][c])):
            hc = (gradient(len(helices[o][c]), start_hue + objid*offset, end_hue + objid*offset, v=int(shade_now)), helices[o][c][i])
        #for hc in zip(gradient(len(helices[o][c]), start_hue, end_hue, v=shade_now), helices[o][c]):
            #a = str(int(helices[o][c][int(hc[1].start)]) - expand)
            #b = str(int(helices[o][c][int(hc[1].end)]) + expand)
            a = str(int(helices[o][c][i].start) - expand)
            b = str(int(helices[o][c][i].end) + expand)
            pymol.cmd.color(hc[0][i], o + " and c. " + c + " and i. " + a + "-" + b)
            #print([hc[0][i], o + " and c. " + c + " and i. " + a + "-" + b])
            if termini:
                pymol.cmd.color("nitrogen", o + " and c. " + c + " and i. " + str(int(a)-1+1))
                pymol.cmd.color("oxygen", o + " and c. " + c + " and i. " + str(int(b)+1-1))
        shade_now *= shade
        objid += 1

pymol.cmd.extend("paint_tmss", paint_tmss)
pymol.cmd.extend("paint_tmss_orig", paint_tmss_orig)
pymol.cmd.extend("tms_paint", paint_test)
