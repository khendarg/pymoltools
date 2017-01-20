#!/usr/bin/env python2
"""
Copyleft 2017 by Kevin Hendargo, all rights reversed
"""
from __future__ import print_function

import pymol, os, re, Bio.PDB, subprocess, math, sys, tempfile

def process_aa(threeletter):
    y = ""
    try: y = Bio.PDB.protein_letters_3to1[threeletter]
    except AttributeError: y = Bio.PDB.to_one_letter_code[threeletter]
    except KeyError: y = ""
    return y

def sele2fa(selection="sele"):
    fa = ">titleneeded\n"
    space = {"seq":[], "code":process_aa}
    pymol.cmd.iterate("sele and n. CA", "seq.append(code(resn))", space=space)
    for r in space["seq"]: fa += r
    return fa

def xclip(selection="sele"): 
    os.system("echo '%s' | xclip" % sele2fa(selection))

def pbcopy(selection="sele"): os.system("echo '%s' | pbcopy" % sele2fa(selection))

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

            f = tempfile.NamedTemporaryFile(delete=False)
            pymol.cmd.save(f.name, o + " and c. " + c)

            #run stride on the resulting pdb file and remove it

            struck = subprocess.check_output(["stride", f.name])
            f.delete()

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

def hmmtop(selection=None):
    """
DESCRIPTION

    hmmtop runs hmmtop on loaded sequences

NOTES

    Many if not most loaded sequences will be incomplete relative to gene product sequences. Use with extreme caution and monitor sodium intake.
    """
    #this turns out to be unexpected behavior
    #if selection not in pymol.cmd.get_names("all", True): selection = "(all)"
    if selection == None: selection = "(all)"

    enabled = pymol.cmd.get_names("all", True)

    pymol.cmd.disable("all")
    pymol.cmd.enable(selection)

    helices = {}

    for o in pymol.cmd.get_names("objects", True):
        
        for c in pymol.cmd.get_chains(selection + " and " + o):
            fasta = ">%s_%s\n" % (o, c)
            space = {"code":process_aa, "seq":[]}
            pymol.cmd.iterate("%s and %s and c. %s and n. CA" % (selection, o, c), "seq.append(resn)", space=space)
            for r in space["seq"]: fasta += process_aa(r)
            
            p = subprocess.Popen(["hmmtop", "-if=--", "-of=--"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            tmss, stderr = p.communicate(fasta)

            tmss = tmss.split()

            try: helices[o][c] = []
            except KeyError: helices[o] = {c:[]}
            for i in range(5, len(tmss)-1, 2):
                helices[o][c].append(Helix(tmss[i], tmss[i+1], c))

    pymol.cmd.disable("all")
    for x in enabled: pymol.cmd.enable(x)

    return helices
    
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

def get_renumber(selection):
    if selection == None: selection = "(all)"

    pymol.cmd.disable("all")
    pymol.cmd.enable(selection)

    for o in pymol.cmd.get_names("objects", True): pass

    return 1

def ph(helxdict, color="red", fix=False):
    num = 0
    for o in helxdict.keys():
        if fix:
            f = tempfile.NamedTemporaryFile(delete=False)
            pymol.cmd.save(f.name, o)
            m = get_fasta_mapping(f.name)
            f.delete()
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

def paint_tmss(selection=None, gray=1, start_hue=0, end_hue=240, expand=0, shade=0.9, saturation=0.7, termini=False, offset=0):
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

    if selection == None: selection = "all"

    enabled = pymol.cmd.get_names("all", True)

    pymol.cmd.disable("all")
    pymol.cmd.enable(selection)

    for o in pymol.cmd.get_names("objects", False): pymol.cmd.enable(selection + " and " + o)

    #run HMMTOP on the sequences

    stuff = hmmtop(selection)

    #make sure arguments are properly typed

    start_hue = int(start_hue)
    end_hue = int(end_hue)
    expand = int(expand)
    shade = float(shade)
    saturation = float(saturation)
    offset = int(offset)
    shade_now = 0.9 * 100
    objid = 0

    for o in pymol.cmd.get_names("objects", True): 
        #gray out everything if requested
        if bool(gray): pymol.cmd.color(hsv2rgb(0, 0, 60), o + " and " + selection)

        #save each object as a temporary PDB and extract sequence from there
        #is what an incompetent coder would do, so let's do something better

        for c in pymol.cmd.get_chains(o, True): 

            v = {"ri":[], "rn":[]}
            pymol.cmd.iterate(o + " and c. " + c + " and " + selection + " and n. CA", "ri.append(resv); rn.append(resn)", space=v)
            fa = ">%s_%s\n" % (o, c)
            for n in v["rn"]: 
                try: fa += Bio.PDB.protein_letters_3to1[n]
                except AttributeError: fa += Bio.PDB.to_one_letter_code[n]
                except KeyError: continue
            
            resi2flat = {}
            flat2resi = {}
            i = 1
            for j in v["ri"]:
                resi2flat[j] = i
                flat2resi[i] = j
                i += 1
            hc = gradient(len(stuff[o][c]), start_hue, end_hue, int(100*saturation), int(100*shade))

            for x in zip(hc, stuff[o][c]):
                pymol.cmd.color(x[0], "%s and %s and c. %s and i. %s-%s" % (selection, o, c, flat2resi[int(x[1].start)], flat2resi[int(x[1].end)]))
                

#next colorer:
#accepts (id1, id2, (helixin1, helixin2), (helixin1, helixin2), (helixin1, helixin2))

def paint_tmss_orig(filename, selection=None, start_hue=0, end_hue=240, expand=0, shade=0.99, termini=False, gray=False, offset=0):
    """
DESCRIPTION

    paint_tmss_orig colors TMSs based on associated UniProt sequences

USAGE

    paint_tmss_orig filename[, selection[, start_hue[, end_hue[, expand[, shade[, termini[, gray[, offset]]]]]]]

ARGUMENTS

    filename = str: File to check for associated UniProt accessions

    selection = str: Object or selection to color (required if ambiguous)

    start_hue = int: First hue in gradient {default:0}

    end_hue = int: Last hue in gradient {default:240}

    expand = int: Number of residues to expand predicted TMSs in each direction {default:0}

    shade = float: How much to shade each successive chain colored {default:0.9}

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

    #attempt to detect targets slightly intelligently
    if not selection:
        #if more than one selection/object exists and selection is unspecified, exit
        if len(pymol.cmd.get_names("all", True)) == 1: selection = pymol.cmd.get_names("all", True)[0]
        else:
            print("Unspecified selection, exiting")
            sys.exit(1)

    #now gray out anything if gray is set
    if str2bool(gray): pymol.cmd.color("gray", selection)

    seqs = {}

    #check for file type
    #TODO: read last line instead of filename

    if "pdb" in filename: 

        #if legacy PDB, download UNP (UniProt) FASTAs to memory

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
            seqs[i] = ">" + filename + ":" + ids[i] + ":" + i + "\n"
            for ss in s[1:]:
                seqs[i] += ss

        p = Bio.PDB.PDBParser(filename)

        struc = p.get_structure(selection, filename)

        atomseqs = {}

        for c in struc.get_chains():
            atomseqs[c.get_id()] = ">" + filename + ":" + "ATOM:" + c.get_id() + "\n"
            for r in c.get_residues():
                try: atomseqs[c.get_id()] += Bio.PDB.protein_letters_3to1[r.get_resname()]
                except KeyError: pass

        for i in sorted(ids):
            f = open(".seq1.tmp", "w")
            f.write(seqs[i])
            f.close()
            f = open(".seq2.tmp", "w")
            f.write(atomseqs[i])

            f = open(".seq1.tmp")
            f.close()
            f = open(".seq2.tmp")
            f.close()

            os.system("needle -gapopen 10 -gapextend 0.5 -asequence .seq1.tmp -bsequence .seq2.tmp -outfile /dev/stdout")

            os.remove(".seq1.tmp")
            os.remove(".seq2.tmp")

    else:

        #if CIF, turn the CIF file into a dict
        #TODO: look for the UNP DBREF equivalents and download too because MMCIF2Dict is horrendously slow

        parseme = Bio.PDB.MMCIF2Dict.MMCIF2Dict(filename)
        #if network disabled...
        raw_seqs = zip(\
parseme["_pdbx_poly_seq_scheme.pdb_strand_id"],\
parseme["_pdbx_poly_seq_scheme.mon_id"],\
parseme["_pdbx_poly_seq_scheme.auth_seq_num"])#,\

        for l in raw_seqs:

            #initiate a basic FASTA header and elongate

            try:
                seqs[l[0]] += Bio.PDB.protein_letters_3to1[l[1]]
            except KeyError:
                seqs[l[0]] = ">" + filename + ":" + l[0] + "\n" + Bio.PDB.protein_letters_3to1[l[1]]

        #relevant = (\
#parseme["_struct_ref.id"],\
#parseme["_struct_ref.db_name"],\
#parseme["_struct_ref.db_code"],\
#parseme["_struct_ref.pdbx_db_accession"],\
#parseme["_struct_ref.entity_id"])

        #for i in relevant: print(i)

        #return
    helices = {}

    o = selection

    for k in sorted(seqs.keys()):

        #prepare to run hmmtop

        p = subprocess.Popen(["hmmtop", "-if=--", "-of=--"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

        #run HMMTOP, piping in the sequence and piping out the output

        tmss = p.communicate(seqs[k])[0]

        tmss = tmss.split()

        c = tmss[2][-1]

        #begin constructing the dicts containing the helix objects

        try: helices[o][c] = []
        except KeyError: helices[o] = {c:[]}

        #parse HMMTOP output

        for i in range(5, len(tmss) - 1, 2):
            helices[o][c].append(Helix(tmss[i], tmss[i+1], c))

    #recurse through the helices, coloring them and stuff

    for c in helices[o].keys():
        for i in range(len(helices[o][c])):
            hc = (gradient(len(helices[o][c]), start_hue + objid*offset, end_hue + objid*offset, v=int(shade_now)), helices[o][c][i])
            a = str(int(helices[o][c][i].start) - expand)
            b = str(int(helices[o][c][i].end) + expand)
            pymol.cmd.color(hc[0][i], o + " and c. " + c + " and i. " + a + "-" + b)
            if termini:
                pymol.cmd.color("nitrogen", o + " and c. " + c + " and i. " + str(int(a)-1+1))
                pymol.cmd.color("oxygen", o + " and c. " + c + " and i. " + str(int(b)+1-1))
        shade_now *= shade
        objid += 1

pymol.cmd.extend("paint_tmss", paint_tmss)
pymol.cmd.extend("pt", paint_tmss)
pymol.cmd.extend("paint_tmss_orig", paint_tmss_orig)
pymol.cmd.extend("hmmtop", hmmtop)

pymol.cmd.extend("xclip", xclip)
pymol.cmd.extend("pbcopy", pbcopy)

