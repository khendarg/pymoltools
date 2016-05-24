#!/usr/bin/env python2
"""
Copyleft 2016 Kevin Hendargo, all rights reversed.
"""
from __future__ import print_function

import pymol
import Bio.PDB
import re, os
import subprocess


def cat(fn):
    fh = open(fn)
    s = ""
    for l in fh:
        s += l
    fh.close()
    return s

def h2rgb(h, shade=1):
    v = shade
    r = int(h/60.)
    c = ()
    #r - y, 0
    if r == 0: c = (1.0, h/60., 0.0)
    #y - g, 60
    elif r == 1: c = (1 - (h - 60)/60., 1.0, 0)
    #g - c, 120
    elif r == 2: c = (0.0, 1.0, (h - 120)/60.)
    #c - b, 180
    elif r == 3: c = (0.0, 1 - (h - 180)/60., 1.0)
    #b - m, 240
    elif r == 4: c = ((h - 240)/60., 0.0, 1.0)
    #m - r, 300
    elif r == 5: c = (1.0, 0.0, 1 - (h - 300)/60.)

    #print(str.format("0x%x%x%x", ((int(c[0]*255*v), int(c[1]*255*v), int(c[2]*255*v)))))
    return "0x{:02x}{:02x}{:02x}".format(int(c[0]*255*v), int(c[1]*255*v), int(c[2]*255*v))
    return (int(c[0]*255*v), int(c[1]*255*v), int(c[2]*255*v))

def dist2rgb(dist, shade=1):
    if dist > 6: return h2rgb(0, shade)
    elif dist < 1: return h2rgb(240, shade)
    #dist = [1, 6]
    #h = [240, 0]
    #144-24dist
    else: return h2rgb(288 - 48.*dist)

def do_ssm(fixed="", mobile="", remove=0, expand=2, color=0, clean=0):
    """
DESCRIPTION

    "do_ssm" performs SSM superpose and provides visualization of the results

USAGE

    do_ssm fixed, mobile[, remove[, expand[, color]]]

ARGUMENTS

    fixed = string: which object or selection to use as the reference structure

    mobile = string: which object or selection to use as the query structure

    remove = float: what to do with unaligned residues

    expand = int: how many residues to include in addition to the aligned residues

    color = int: whether to color the results by distance

NOTES

    The remove argument controls what is done to unaligned regions as follows:

        0: Hide unaligned regions
        1: Remove unaligned regions
        2: Gray out unaligned regions
        3-4: Fade out unaligned regions. The decimal component controls the transparency of the unaligned residues.

TODO

    Actually use the resulting transformation matrix to move the mobile structure (and print out the RMSD and coverage)
    """
    if not (fixed and mobile) and len(pymol.cmd.get_names("objects", True)) == 2:
        fixed = pymol.cmd.get_names("objects", True)[0]
        mobile = pymol.cmd.get_names("objects", True)[1]
    elif not (fixed and mobile):
        raise ValueError
    #pymol.cmd.save(fixed + ".tmp.pdb", fixed)
    #pymol.cmd.save(mobile + ".tmp.pdb", mobile)
    pymol.cmd.save(fixed + ".pdb", fixed)
    pymol.cmd.save(mobile + ".pdb", mobile)
    #cmd = ["superpose", fixed + ".tmp.pdb", mobile + ".tmp.pdb", "-o", fixed + "_fix_" + mobile + "_mov.pdb"]
    #cmd = ["superpose", mobile + ".tmp.pdb", fixed + ".tmp.pdb", "-o", fixed + "_fix_" + mobile + "_mov.pdb"]
    cmd = ["superpose", mobile + ".pdb", fixed + ".pdb", "-o", fixed + "_fix_" + mobile + "_mov.pdb"]
    out = subprocess.check_output(cmd)
    if not clean:
        fh = open(fixed + "_fix_" + mobile + "_mov.rmsd", "w")
        fh.write(out)
        fh.close()
    pymol.cmd.load(fixed + "_fix_" + mobile + "_mov.pdb")
    #os.remove(fixed + ".tmp.pdb")
    #os.remove(mobile + ".tmp.pdb")
    results = parse_ssm(out, remove, expand, color)
    print(results[0])
    print(results[1])
    return results

def parse_ssm_file(filename, remove=0, expand=2, color=0):
    """
DESCRIPTION

    "parse_ssm_file" parses SSM superpose output for visualization

USAGE

    parse_ssm_file filename[, remove[, expand[, color]]]

ARGUMENTS

    filename = string: file to read SSM output from

    remove = float: what to do with unaligned residues

    expand = int: how many residues to include in addition to the aligned residues

    color = int: whether to color the results by distance

NOTES

    The remove argument controls what is done to unaligned regions as follows:

        0: Hide unaligned regions
        1: Remove unaligned regions
        2: Gray out unaligned regions
        3-4: Fade out unaligned regions. The decimal component controls the transparency of the unaligned residues.
    """
    return parse_ssm(cat(filename), remove, expand, color, loaded=0)

def parse_ssm(output, remove=0, expand=2, color=0, loaded=1):
    #remove levels:
    #0: hide unaligned regions
    #1: remove unaligned regions
    #2: gray out unaligned regions
    #3: fade out unaligned regions
    
    expand = int(expand)

    remove = float(remove)
    transparency = 0.5
    if 3 < remove < 4:
        transparency = remove - 3
    remove= int(remove)

    rec_matrix = 0
    querytgt = 0
    strucs = []
    found_querytgt = False

    rmsd = ""
    cov = ""
    movmatrix = []
    found_rmsd = False
    found_cov = False
    found_matrix = False

    resalig = 0
    q_pos = []
    q_neg = []
    t_pos = []
    t_neg = []
    found_alig = False

    q_cols = {}
    t_cols = {}

    dist = 100
    #qbuffer=[]
    #tbuffer=[]

    for l in output.split("\n"):
        #get fixed and mobile structure names
        if not found_querytgt and "Alignment results" in l: 
            querytgt = 2
            continue
        if querytgt: 
            strucs.append(l.strip().split("/")[-1].split()[-1])
            strucs[-1] = strucs[-1][0:strucs[-1].index(".")]
            querytgt -= 1
            if not querytgt:
                strucs[0] = strucs[1] + "_fix_" + strucs[0] + "_mov"
            found_querytgt = True

        #get matrix
        #if not found_matrix and "Rx" and "Ry" and "Rz" in l: 
        #    rec_matrix = 3
        #    continue
        #
        #if rec_matrix:
        #    movmatrix.append(l.strip().split())
        #    rec_matrix -= 1
        #    found_matrix = True

        #get rmsd
        if not found_rmsd and "r.m.s." in l: 
            #rmsd = float(re.findall("\d+.\d+", l)[0])
            #print(l)
            rmsd = l.strip()
            found_rmsd = True

        #get coverage
        if not found_cov and "Nalign" in l:
            #cov = int(re.findall("\d+", l)[0])
            #print(l)
            cov = l.strip()
            found_cov = True

        #get alignment
        if not found_alig and "Residue alignment" in l:
            found_alig = True

        if found_alig:
            if "Notations" in l: break

            results = l.strip().split("|")
            while '' in results: results.remove('')
            if len(results) == 3:
                q = re.findall("\d+", results[0])
                t = re.findall("\d+", results[2])
                if results[1].strip(): 
                    try:
                        dist = float(re.findall("[0-9]+.[0-9]+", results[1])[0])
                    except IndexError: pass
                    if q:
                        q_pos.append(int(q[0]))
                        if color:
                            try: q_cols[int(q[0])] = dist2rgb(dist)
                            except ValueError: pass
                            try: t_cols[int(t[0])] = dist2rgb(dist)
                            except ValueError: pass
                    if t: 
                        t_pos.append(int(t[0]))
                        if color:
                            try: q_cols[int(q[0])] = dist2rgb(dist, 0.8)
                            except ValueError: pass
                else:
                    try: q_neg.append(re.findall("\d+", results[0])[0])
                    except IndexError: pass
                    try: t_neg.append(re.findall("\d+", results[2])[0])
                    except IndexError: pass
    #expansion
    if expand:
        newq = set(q_pos)
        newt = set(t_pos)
        for i in q_pos: 
            for x in range(i - expand, i + expand + 1):
                newq = newq.union({x})
                try: q_neg.remove(x)
                except ValueError: pass
        for i in t_pos:
            for x in range(i - expand, i + expand + 1): 
                newt = newt.union({x})
                try: t_neg.remove(x)
                except ValueError: pass
        newq = list(newq)
        newt = list(newt)
    else:
        newq = q_pos
        newt = t_pos

    #basic shows/removals
    if not loaded:
        pymol.cmd.load(strucs[0] + ".pdb")
        pymol.cmd.load(strucs[1] + ".pdb")

    if remove == 0: pymol.cmd.hide("everything")

    for i in newq: pymol.cmd.show("cartoon", strucs[0] + " and i. " + str(i))
    for i in newt: pymol.cmd.show("cartoon", strucs[1] + " and i. " + str(i))

    if remove == 1:
        for i in q_neg: pymol.cmd.remove(strucs[0] + " and i. " + str(i))
        for i in t_neg: pymol.cmd.remove(strucs[1] + " and i. " + str(i))

    if remove == 2:
        for i in q_neg: pymol.cmd.color("gray", strucs[0] + " and i. " + str(i))
        for i in t_neg: pymol.cmd.color("gray", strucs[1] + " and i. " + str(i))

    if remove == 3:
        pymol.cmd.create(strucs[0] + "_fade", strucs[0])
        pymol.cmd.create(strucs[1] + "_fade", strucs[1])

        pymol.cmd.hide("everything", strucs[0] + "_fade")
        pymol.cmd.hide("everything", strucs[1] + "_fade")
        pymol.cmd.color("gray", strucs[0] + "_fade")
        pymol.cmd.color("gray", strucs[1] + "_fade")

        pymol.cmd.set("cartoon_transparency", transparency, strucs[0] + "_fade")
        pymol.cmd.set("cartoon_transparency", transparency, strucs[1] + "_fade")

        for i in q_neg: pymol.cmd.show("cartoon", strucs[0] + "_fade and i. " + str(i))
        for i in t_neg: pymol.cmd.show("cartoon", strucs[1] + "_fade and i. " + str(i))

    #distance-based coloring
    if color:
        pymol.cmd.color("gray")
        pymol.cmd.color("0x80c080", strucs[0])
        pymol.cmd.color("0x8080c0", strucs[1])
        #for i in q_pos:
        #    pymol.cmd.color("0xc0c0c0", strucs[0] + " and i. " + str(i))
        #for i in t_pos:
        #    pymol.cmd.color("0x404040", strucs[1] + " and i. " + str(i))
        for i in q_cols.keys(): 
            pymol.cmd.color(q_cols[i], strucs[0] + " and i. " + str(i))
        for i in t_cols.keys(): 
            pymol.cmd.color(t_cols[i], strucs[1] + " and i. " + str(i))

    return (rmsd, cov)

def load_and_ssm(file1, file2):
    """
DESCRIPTION

    "load_and_ssm" loads two files and conducts superpose via SSM on them

USAGE

    load_and_ssm file1, file2

ARGUMENTS

    file1 = string: First file to load, fixed/target structure

    file2 = string: Second file to load, mobile structure

NOTES
    """
    pymol.cmd.load(file1)
    pymol.cmd.load(file2)

    base1 = file1.split("/")[-1].split(".")[0]
    base2 = file2.split("/")[-1].split(".")[0]
    return do_ssm(base1, base2, remove=2, expand=5)

pymol.cmd.extend("do_ssm", do_ssm)
pymol.cmd.extend("parse_ssm_file", parse_ssm_file)
pymol.cmd.extend("load_and_ssm", load_and_ssm)
