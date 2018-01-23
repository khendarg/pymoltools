#!/usr/bin/env python2
"""
Copyleft 2016 by Kevin Hendargo, all rights reversed
"""
import pymol



def focus(selection="sele", cutoff=5):
    """
DESCRIPTION

    Hide everything beyond a certain distance away from a target.

USAGE

    focus [selection [, cutoff]]

ARGUMENTS

    selection = string: selection of interest {default: sele}

    cutoff = float: distance to start hiding at {default: 5}
    """
    if selection not in pymol.cmd.get_names("all", True): selection = "(all)"
    pymol.cmd.hide("everything", "(all) beyond " + str(cutoff) + " of " + selection)

def chview(view="cartoon", selection="sele", hidelig=False):
    """
DESCRIPTION

    Switch view for selection in one command.

USAGE

    cv [view [, selection]]

ARGUMENTS

    view = string: view to switch to {default: cartoon}

    selection = string: selection to switch views for {default: sele}
    """
    if selection not in pymol.cmd.get_names("all", True): selection = "(all)"
    if not hidelig: selection += ' and not het'
    print('test')
    pymol.cmd.hide("everything", selection)
    pymol.cmd.show(view, selection)

def slab(thickness=10., selection='all'):
	"""
DESCRIPTION

    Change slab size. Useful for viewing pores
	"""
    try: thickness = float(thickness)
    except ValueError: 
        print("Could not parse %s as a float" % thickness)
        return

    pymol.cmd.clip(mode="slab", distance=thickness, selection=selection)

pymol.cmd.extend("focus", focus)
pymol.cmd.extend("cv", chview)
pymol.cmd.extend("cl", slab)
