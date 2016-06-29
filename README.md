# pymoltools

To use, simply unpack tools.tar.gz to your preferred Pymol plugins directory.

Contents:

colorhelices.py

    paint_tmss : Requires HMMTOP, paints helices based on HMMTOP-predicted TMSs

goodies.py

    cv : Shortcut for changing representations

    focus : Shortcut for hiding everything beyond a certain distance away from a
    target

superposter.py

    do_ssm : Attempt SSM superpose on two different objects, returning a new
    transformed object

    load_and_ssm : Load two coordinates files and attempt to superpose them,
    also by SSM

    parse_ssm_file : Load superpose output and display the results

seqtools.py

    seqload : Load extended sequence data from coordinates files (if possible)
