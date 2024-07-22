#!/usr/bin/env python3
import os
from os.path import join
from collections import defaultdict


# ~~~ Common helper functions for blocking or controls
def test_for_block(groupdata, contrast, blocks):
    """ only want to run blocking on contrasts where all
    individuals are on both sides of the contrast """
    contrastBlock = [ ]
    for con in contrast:
        group1 = con[0]
        group2 = con[1]
        block1 = [ blocks[sample] for sample in groupdata[group1] ]
        block2 = [ blocks[sample] for sample in groupdata[group2] ]
        if len(block1) == len(block2):
            if len(set(block1).intersection(block2)) == len(block1):
                contrastBlock.append(con)
    return contrastBlock


def ctrl_test(ctrl_dict, input_name, in_dir, mode=None):
    sample = join(in_dir, f"{input_name}.Q5DD.RPGC.bw")
    assert mode in ('chip', 'ctrl'), 'Unrecognized input file mode.'
    # assert os.path.exists(sample), f'{sample} sample does not exist!'
    
    if input_name in ctrl_dict:
        norm = join(in_dir, ctrl_dict[input_name] + ".Q5DD.RPGC.bw")
        # assert os.path.exists(norm), f'{norm} control does not exist!'
    else:
        raise ValueError(f'ChIP sample {input_name} missing from input lookup: \n{str(ctrl_dict)}')
    outs = {'chip': sample, 'ctrl': norm}
    return outs[mode]