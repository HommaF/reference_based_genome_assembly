#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

tmp_blocks = sys.argv[1]
outfile_blocks = sys.argv[2]
outfile_super = sys.argv[3]

df = pd.read_csv(tmp_blocks, engine='python', sep='\t', names=['chr', 'type', 'start', 'end', 'len'])

new = []
counter = 0
min_len = 12000

for chrom in np.unique(df.chr.values):
    
    start = 1
    end_counter = 0

    sblock_len = 0
    sblock_nb = 0

    tmp = df[df.chr == chrom]
    tmp = tmp[tmp.type == "block"]
        

    for i in tmp.index:

        end_counter += 1

        if start == 1:
            sblock_start = tmp.loc[i, 'start']
            sblock_end = tmp.loc[i, 'end']
            sblock_len = tmp.loc[i, 'len']

            sblock_nb += 1
            start = 0

        elif ((sblock_len < min_len) | (sblock_nb == 1)):

            sblock_end = tmp.loc[i, 'end']
            sblock_len += tmp.loc[i, 'len']
            sblock_nb += 1

            if ((sblock_len > min_len) & (sblock_nb > 1)):
                new.append([chrom, sblock_start, sblock_end, sblock_len, sblock_nb, 'superblock_{}'.format(counter), 'results/blocks/'])
                counter += 1
                                      
                sblock_start = tmp.loc[i, 'start']
                sblock_len = tmp.loc[i, 'end'] - tmp.loc[i, 'start']
                sblock_nb = 1

        if ((end_counter == len(tmp)) & (sblock_nb > 1)):
            new.append([chrom, sblock_start, sblock_end, sblock_len, sblock_nb, 'superblock_{}'.format(counter), 'results/blocks/'])
            counter += 1
                          
                                  
superblocks = pd.DataFrame(new, columns=['chrom', 'start', 'end', 'len', 'number_of_blocks', 'name', 'wdir'])


df[df.type == 'block'].to_csv(outfile_blocks, sep='\t', index=None)
superblocks.to_csv(outfile_super, sep='\t', index=None)
