#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import pandas as pd

tmp_blocks = sys.argv[1]
outfile_blocks = sys.argv[2]
outfile_super = sys.argv[3]


df = pd.read_csv(tmp_blocks, engine='python', sep='\t', names=['chr', 'type', 'start', 'end', 'len'])
blocks = df[df.type == 'block']
blocks.index = range(len(blocks))


## break up large blocks (> 100kb) into blocks between 50 - 99 kb


max_len = 100000
inter = max_len/2
new = []
blocks_small = blocks.copy()


for i in blocks_small.index:
    if blocks_small.loc[i, 'len'] > max_len:
        splits = int(blocks_small.loc[i, 'len'] / (inter))
        rest = int(blocks_small.loc[i, 'len'] % inter)

        for j in range(splits):
            
            if j != splits - 1:
                new.append([blocks_small.loc[i, 'chr'], blocks_small.loc[i, 'type'], blocks_small.loc[i, 'start'] + j * inter, blocks_small.loc[i, 'start'] + (j+1) * inter, inter])
            
            else:
                new.append([blocks_small.loc[i, 'chr'], blocks_small.loc[i, 'type'], blocks_small.loc[i, 'start'] + j * inter, blocks_small.loc[i, 'start'] + (j+1) * inter + rest, inter+rest])
        
        blocks_small = blocks_small.drop(i)
        
blocks_small = blocks_small.append(pd.DataFrame(new, columns=['chr', 'type', 'start', 'end', 'len']), ignore_index=True)


# define superblocks overlapping 300 - 900 bp


new = []
counter = 0

min_len = 12000

min_overlap = 300
max_overlap = 3*min_overlap

if (max_overlap/2) > (max_overlap - min_overlap):
    print('ERROR: Logic downstream of here is not valid any more - increase min_overlap size')


for chrom in np.unique(blocks_small.chr.values):
            
    start = 1
    sblock_len = 0
    sblock_nb = 0

    
    tmp = blocks_small[blocks_small.chr == chrom]
    tmp = tmp.sort_values(by='start')
    tmp.index = range(len(tmp))
    
    end_counter = 0
            
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

                
                if ((tmp.loc[i, 'len'] >= min_overlap) & (tmp.loc[i, 'len'] <= max_overlap)):
                    sblock_start = tmp.loc[i, 'start']
                    sblock_len = tmp.loc[i, 'len']
                    sblock_nb = 1
                    
                elif tmp.loc[i, 'len'] > max_overlap:
                    sblock_start = tmp.loc[i, 'end'] - max_overlap/2
                    sblock_len = max_overlap/2
                    sblock_nb = 1
                    
                elif tmp.loc[i, 'len'] < min_overlap:
                                        
                    sblock_len = tmp.loc[i, 'len']
                    sblock_number = 1
                    travel = i-1
                    
                    while sblock_len < min_overlap:
                        
                        if sblock_len + tmp.loc[travel, 'len'] < min_overlap:
                            sblock_len += tmp.loc[travel, 'len']
                            sblock_nb += 1
                            
                            travel -= 1
                        
                        elif sblock_len + tmp.loc[travel, 'len'] <= max_overlap:
                            sblock_start = tmp.loc[travel, 'start'] 
                            sblock_len += tmp.loc[travel, 'len']
                            sblock_nb += 1
                            
                        elif sblock_len + tmp.loc[travel, 'len'] > max_overlap:                            
                            sblock_start = tmp.loc[travel, 'end'] - max_overlap/2
                            sblock_len += max_overlap/2
                            sblock_nb += 1
                                                        
                            
                        else:
                            print('ERROR: Scenario happened at index {} that I didnt account for...'.format(i))
                                            
                
        if ((end_counter == len(tmp)) & (sblock_nb > 1)):
            new.append([chrom, sblock_start, sblock_end, sblock_len, sblock_nb, 'superblock_{}'.format(counter), 'results/blocks/'])
            counter += 1
            
superblocks = pd.DataFrame(new, columns=['chr', 'start', 'end', 'len', 'number_of_blocks', 'name', 'wdir'])

blocks_small.to_csv(outfile_blocks, sep='\t', index=None)
superblocks.to_csv(outfile_super, sep='\t', index=None)
