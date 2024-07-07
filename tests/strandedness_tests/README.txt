F1R2_pair.bam (edit at chr17 43044352)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Paired ends where the ends overlap. They both contain the same G>A conversion, so this edit should only be counted once,
and the read should only be counted once if MARINE is run in --paired_end mode.
So at chr17 43044352 we should see 1 edit and 1 read depth coverage, as a G>A. If --paired_end mode is not used, the
number of edits will not be double-counted but the coverage will be, so coverage at this position will be 2. 

F2R1_end.bam (edits at chr 17 43001715, 43001716, 43001717)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Paired ends, but F2R1 so that the second read in the pair is actually mapped to the forward strand. This second read
should then have the edits A>G, G>A, and A>G all right next to each other, with coverage 1 only.

same_pos_dif_reads.bam (edits at chr17 83199872)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Same edited position, but across multiple reads with different IDs. There should only be 9
unique C>G edits with 9 total coverage, because many of these are pairs.


tax1bp3_chr17_3665556.bam and tax1bp3_chr17_3665556_read.bam (edits at chr17 3665556)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Should be one G>A edit on one read at this position. But since it is a negative stranded gene, TAX1BP3, ultimate
conversion would be C>T.


pair_test (pair_example_18_49488551_49590000.sorted.bam) (edits in RPL17 (-) at 49491556, edits in LIPG (+) at 49567494)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In RPL17, we have 5 reads but actually have just two different pairs. One pair has both ends overlapping at this 
location, and is duplicated, so there are 4 reads with the edit that are all actually just 1 read. The other pair has
just one end at this location, so it contributes 1 read. So there are 2 edited uniqe pairs and a coverage of 2 for 
RPL17 at position 18:49491556.

In LIPG we have 3 reads but just two different pairs. One pair has both ends overlapping at this location, so there are 
2 reads with the edit that are both actually just different ends of 1 read. The other pair has just one end at this 
location, so it contributes 1 read. So there are 2 edited uniqe pairs and a coverage of 2 for LIPG at position 18:49567494.

