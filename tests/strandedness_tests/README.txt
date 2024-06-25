F1R2_pair.bam (edit at chr17 43044352)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Paired ends where the ends overlap. They both contain the same G>A conversion, so this edit should only be counted once,
and the read should only be counted once. So at chr17 43044352 we should see 1 edit and 1 read depth coverage, as a G>A. 

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
