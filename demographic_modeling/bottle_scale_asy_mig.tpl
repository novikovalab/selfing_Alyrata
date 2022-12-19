//Number of population samples (demes)
2
//Population effective sizes (number of genes)
TWO_NE_PRESENT_POP1
TWO_NE_PRESENT_POP2
//Sample sizes
8
20
//Growth rates
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 migr21
migr12 0
//Migration matrix 1: No migration
0 0
0 0
//Historical events(time, source, sink, migrants, new size, growth rate, migr. matrix)
3 historical event
T_BOT_POP1_END 0 0 0 RESIZE_POP1_END 0 0
T_BOT_POP1_BEG 0 0 0 RESIZE_POP1_BEG 0 0
T_SPL 0 1 1 RESIZE_ANC 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 7.1e-9 OUTEXP
