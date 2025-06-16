This codes work together with the structures obtained from Ceona's code.
It takes the events found by her algorithm and runs a manual checking for them. Then, it recalculates the max intensity (it was wrong in her calculation) and saves the strucutres with the new values of only the events that are decided to be aurora. The position of the maximum can also be adjusted. This is done in AA_Clean_events_newI.m
The AB_newalt.m recalculates the new altitude of the newI, which is necessary if the position is adjusted in the previous code.

The codes starting with "P_" plot things.
The codes in the folder "Figures Paper" plot the figures in the paper submitted "A statistical study of the O2Atm(0-0) aurora observed by the Swedish satellite MATS".

