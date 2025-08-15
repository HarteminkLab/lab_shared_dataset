Okay!! This is what Nhat and Trung did to generate these files that have Rossi Strand info


Basically in the Rossi chexmix repo, there are peak calls for most TFs but no strand info associated
But Rossi must have the strand info SOMEWHERE since they have a website that orient the motifs
https://yeastepigenome.org/

so we found that here
https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/

but these are each a replicate of any TF that Rossi did and in each zip file there are motif 1 2 or 3 or more or w/e
so if we go into each of these zip file there would be a motif1_fourcolor.bed file. These are the 60bp windows to make the 
colorful motif logos for i.e Abf1 in that yeastepigenome.org website 



But we basically need peakVal from Chexmix AND also strand info.......

So we decided to merge them. Basically, take Abf1
Chexmix has peak calls and peakval for motif 1, 2 and 3 together in just 1 file (~700 sites)
We merge that with motif1_fourcolor.bed, then merge that with motif2_fourcolor.bed, then merge with motif3_fourcolor.bed
Then we take the unique of all of them and put them all together in a file with peakval and strand info
final file for abf1 has like 500 sites. Which is good enough I guess



so.... there is a python script for it. then a batch script to submit the python script as a job

I ran it for all of the TFs available in the supplementary-data-4-210306.xlsx


BAM DONE! Push to github and call it a day! 




Work done on Aug 14th by Nhat and Trung

If you use this info, you owe us coffee :) 


