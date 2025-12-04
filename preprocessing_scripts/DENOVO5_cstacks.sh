# In this script we combine the stacks of all the samples into one catalog for the whole species
exec > >(tee -i DENOVO5_cstacks.log) 2>&1

cstacks -P ./DENOVO45_stacks -M ./DENOVO5_PL_Map.txt -n 4 -p 20
