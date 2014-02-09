# Niels Hanson
# MetaPathways Tutorial: Pathway Analysis

# 1. Move processed ePGDBs to the ptools-local/pgdbs/user/ directory
# - in each, make sure /1.0/input/organism.dat file is updated
# to avoid the unclassified sequences problem

# 2. In another terminal start Pathway Tools in -api mode
~/pathway-tools/pathway-tools -api

# 3. List Available ePGDBs in Pathway Tools
perl extract_pathway_table_from_pgdb.pl -l

# 4. Extract pathways from pathway tools 
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic -out 1_upper_euphotic.lookup.txt -t lookup

# - extract but excude coverage within a certain coverage -c  and support -s 7
# e.g half the reaction covered -c 0.5 and at least 7 ORFs in pathway overall -s 7
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic -out 1_upper_euphotic.lookup.c05.s7.txt -t lookup -c 0.5 -s 7

# - long table format displays each each ORF
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic -out 1_upper_euphotic.long.txt -t long 

# - wide table format of pathways from multiple samplpes
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic 6_upper_euphotic 2_chlorophyllmax 3_below_euphotic 5_uppermesopelagic 7_omz 4_deepabyss -out HOT_Sanger_pwy.wide.txt -t wide 

# - wide table format of rxns 1_upper_euphotic
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic -out 1_upper_euphotic_rxn.wide.txt -t wide -rxn

# - wide table format of rxns from multiple samplpes
perl extract_pathway_table_from_pgdb.pl -f 1_upper_euphotic 6_upper_euphotic 2_chlorophyllmax 3_below_euphotic 5_uppermesopelagic 7_omz 4_deepabyss -out HOT_Sanger_rxn.wide.txt -t wide -rxn
