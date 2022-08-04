############# READ BEFORE RUN SCRIPT ###########################

- "Biostrings", "tidyr", "dplyr", "stringr" packages must be installed.

Biostrings package can report a error in installing process, if it is the case run the following code: 


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

for install the other packages run this code: 

install.packages("tidyr")
install.packages("dplyr")
install.packages("stringr")

############# RUNNING THE SCRIPT ###########################

For execute this code, your table must be in csv format "your_tbl.csv".

your_tbl must have the following columns "ID", "sequencep", "sequencen". 

sequencep = columns containning the aminoacid sequence with "p" before the modified aminoacid. 
sequencen = aminoacid sequence. 

you have to rename your FASTA file to "human.FASTA", or rename the code in script. 

IMPORTANT: you must have "your_tbl", "human/mouse/rat.fasta" in the same folder in R project. 
