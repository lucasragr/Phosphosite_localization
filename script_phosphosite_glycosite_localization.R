library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)

#FASTA File import
database_fasta = readAAStringSet("human.fasta")

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                #Transform FASTA into a dataframe
seq_name = names(database_fasta)
sequence = paste(database_fasta)
df <- data.frame(seq_name, sequence)

colnames(df)

#Cleanning the table for (ID, sequence)


df1 = separate(data = df, col = seq_name, into = c("trash1", "ID"), sep = "\\|")
tbl = subset(df1, select = -trash1)


#Upload my table

#glycoproteome table

mytbl = read.csv("~/phosphosite_localization/Phosphosite_localization/glyco_tbl.csv", sep=",")

#phosphoproteome table

#mytbl = read.csv("~/phosphosite_localization/Phosphosite_localization/phospho_tbl.csv", sep=",")

#filter glycopeptides 

mytbl = mytbl %>% filter(grepl("d", sequencep))

#filter phosphopeptides 

#mytbl = mytbl %>% filter(grepl("p", sequencep))

#Merge df and mytbl



tbl_merged = right_join(tbl, mytbl, by = "ID")

##locate returns a matrix and we are just interested in the start of keyword.

tbl_merged$position_stringr = str_locate(tbl_merged$sequence, tbl_merged$sequencen)


#calculates the position of phosphorylation into phosphopeptide

#for phosphoproteome
tbl_merged$p_loc = (gregexpr(pattern = "p", tbl_merged$sequencep))

#for glycoproteome
#tbl_merged$p_loc = (gregexpr(pattern = "d", tbl_merged$sequencep))
#tbl_merged$p_S_loc = gregexpr(pattern = "pS", tbl_merged$Sequencep)
#tbl_merged$p_T_loc = gregexpr(pattern = "pT", tbl_merged$Sequencep)
#tbl_merged$p_Y_loc = gregexpr(pattern = "pY", tbl_merged$Sequencep)

#split phosphosites that were phosphorylated or glycosilated in each peptide 

tbl_merged1= separate(data = tbl_merged, col = p_loc, into = c("phosphosite1", "phosphosite2", "phosphosite3", "phosphosite4", "phosphosite5", "phosphosite6", "phosphosite7"), sep = "\\,")


#transform phosphosites in numeric

tbl_merged1 = tbl_merged1  %>%
  mutate_at("phosphosite1", str_replace, "c", "")

tbl_merged1$phosphosite1 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite1)
tbl_merged1$phosphosite2 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite2)
tbl_merged1$phosphosite3 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite3)
tbl_merged1$phosphosite4 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite4)
tbl_merged1$phosphosite5 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite5)
tbl_merged1$phosphosite6 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite6)
tbl_merged1$phosphosite7 = gsub("[[:punct:]]", "", tbl_merged1$phosphosite7)


#Transform the phosphorylation sites to a compatible type of class 
tbl_mcl = NULL
tbl_mcl$pr = as.numeric(tbl_merged1$phosphosite1) 
tbl_mcl$sc = as.numeric(tbl_merged1$phosphosite2) 
tbl_mcl$tr = as.numeric(tbl_merged1$phosphosite3) 
tbl_mcl$fr = as.numeric(tbl_merged1$phosphosite4) 
tbl_mcl$fv = as.numeric(tbl_merged1$phosphosite5) 
tbl_mcl$sx = as.numeric(tbl_merged1$phosphosite6)
tbl_mcl$sv = as.numeric(tbl_merged1$phosphosite7)
tbl_mcl$sequencen = as.character(tbl_merged1$sequencen)
tbl_mcl$position_init = as.numeric(tbl_merged1$position_stringr[,"start"])
tbl_mcl$sequencep = as.character(tbl_merged1$sequencep)
tbl_mcl$ID = as.character(tbl_merged1$ID)
tbl_mcl$sequence = as.character(tbl_merged1$sequence)

#transform to datafram 

tbl_mcl$phosphosite1 = tbl_mcl$pr + tbl_mcl$position_init
tbl_mcl$phosphosite2 = tbl_mcl$sc + tbl_mcl$position_init
tbl_mcl$phosphosite3 = tbl_mcl$tr + tbl_mcl$position_init
tbl_mcl$phosphosite4 = tbl_mcl$fr + tbl_mcl$position_init
tbl_mcl$phosphosite5 = tbl_mcl$fv + tbl_mcl$position_init
tbl_mcl$phosphosite6 = tbl_mcl$sx + tbl_mcl$position_init
tbl_mcl$phosphosite7 = tbl_mcl$sv + tbl_mcl$position_init

#acha a posição de fosforilaçao na proteina

tbl_mcl$phosphosite1 = tbl_mcl$phosphosite1 -1
tbl_mcl$phosphosite2 = tbl_mcl$phosphosite2 -2
tbl_mcl$phosphosite3 = tbl_mcl$phosphosite3 -3
tbl_mcl$phosphosite4 = tbl_mcl$phosphosite4 -4
tbl_mcl$phosphosite5 = tbl_mcl$phosphosite5 -5
tbl_mcl$phosphosite5 = tbl_mcl$phosphosite6 -6
tbl_mcl$phosphosite5 = tbl_mcl$phosphosite7 -7

#fala em qual aminoacido
tbl_mcl$aa_pos1 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite1,tbl_mcl$phosphosite1)
tbl_mcl$aa_pos2 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite2,tbl_mcl$phosphosite2)
tbl_mcl$aa_pos3 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite3,tbl_mcl$phosphosite3)
tbl_mcl$aa_pos4 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite4,tbl_mcl$phosphosite4)
tbl_mcl$aa_pos5 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite5,tbl_mcl$phosphosite5)
tbl_mcl$aa_pos6 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite6,tbl_mcl$phosphosite6)
tbl_mcl$aa_pos7 = substr(tbl_mcl$sequence, tbl_mcl$phosphosite7,tbl_mcl$phosphosite7)

#concatenate aminoacid and position

tbl_mcl$phosphosite1 <- paste(tbl_mcl$aa_pos1, tbl_mcl$phosphosite1)
tbl_mcl$phosphosite2 <- paste(tbl_mcl$aa_pos2, tbl_mcl$phosphosite2)
tbl_mcl$phosphosite3 <- paste(tbl_mcl$aa_pos3, tbl_mcl$phosphosite3)
tbl_mcl$phosphosite4 <- paste(tbl_mcl$aa_pos4, tbl_mcl$phosphosite4)
tbl_mcl$phosphosite5 <- paste(tbl_mcl$aa_pos5, tbl_mcl$phosphosite5)
tbl_mcl$phosphosite6 <- paste(tbl_mcl$aa_pos6, tbl_mcl$phosphosite6)
tbl_mcl$phosphosite7 <- paste(tbl_mcl$aa_pos7, tbl_mcl$phosphosite7)

tbl_mcl = as.data.frame(tbl_mcl)

tbl_mcl1 = as.data.frame(tbl_mcl)


#rename phosphosites columns to site

tbl_mcl1 = tbl_mcl1 %>% 
  rename(site1 = phosphosite1,
         site2 = phosphosite2,
         site3 = phosphosite3,
         site4 = phosphosite4,
         site5 = phosphosite5,
         site6 = phosphosite6,
         site7 = phosphosite7)

###Cleaning table 

tbl_final = select(tbl_mcl1, ID, sequencep, sequencen, site1,
                   site2, site3, site4, site5, site6, site7)




write.csv(tbl_final,file = "glycosite_tbl.csv")




    