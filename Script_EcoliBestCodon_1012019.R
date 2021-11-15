#VipulMadahar 1012019
#Downloaded and Installed Libraries-----------------------------------------------------
library(stringr);
library(seqinr);
#---------------------------------------------------------------------------------------

#Work Directory-------------------------------------------------------------------------
setwd("C:/Users/vipul/OneDrive/Documents/Research/PIAS/Optimized");
#---------------------------------------------------------------------------------------

##Read in DNA Sequences------------------------------------------------
DNASeqRef = toupper(readLines("PIAS.txt"));
##----------------------------------------------------------------------

##Make a DNA Vector-----------------------------------------------------
DNAVectorSt <- s2c(paste(DNASeqRef, collapse = ""));
##----------------------------------------------------------------------

##Make a Vector of Codons-----------------------------------------------
CodonSeqSt = splitseq(DNAVectorSt, frame = 0, word = 3);
#Create Reference Sequence
CodonSeqXX = CodonSeqSt;
##-----------------------------------------------------------------------

#CodonOptimization Loop - All set to best codon for E.Coli --------------
LengthOfLoop <- length(CodonSeqXX);
i = 2;#Start on the second codon as the first should be start codon

for (i in 2:LengthOfLoop){
	#Ala
	if (CodonSeqSt[i] == "GCG" | CodonSeqSt[i] == "GCA" | CodonSeqSt[i] == "GCT" | CodonSeqSt[i] == "GCC"){
		CodonSeqXX[i] = "GCG"
	}
	#Gly
	if (CodonSeqSt[i] == "GGG" | CodonSeqSt[i] == "GGA" | CodonSeqSt[i] == "GGT" | CodonSeqSt[i] == "GGC"){
		CodonSeqXX[i] = "GGC"
	}
	#Arg
	if (CodonSeqSt[i] == "CGT" | CodonSeqSt[i] == "CGC" | CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
		CodonSeqXX[i] = "CGC"
	}
	#His
	if (CodonSeqSt[i] == "CAT" | CodonSeqSt[i] == "CAC"){
		CodonSeqXX[i] = "CAT"
	}
	#Lys
	if (CodonSeqSt[i] == "AAA" | CodonSeqSt[i] == "AAG"){
		CodonSeqXX[i] = "AAA"
	}
	#Glu
	if (CodonSeqSt[i] == "GAG" | CodonSeqSt[i] == "GAA"){
		CodonSeqXX[i] = "GAA"
	}
	#Ser
	if (CodonSeqSt[i] == "AGT" | CodonSeqSt[i] == "AGC" | CodonSeqSt[i] == "TCT" | CodonSeqSt[i] == "TCC" | CodonSeqSt[i] == "TCA" | CodonSeqSt[i] == "TCG"){
		CodonSeqXX[i] = "AGC"
	}
	#Thr
	if (CodonSeqSt[i] == "ACT" | CodonSeqSt[i] == "ACC" | CodonSeqSt[i] == "ACA" | CodonSeqSt[i] == "ACG"){
		CodonSeqXX[i] = "ACC"
	}
	#Asn
	if (CodonSeqSt[i] == "AAT" | CodonSeqSt[i] == "AAC"){
		CodonSeqXX[i] = "AAT"
	}
	#Gln
	if (CodonSeqSt[i] == "CAA" | CodonSeqSt[i] == "CAG"){
		CodonSeqXX[i] = "CAG"
	}
	#Cys
	if (CodonSeqSt[i] == "TGT" | CodonSeqSt[i] == "TGC"){
		CodonSeqXX[i] = "TGC"
	}
	#Pro
	if (CodonSeqSt[i] == "CCT" | CodonSeqSt[i] == "CCC" | CodonSeqSt[i] == "CCA" | CodonSeqSt[i] == "CCG"){
		CodonSeqXX[i] = "CCG"
	}
	#Val
	if (CodonSeqSt[i] == "GTG" | CodonSeqSt[i] == "GTA" | CodonSeqSt[i] == "GTT" | CodonSeqSt[i] == "GTC"){
		CodonSeqXX[i] = "GTG"
	}
	#Ile
	if (CodonSeqSt[i] == "ATT" | CodonSeqSt[i] == "ATC" | CodonSeqSt[i] == "ATA"){
		CodonSeqXX[i] = "ATT"
	}
	#Leu
	if (CodonSeqSt[i] == "TTA" | CodonSeqSt[i] == "TTG" | CodonSeqSt[i] == "CTT" | CodonSeqSt[i] == "CTC" | CodonSeqSt[i] == "CTA" | CodonSeqSt[i] == "CTG"){
		CodonSeqXX[i] = "CTG"
	}
	#Phe
	if (CodonSeqSt[i] == "TTT" | CodonSeqSt[i] == "TTC"){
		CodonSeqXX[i] = "TTT"
	}
	#Tyr
	if (CodonSeqSt[i] == "TAT" | CodonSeqSt[i] == "TAC"){
		CodonSeqXX[i] = "TAT"
	}
	#Trp
	if (CodonSeqSt[i] == "TGG"){
		CodonSeqXX[i] = "TGG"
	}
	#StopCodon
	if (CodonSeqSt[i] == "TGA" | CodonSeqSt[i] == "TAG"){
		CodonSeqXX[i] = "TAA"

	i = i + 1;
}

#Concatenante optimized codon vector & output into work folder------------------------------------------------------ 
OpSeq = paste(CodonSeqXX,collapse="")
write(OpSeq, file = "Ecoli_BestCodon_Op.txt", sep = "")
#-------------------------------------------------------------------------------------------------------------------
#VipulMadahar1012019