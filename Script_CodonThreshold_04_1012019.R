#Vipul Madahar 2018/2019
#Downloaded and Installed Libraries-----------------------------------------------------
library(stringr);
library(seqinr);
#---------------------------------------------------------------------------------------

#Work Directory-------------------------------------------------------------------------
setwd("C:/Users/vipul/OneDrive/Documents/Research/CodonOptimization/RNF4");
#---------------------------------------------------------------------------------------

##Read in DNA Sequences-----------------------------------------------------------------
#Creating the variable DNASeqRef to read in the text file with your cDNA
#The variable will be the cDNA sequence of your text file, all in upper case

DNASeqRef = toupper(readLines("RNF4Wt.txt"));
##--------------------------------------------------------------------------------------

##Make a DNA Vector---------------------------------------------------------------------
#Takes in the DNASeqRef and creates a new variable, DNAVectorSt
#which is a 1 by (n) matrix, where "n" is the number of bases

DNAVectorSt <- s2c(paste(DNASeqRef, collapse = ""));
##--------------------------------------------------------------------------------------

##Make a Vector of Codons---------------------------------------------------------------
#Takes in teh DNAVectorSt and now creates it into codons based on +1, as we assume it is 
#cDNA and every 3 base is the codon we will translate. 

CodonSeqSt = splitseq(DNAVectorSt, frame = 0, word = 3);
CodonSeqXX = CodonSeqSt;
##--------------------------------------------------------------------------------------

##Codon Selection based on codon BIAS with CAI threshold of 0.4
#Loop Conditions------------------------------------------------------------------------
LengthOfLoop = length(CodonSeqXX);
i = 2;#Start on the second codon as the first should be start codon

for (i in 2:LengthOfLoop){
	#Leu Check
	if (CodonSeqSt[i] == "TTA" | CodonSeqSt[i] == "TTG" | CodonSeqSt[i] == "CTT" | CodonSeqSt[i] == "CTC" | CodonSeqSt[i] == "CTA"){
		CodonSeqXX[i] = "CTG"
	}	
	#Ile Check
	if (CodonSeqSt[i] == "ATA"){
		CodonSeqXX[i] = "ATC"
		b = i - 1
		if (CodonSeqSt[b] == "ATT" | CodonSeqSt[b] == "ATA" | CodonSeqSt[b] == "ATC"){
			CodonSeqXX[i] = "ATT"
			CodonSeqXX[b] = "ATT"	
		}
	}
	#His Check
	if (CodonSeqSt[i] == "CAC" | CodonSeqSt[i] == "CAT"){
		b = i - 1;
		if(CodonSeqSt[b] == "CAT" | CodonSeqSt[b] == "CAC"){
			CodonSeqXX[i] = "CAC"
			CodonSeqXX[b] = "CAC" 
		}
	}
	#Gly Check
	if (CodonSeqSt[i] == "GGG" | CodonSeqSt[i] == "GGA"){
		CodonSeqSt[i] = "GGT";
		b = i - 1;
		if(CodonSeqSt[b] == "GGG" | CodonSeqSt[b] == "GGA" | CodonSeqSt[b] == "GGT" | CodonSeqSt[b] == "GGC"){
			CodonSeqXX[i] = "GGC"
			CodonSeqXX[b] = "GGC" 
		}
	}
	#Asp Check
	if (CodonSeqSt[i] == "GAT" | CodonSeqSt[i] == "GAC"){
		b = i - 1;
		if(CodonSeqSt[b] == "GAT" | CodonSeqSt[b] == "GAC"){
			CodonSeqXX[i] = "GAC"
			CodonSeqXX[b] = "GAC" 
		}
	}
	#Lys 
	if (CodonSeqSt[i] == "AAG"){
		CodonSeqXX[i] = "AAA"
		b = i - 1;
		if(CodonSeqSt[b] == "AAA" | CodonSeqSt[b] == "AAG"){
			CodonSeqXX[i] = "AAG"
			CodonSeqXX[b] = "AAG" 
		}
	}
	#Glu 
	if (CodonSeqSt[i] == "GAG" | CodonSeqSt[i] == "GAA"){
		b = i - 1;
		if(CodonSeqSt[b] == "GAG" | CodonSeqSt[b] == "GAG"){
			CodonSeqXX[i] = "GAG"
			CodonSeqXX[b] = "GAG" 
		}
	}
	#Ser TC
	if (CodonSeqSt[i] == "TCA" | CodonSeqSt[i] == "TCG"){
		CodonSeqXX[i] = "TCT"
		b = i - 1;
		if(CodonSeqSt[b] == "AGT" | CodonSeqSt[b] == "AGC" | CodonSeqSt[b] == "TCT" | CodonSeqSt[b] == "TCC" | CodonSeqSt[b] == "TCA" | CodonSeqSt[b] == "TCG"){
			CodonSeqXX[i] = "AGC"
			CodonSeqXX[b] = "AGC" 
		}
	}
	#Ser AG
	if (CodonSeqSt[i] == "AGT"){
		CodonSeqXX[i] = "AGC"
		b = i - 1;
		if(CodonSeqSt[b] == "AGT" | CodonSeqSt[b] == "AGC" | CodonSeqSt[b] == "TCT" | CodonSeqSt[b] == "TCC" | CodonSeqSt[b] == "TCA" | CodonSeqSt[b] == "TCG"){
			CodonSeqXX[i] = "AGC"
			CodonSeqXX[b] = "AGC" 
		}
	}
	#Arg Check
	if (CodonSeqSt[i] == "CGT" | CodonSeqSt[i] == "CGC" | CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
		if (CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
			CodonSeqXX[i] = "CGT"
		}
		b = i - 1
		if(CodonSeqSt[b] == "CGT" | CodonSeqSt[b] == "CGC" | CodonSeqSt[b] == "CGA" | CodonSeqSt[b] == "CGG"|CodonSeqSt[b] == "AGA" | CodonSeqSt[b] == "AGG"){
			CodonSeqXX[i] = "CGC"
			CodonSeqXX[b] = "CGC" 
		}
	}
	#Asn Check
	if (CodonSeqSt[i] == "AAT" | CodonSeqSt[i] == "AAC"){
		b = i - 1
		if (CodonSeqSt[b] == "AAT" | CodonSeqSt[b] == "AAC"){
			CodonSeqXX[b] = "AAT"
			CodonSeqXX[i] = "AAT"
		}	
	}
	#Thr Check
	if (CodonSeqSt[i] == "ACT" | CodonSeqSt[i] == "ACC" | CodonSeqSt[i] == "ACA" | CodonSeqSt[i] == "ACG"){
		if (CodonSeqSt[i] == "ACA" | CodonSeqSt[i] == "ACT"){
			CodonSeqXX[i] = "ACC"
		}
		#No repeat check for Thr
	}
	#Gln Check
	if (CodonSeqSt[i] == "CAA" | CodonSeqSt[i] == "CAG"){
		b = i - 1
		if (CodonSeqSt[b] == "CAA" | CodonSeqSt[b] == "CAG"){
			CodonSeqXX[b] = "CAA"
			CodonSeqXX[i] = "CAA"
		}
	}
	#Pro Check
	if (CodonSeqSt[i] == "CCT" | CodonSeqSt[i] == "CCC" | CodonSeqSt[i] == "CCA" | CodonSeqSt[i] == "CCG"){
		CodonSeqXX[i] = "CCG"
		#No repeat check, Proline BIAS too high
	}
	#Cys Check
	if (CodonSeqSt[i] == "TGT" | CodonSeqSt[i] == "TGC"){
		b = i - 1
		if (CodonSeqSt[b] == "TGT" | CodonSeqSt[b] == "TGC"){
			CodonSeqXX[i] = "TGT"
			CodonSeqXX[b] = "TGT"
		}
	}
	#Ala Check
	if (CodonSeqSt[i] == "GCG" | CodonSeqSt[i] == "GCA" | CodonSeqSt[i] == "GCT" | CodonSeqSt[i] == "GCC"){
		if (CodonSeqSt[i] == "GCT"){
			CodonSeqXX[i] = "GCG"
		}
		b = i - 1
		if (CodonSeqSt[b] == "GCG" | CodonSeqSt[b] == "GCA" | CodonSeqSt[b] == "GCT" | CodonSeqSt[b] == "GCC"){
			CodonSeqXX[b] = "GCG"
			CodonSeqXX[i] = "GCG"
		}	
	}
	#Val Check
	if (CodonSeqSt[i] == "GTG" | CodonSeqSt[i] == "GTA" | CodonSeqSt[i] == "GTT" | CodonSeqSt[i] == "GTC"){
		if (CodonSeqSt[i] == "GTC" | CodonSeqSt[i] == "GTA"){
			CodonSeqXX[i] = "GTG"
		}
		b = i - 1
		if (CodonSeqSt[b] == "GTG" | CodonSeqSt[b] == "GTA" | CodonSeqSt[b] == "GTT" | CodonSeqSt[b] == "GTC"){
			CodonSeqXX[b] = "GTG"
			CodonSeqXX[i] = "GTG"
		}
	}
	#Phe Check
	if (CodonSeqSt[i] == "TTT" | CodonSeqSt[i] == "TTC"){
		b = i - 1
		if (CodonSeqSt[b] == "TTT" | CodonSeqSt[b] == "TTC"){
			CodonSeqXX[b] = "TTT"
			CodonSeqXX[i] = "TTT"
		}
	}
	#Tyr Check
	if (CodonSeqSt[i] == "TAT" | CodonSeqSt[i] == "TAC"){
		CodonSeqXX[i] = "TAT"
		b = i - 1
		if (CodonSeqSt[b] == "TAT" | CodonSeqSt[b] == "TAC"){	
			CodonSeqXX[i] = "TAT"
			CodonSeqXX[b] = "TAT"
		}
	}
	#Trp Check
	if (CodonSeqSt[i] == "TGG"){
		CodonSeqXX[i] = "TGG"
	}
	#StopCodon
	if (CodonSeqSt[i] == "TGA" | CodonSeqSt[i] == "TAG"){
		CodonSeqXX[i] = "TAA"
	}
	i = i + 1;
}
#------------------------------------------------------------------------


#CodonFrequencyOptimized-------------------------------------------------
OpSeq = paste(CodonSeqXX,collapse="")
write(OpSeq, file = "OPT_CAI04_threshold.txt", sep = "")
#------------------------------------------------------------------------


#Loop Conditions Higher AT% at the start---------------------------------
CodonRNA5P = CodonSeqXX
LengthOfLoop = 16;
i = 2;#Start on the second codon as the first should be start codon

for (i in 2:LengthOfLoop){
	#Ala
	if (CodonRNA5P[i] == "GCG" | CodonRNA5P[i] == "GCA" | CodonRNA5P[i] == "GCT" | CodonRNA5P[i] == "GCC"){
		CodonRNA5P[i] = "GCA"
	}
	#Gly
	if (CodonRNA5P[i] == "GGG" | CodonRNA5P[i] == "GGA" | CodonRNA5P[i] == "GGT" | CodonRNA5P[i] == "GGC"){
		CodonRNA5P[i] = "GGT" 
	}
	#Arg
	if (CodonRNA5P[i] == "CGT" | CodonRNA5P[i] == "CGC" | CodonRNA5P[i] == "CGA" | CodonRNA5P[i] == "CGG" | CodonRNA5P[i] == "AGA" | CodonRNA5P[i] == "AGG"){
		CodonRNA5P[i] = "CGT"
	}
	#His
	if (CodonRNA5P[i] == "CAT" | CodonRNA5P[i] == "CAC"){
		CodonRNA5P[i] = "CAT"
	}
	#Asp
	if (CodonRNA5P[i] == "GAT" | CodonRNA5P[i] == "GAC"){
		CodonRNA5P[i] = "GAT"
	}
	#Lys
	if (CodonRNA5P[i] == "AAA" | CodonRNA5P[i] == "AAG"){
		CodonRNA5P[i] = "AAA"
	}
	#Glu
	if (CodonRNA5P[i] == "GAG" | CodonRNA5P[i] == "GAA"){
		CodonRNA5P[i] = "GAA"
	}
	#Ser TCT, TCC, TCA, TCG, AGT, AGC
	if (CodonRNA5P[i] == "AGT" | CodonRNA5P[i] == "AGC" | CodonRNA5P[i] == "TCT" | CodonRNA5P[i] == "TCC" | CodonRNA5P[i] == "TCA" | CodonRNA5P[i] == "TCG"){
		CodonRNA5P[i] = "AGT"
	}
	#Thr 	#Thr ACT, ACC, ACA, ACG
	if (CodonRNA5P[i] == "ACT" | CodonRNA5P[i] == "ACC" | CodonRNA5P[i] == "ACA" | CodonRNA5P[i] == "ACG"){
		CodonRNA5P[i] = "ACG"
	}
	#Asn
	if (CodonRNA5P[i] == "AAT" | CodonRNA5P[i] == "AAC"){
		CodonRNA5P[i] = "AAT"
	}

	#Gln
	if (CodonRNA5P[i] == "CAA" | CodonRNA5P[i] == "CAG"){
		CodonRNA5P[i] = "CAA"
	}
	#Cys
	if (CodonRNA5P[i] == "TGT" | CodonRNA5P[i] == "TGC"){
		CodonRNA5P[i] = "TGT"
	}
	#Pro
	if (CodonRNA5P[i] == "CCT" | CodonRNA5P[i] == "CCC" | CodonRNA5P[i] == "CCA" | CodonRNA5P[i] == "CCG"){
		CodonRNA5P[i] = "CCG"
	}
	#Val
	if (CodonRNA5P[i] == "GTG" | CodonRNA5P[i] == "GTA" | CodonRNA5P[i] == "GTT" | CodonRNA5P[i] == "GTC"){
		CodonRNA5P[i] = "GTT"
	}
	#Ile
	if (CodonRNA5P[i] == "ATT" | CodonRNA5P[i] == "ATC" | CodonRNA5P[i] == "ATA"){
		CodonRNA5P[i] = "ATT"
	}
	#Leu
	if (CodonRNA5P[i] == "TTA" | CodonRNA5P[i] == "TTG" | CodonRNA5P[i] == "CTT" | CodonRNA5P[i] == "CTC" | CodonRNA5P[i] == "CTA" | CodonRNA5P[i] == "CTG"){
		CodonRNA5P[i] = "CTG"
	}
	#Phe
	if (CodonRNA5P[i] == "TTT" | CodonRNA5P[i] == "TTC"){
		CodonRNA5P[i] = "TTT"
	}
	#Tyr
	if (CodonRNA5P[i] == "TAT" | CodonRNA5P[i] == "TAC"){
		CodonRNA5P[i] = "TAT"
	}
	#Trp
	if (CodonRNA5P[i] == "TGG"){
		CodonRNA5P[i] = "TGG"
	}
	i = i + 1;
}
#------------------------------------------------------------------------

#------------------------------------------------------------------------
OpSeq = paste(CodonRNA5P,collapse="")
write(OpSeq, file = "5POnly_Opt.txt", sep = "")
#------------------------------------------------------------------------

#Loop Conditions Higher GC% at the end---------------------------------
CodonRNA3P = CodonRNA5P
LengthOfLoop = 16;
StartCodon = length(CodonRNA3P) - 16
i = StartCodon;#Start on the second codon as the first should be start codon

for (i in 1:LengthOfLoop){
	#Ala
	if (CodonRNA3P[i] == "GCG" | CodonRNA3P[i] == "GCA" | CodonRNA3P[i] == "GCT" | CodonRNA3P[i] == "GCC"){
		CodonRNA3P[i] = "GCC"
	}
	#Gly
	if (CodonRNA3P[i] == "GGG" | CodonRNA3P[i] == "GGA" | CodonRNA3P[i] == "GGT" | CodonRNA3P[i] == "GGC"){
		CodonRNA3P[i] = "GGC" 
	}
	#Arg
	if (CodonRNA3P[i] == "CGT" | CodonRNA3P[i] == "CGC" | CodonRNA3P[i] == "CGA" | CodonRNA3P[i] == "CGG" | CodonRNA3P[i] == "AGA" | CodonRNA3P[i] == "AGG"){
		CodonRNA3P[i] = "CGT"
	}
	#His
	if (CodonRNA3P[i] == "CAT" | CodonRNA3P[i] == "CAC"){
		CodonRNA3P[i] = "CAC"
	}
	#Asp
	if (CodonRNA3P[i] == "GAT" | CodonRNA3P[i] == "GAC"){
		CodonRNA3P[i] = "GAC"
	}
	#Lys
	if (CodonRNA3P[i] == "AAA" | CodonRNA3P[i] == "AAG"){
		CodonRNA3P[i] = "AAG"
	}
	#Glu
	if (CodonRNA3P[i] == "GAG" | CodonRNA3P[i] == "GAA"){
		CodonRNA3P[i] = "GAG"
	}
	#Ser TCT, TCC, TCA, TCG, AGT, AGC
	if (CodonRNA3P[i] == "AGT" | CodonRNA3P[i] == "AGC" | CodonRNA3P[i] == "TCT" | CodonRNA3P[i] == "TCC" | CodonRNA3P[i] == "TCA" | CodonRNA3P[i] == "TCG"){
		CodonRNA3P[i] = "TCG"
	}
	#Thr ACT, ACC, ACA, ACG
	if (CodonRNA3P[i] == "ACT" | CodonRNA3P[i] == "ACC" | CodonRNA3P[i] == "ACA" | CodonRNA3P[i] == "ACG"){
		CodonRNA3P[i] = "ACG"
	}
	#Asn
	if (CodonRNA3P[i] == "AAT" | CodonRNA3P[i] == "AAC"){
		CodonRNA3P[i] = "AAC"
	}

	#Gln
	if (CodonRNA3P[i] == "CAA" | CodonRNA3P[i] == "CAG"){
		CodonRNA3P[i] = "CAG"
	}
	#Cys
	if (CodonRNA3P[i] == "TGT" | CodonRNA3P[i] == "TGC"){
		CodonRNA3P[i] = "TGC"
	}
	#Pro
	if (CodonRNA3P[i] == "CCT" | CodonRNA3P[i] == "CCC" | CodonRNA3P[i] == "CCA" | CodonRNA3P[i] == "CCG"){
		CodonRNA3P[i] = "CCG"
	}
	#Val
	if (CodonRNA3P[i] == "GTG" | CodonRNA3P[i] == "GTA" | CodonRNA3P[i] == "GTT" | CodonRNA3P[i] == "GTC"){
		CodonRNA3P[i] = "GTG"
	}
	#Ile
	if (CodonRNA3P[i] == "ATT" | CodonRNA3P[i] == "ATC" | CodonRNA3P[i] == "ATA"){
		CodonRNA3P[i] = "ATC"
	}
	#Leu
	if (CodonRNA3P[i] == "TTA" | CodonRNA3P[i] == "TTG" | CodonRNA3P[i] == "CTT" | CodonRNA3P[i] == "CTC" | CodonRNA3P[i] == "CTA" | CodonRNA3P[i] == "CTG"){
		CodonRNA3P[i] = "CTG"
	}
	#Phe
	if (CodonRNA3P[i] == "TTT" | CodonRNA3P[i] == "TTC"){
		CodonRNA3P[i] = "TTC"
	}
	#Tyr
	if (CodonRNA3P[i] == "TAT" | CodonRNA3P[i] == "TAC"){
		CodonRNA3P[i] = "TAC"
	}
	#Trp
	if (CodonRNA3P[i] == "TGG"){
		CodonRNA3P[i] = "TGG"
	}
	i = i + 1;
}
#------------------------------------------------------------------------

#Concatenante optimized codon vector & output into work folder----------- 
OpSeq = paste(CodonRNA3P,collapse="")
write(OpSeq, file = "CodonRNA3P_5P.txt", sep = "")
#------------------------------------------------------------------------



#Vipul Madahar 2018/2019
