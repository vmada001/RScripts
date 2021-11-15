## Script Objective: Codon BIAS and AT/GC Content at 3' and 5' of RNA

#Downloaded and Installed Libraries-----------------------------------------------------
library(stringr);
library(seqinr);
#---------------------------------------------------------------------------------------

#Work Directory-------------------------------------------------------------------------
setwd("C:/Users/vipul/OneDrive/Documents/Research/CodonOptimization/PIAS1");
#---------------------------------------------------------------------------------------

##Read in DNA Sequences------------------------------------------------
DNASeqRef = toupper(readLines("PIAS1wt.txt"));
##---------------------------------------------------------------------

##Make a DNA Vector-----------------------------------------------------
DNAVectorSt <- s2c(paste(DNASeqRef, collapse = ""));
##----------------------------------------------------------------------

##Make a Vector of Codons-----------------------------------------------
CodonSeqSt = splitseq(DNAVectorSt, frame = 0, word = 3);
CodonSeqXX = CodonSeqSt;
##-----------------------------------------------------------------------


##Codon Count for harmonization------------------------------------------
#Loop Conditions----------------------
LengthOfLoop = length(CodonSeqXX);
i = 2;#Start on the second codon as the first should be start codon
ALA = 0
GLY = 0
ARG = 0
HIS = 0
ASP = 0
GLU = 0
SER = 0
THR = 0 
ASN = 0
GLN = 0
CYS = 0 
PRO = 0 
VAL = 0 
ILE = 0
LEU = 0
PHE = 0 
TYR = 0 
TRP = 0 
LYS = 0 

for (i in 2:LengthOfLoop){

	#Ala
	if (CodonSeqSt[i] == "GCG" | CodonSeqSt[i] == "GCA" | CodonSeqSt[i] == "GCT" | CodonSeqSt[i] == "GCC"){
		
		ALA = ALA + 1	
	}
	#Gly
	if (CodonSeqSt[i] == "GGG" | CodonSeqSt[i] == "GGA" | CodonSeqSt[i] == "GGT" | CodonSeqSt[i] == "GGC"){
		GLY = GLY + 1
	}
	#Arg
	if (CodonSeqSt[i] == "CGT" | CodonSeqSt[i] == "CGC" | CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
		ARG = ARG + 1
	}
	#His
	if (CodonSeqSt[i] == "CAT" | CodonSeqSt[i] == "CAC"){
		HIS = HIS + 1
	}
	#Asp
	if (CodonSeqSt[i] == "GAT" | CodonSeqSt[i] == "GAC"){
		ASP = ASP + 1
	}
	#Lys
	if (CodonSeqSt[i] == "AAA" | CodonSeqSt[i] == "AAG"){
		LYS = LYS + 1
	}
	#Glu
	if (CodonSeqSt[i] == "GAG" | CodonSeqSt[i] == "GAA"){
		GLU = GLU + 1
	}
	#Ser
	if (CodonSeqSt[i] == "AGT" | CodonSeqSt[i] == "AGC" | CodonSeqSt[i] == "TCT" | CodonSeqSt[i] == "TCC" | CodonSeqSt[i] == "TCA" | CodonSeqSt[i] == "TCG"){
		SER = SER + 1
	}
	#Thr
	if (CodonSeqSt[i] == "ACT" | CodonSeqSt[i] == "ACC" | CodonSeqSt[i] == "ACA" | CodonSeqSt[i] == "ACG"){
		THR = THR + 1	
	}
	#Asn
	if (CodonSeqSt[i] == "AAT" | CodonSeqSt[i] == "AAC"){
		ASN = ASN + 1
	}
	#Gln
	if (CodonSeqSt[i] == "CAA" | CodonSeqSt[i] == "CAG"){
		GLN = GLN + 1
	}
	#Cys
	if (CodonSeqSt[i] == "TGT" | CodonSeqSt[i] == "TGC"){
		CYS = CYS + 1
	}
	#Pro
	if (CodonSeqSt[i] == "CCT" | CodonSeqSt[i] == "CCC" | CodonSeqSt[i] == "CCA" | CodonSeqSt[i] == "CCG"){
		PRO = PRO + 1
	}
	#Val
	if (CodonSeqSt[i] == "GTG" | CodonSeqSt[i] == "GTA" | CodonSeqSt[i] == "GTT" | CodonSeqSt[i] == "GTC"){
		VAL = VAL + 1
	}
	#Ile
	if (CodonSeqSt[i] == "ATT" | CodonSeqSt[i] == "ATC" | CodonSeqSt[i] == "ATA"){
		ILE = ILE + 1
	}
	#Leu
	if (CodonSeqSt[i] == "TTA" | CodonSeqSt[i] == "TTG" | CodonSeqSt[i] == "CTT" | CodonSeqSt[i] == "CTC" | CodonSeqSt[i] == "CTA" | CodonSeqSt[i] == "CTG"){
		LEU = LEU + 1
	}
	#Phe
	if (CodonSeqSt[i] == "TTT" | CodonSeqSt[i] == "TTC"){
		PHE = PHE + 1
	}
	#Tyr
	if (CodonSeqSt[i] == "TAT" | CodonSeqSt[i] == "TAC"){
		TYR = TYR + 1
	}
	#Trp
	if (CodonSeqSt[i] == "TGG"){
		TRP = TRP + 1
	}
	i = i + 1;
}


##Codon Selection based on BIAS and Repeating AAs
#Loop Conditions---------------------------------------------------------
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
	if (CodonSeqSt[i] == "CAC" | CodonSeqSt[b] == "CAT"){
		b = i - 1;
		if(CodonSeqSt[b] == "CAT" | CodonSeqSt[b] == "CAC"){
			CodonSeqXX[i] == "CAC"
			CodonSeqXX[b] == "CAC" 
		}
	}

	#Gly Check
	if (CodonSeqSt[i] == "GGG" | CodonSeqSt[i] == "GGA"){
		CodonSeqSt[i] == "GGT";
		b = i - 1;
		if(CodonSeqSt[b] == "GGG" | CodonSeqSt[b] == "GGA" | CodonSeqSt[b] == "GGT" | CodonSeqSt[b] == "GGC"){
			CodonSeqXX[i] == "GGC"
			CodonSeqXX[b] == "GGC" 
		}
	}

	#His Check
	if (CodonSeqSt[i] == "CAC"){
		CodonSeqXX[i] = "CAT"
		b = i - 1;
		if(CodonSeqSt[b] == "CAT" | CodonSeqSt[b] == "CAC"){
			CodonSeqXX[i] == "CAC"
			CodonSeqXX[b] == "CAC" 
		}
	}

	#Asp Check
	if (CodonSeqSt[i] == "GAT" | CodonSeqSt[i] == "GAC"){
		b = i - 1;
		if(CodonSeqSt[b] == "GAT" | CodonSeqSt[b] == "GAC"){
			CodonSeqXX[i] == "GAC"
			CodonSeqXX[b] == "GAC" 
		}
	}

	#Lys Check
	if (CodonSeqSt[i] == "AAG"){
		b = i - 1;
		if(CodonSeqSt[b] == "AAA" | CodonSeqSt[b] == "AAG"){
			CodonSeqXX[i] == "AAG"
			CodonSeqXX[b] == "AAG" 
		}
	}

	#Glu Check
	if (CodonSeqSt[i] == "GAG" | CodonSeqSt[i] == "GAA"){
		b = i - 1;
		if(CodonSeqSt[b] == "GAG" | CodonSeqSt[b] == "GAG"){
			CodonSeqXX[i] == "GAG"
			CodonSeqXX[b] == "GAG" 
		}
	}

	#Ser TC
	if (CodonSeqSt[i] == "TCA" | CodonSeqSt[i] == "TCG"){
		CodonSeqXX[i] = "TCT"
		b = i - 1;
		if(CodonSeqSt[b] == "AGT" | CodonSeqSt[b] == "AGC" | CodonSeqSt[b] == "TCT" | CodonSeqSt[b] == "TCC" | CodonSeqSt[b] == "TCA" | CodonSeqSt[b] == "TCG"){
			CodonSeqXX[i] == "AGC"
			CodonSeqXX[b] == "AGC" 
		}
	}
	#Ser AG Check
	if (CodonSeqSt[i] == "AGT"){
		CodonSeqXX[i] = "AGC"
		b = i - 1;
		if(CodonSeqSt[b] == "AGT" | CodonSeqSt[b] == "AGC" | CodonSeqSt[b] == "TCT" | CodonSeqSt[b] == "TCC" | CodonSeqSt[b] == "TCA" | CodonSeqSt[b] == "TCG"){
			CodonSeqXX[i] == "AGC"
			CodonSeqXX[b] == "AGC" 
		}
	}

	#Arg Check
	if (CodonSeqSt[i] == "CGT" | CodonSeqSt[i] == "CGC" | CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
		if (CodonSeqSt[i] == "CGA" | CodonSeqSt[i] == "CGG" | CodonSeqSt[i] == "AGA" | CodonSeqSt[i] == "AGG"){
			CodonSeqXX[i] = "CGT"
		}
		b = i - 1
		if(CodonSeqSt[b] == "CGT" | CodonSeqSt[b] == "CGC" | CodonSeqSt[b] == "CGA" | CodonSeqSt[b] == "CGG"|CodonSeqSt[b] == "AGA" | CodonSeqSt[b] == "AGG"){
			CodonSeqXX[i] == "CGC"
			CodonSeqXX[b] == "CGC" 
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
		if (CodonSeqSt[i] == "ACA" | CodonSeqSt[i] == "ACG"){
			CodonSeqXX[i] = "ACC"
		}
		b = i - 1
		if (CodonSeqSt[b] == "ACT" | CodonSeqSt[b] == "ACC" | CodonSeqSt[b] == "ACA" | CodonSeqSt[b] == "ACG"){
			CodonSeqXX[b] = "ACT"
			CodonSeqXX[i] = "ACT"
		}
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
		b = i - 1
		if (CodonSeqSt[b] == "CCT" | CodonSeqSt[b] == "CCC" | CodonSeqSt[b] == "CCA" | CodonSeqSt[b] == "CCG"){
			CodonSeqXX[b] = "CCA"
		}
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
		if (CodonSeqSt[i] == "GCC"){
			CodonSeqXX[i] = "GCG"
		}
		b = i - 1
		if (CodonSeqSt[b] == "GCG" | CodonSeqSt[b] == "GCA" | CodonSeqSt[b] == "GCT" | CodonSeqSt[b] == "GCC"){
			CodonSeqXX[b] = "GCT"
			CodonSeqXX[i] = "GCT"
		}	
	}

	#Val Check
	if (CodonSeqSt[i] == "GTG" | CodonSeqSt[i] == "GTA" | CodonSeqSt[i] == "GTT" | CodonSeqSt[i] == "GTC"){
		if (CodonSeqSt[i] == "GTC"){
			CodonSeqXX[i] = "GTT"
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
		CodonSeqXX[i] = "TAC"
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
	if (CodonSeqSt[i] == "TGA"){
		CodonSeqXX[i] = "TAA"
	}
	i = i + 1;
}
#------------------------------------------------------------------------

#CodonFrequencyOptimized-------------------------------------------------
OpSeq = paste(CodonSeqXX,collapse="")
write(OpSeq, file = "PIAS1_01222019.txt", sep = "")
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
	#Ser
	if (CodonRNA5P[i] == "AGT" | CodonRNA5P[i] == "AGC" | CodonRNA5P[i] == "TCT" | CodonRNA5P[i] == "TCC" | CodonRNA5P[i] == "TCA" | CodonRNA5P[i] == "TCG"){
		CodonRNA5P[i] = "AGT"
	}
	#Thr
	if (CodonRNA5P[i] == "ACT" | CodonRNA5P[i] == "ACC" | CodonRNA5P[i] == "ACA" | CodonRNA5P[i] == "ACG"){
		CodonRNA5P[i] = "AGC"
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
write(OpSeq, file = "PIAS1_5P_01232019.txt", sep = "")
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
	#Ser
	if (CodonRNA3P[i] == "AGT" | CodonRNA3P[i] == "AGC" | CodonRNA3P[i] == "TCT" | CodonRNA3P[i] == "TCC" | CodonRNA3P[i] == "TCA" | CodonRNA3P[i] == "TCG"){
		CodonRNA3P[i] = "TCG"
	}
	#Thr
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
write(OpSeq, file = "PIAS1_RNA3p_5p_09032018.txt", sep = "")
#------------------------------------------------------------------------


#Concatenante optimized codon vector & output into work folder----------- 
OpSeq = paste(CodonRNA5P,collapse="")
RNA5pSeq = paste(CodonRNA5P,collapse="")
RNA3pSeq = paste(CodonRNA5P,collapse="")
write(OpSeq, file = "Codon_RNA3p09032018.txt", sep = "")
#------------------------------------------------------------------------
