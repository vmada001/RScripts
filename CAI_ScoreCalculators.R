
#Downloaded and Installed Libraries-----------------------------------------------------
library(stringr);
library(seqinr);
#---------------------------------------------------------------------------------------

#Work Directory-------------------------------------------------------------------------
setwd("C:/Users/vipul/OneDrive/Documents/Research/CodonOptimization/RNF4");
#---------------------------------------------------------------------------------------

##Read in DNA Sequences------------------------------------------------
DNASeqRef = toupper(readLines("RNF4.txt"));
##---------------------------------------------------------------------

##Make a DNA Vector-----------------------------------------------------
DNAVectorSt <- s2c(paste(DNASeqRef, collapse = ""));
##----------------------------------------------------------------------

##Make a Vector of Codons-----------------------------------------------
CodonSeqSt = splitseq(DNAVectorSt, frame = 0, word = 3);
CodonSeqXX = CodonSeqSt;
##-----------------------------------------------------------------------


#Loop Conditions CAI Calculations-----------------------------------------
CodonCAI = 0;
CodonCAIFinal = 0;
CodonSeqRef = CodonSeqSt;
LengthOfLoop <- length(CodonSeqRef);
i = 1;#Start on the second codon as the first should be start codon

for (i in 1:LengthOfLoop){
	#Start/Stop Codon
	if (CodonSeqRef[i] == "ATG" | CodonSeqRef[i] == "TGA" | CodonSeqRef[i] == "TAA"| CodonSeqRef[i] == "TAG"){CodonCAI[i] = 1.000}

	#Ala
	if (CodonSeqRef[i] == "GCG"){CodonCAI[i] = 1.000}
      if (CodonSeqRef[i] == "GCA"){CodonCAI[i] = 0.601}
	if (CodonSeqRef[i] == "GCT"){CodonCAI[i] = 0.423}
	if (CodonSeqRef[i] == "GCC"){CodonCAI[i] = 0.781}

	#Gly
	if (CodonSeqRef[i] == "GGG"){CodonCAI[i] = 0.432}
	if (CodonSeqRef[i] == "GGA"){CodonCAI[i] = 0.248}
	if (CodonSeqRef[i] == "GGT"){CodonCAI[i] = 0.736}
	if (CodonSeqRef[i] == "GGC"){CodonCAI[i] = 1.000}

	#Pro
	if (CodonSeqRef[i] == "CCT"){CodonCAI[i] = 0.234}
	if (CodonSeqRef[i] == "CCC"){CodonCAI[i] = 0.096}
	if (CodonSeqRef[i] == "CCA"){CodonCAI[i] = 0.298}
	if (CodonSeqRef[i] == "CCG"){CodonCAI[i] = 1.000}

	#Arg
	if (CodonSeqRef[i] == "CGT"){CodonCAI[i] = 0.873}
	if (CodonSeqRef[i] == "CGC"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "CGA"){CodonCAI[i] = 0.127}
	if (CodonSeqRef[i] == "CGG"){CodonCAI[i] = 0.268}
	if (CodonSeqRef[i] == "AGA"){CodonCAI[i] = 0.127}
	if (CodonSeqRef[i] == "AGG"){CodonCAI[i] = 0.113}
	
	#His
	if (CodonSeqRef[i] == "CAT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "CAC"){CodonCAI[i] = 0.772}

	#Lys
	if (CodonSeqRef[i] == "AAA"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "AAG"){CodonCAI[i] = 0.300}

	#Glu
	if (CodonSeqRef[i] == "GAG"){CodonCAI[i] = 0.612}
	if (CodonSeqRef[i] == "GAA"){CodonCAI[i] = 1.000}

	#Ser
	if (CodonSeqRef[i] == "AGT"){CodonCAI[i] = 0.630}
	if (CodonSeqRef[i] == "AGC"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "TCT"){CodonCAI[i] = 0.593}
	if (CodonSeqRef[i] == "TCC"){CodonCAI[i] = 0.556}
	if (CodonSeqRef[i] == "TCA"){CodonCAI[i] = 0.426}
	if (CodonSeqRef[i] == "TCG"){CodonCAI[i] = 0.796}

	#Asp
	if (CodonSeqRef[i] == "GAT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "GAC"){CodonCAI[i] = 0.505}

	#Thr
	if (CodonSeqRef[i] == "ACT"){CodonCAI[i] = 0.305} 
	if (CodonSeqRef[i] == "ACC"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "ACA"){CodonCAI[i] = 0.242}
	if (CodonSeqRef[i] == "ACG"){CodonCAI[i] = 0.579}

	#Asn
	if (CodonSeqRef[i] == "AAT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "AAC"){CodonCAI[i] = 0.750}

	#Gln
	if (CodonSeqRef[i] == "CAA"){CodonCAI[i] = 0.548}
	if (CodonSeqRef[i] == "CAG"){CodonCAI[i] = 1.000}

	#Cys
	if (CodonSeqRef[i] == "TGT"){CodonCAI[i] = 0.727}
	if (CodonSeqRef[i] == "TGC"){CodonCAI[i] = 1.000}

	#Val
	if (CodonSeqRef[i] == "GTG"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "GTA"){CodonCAI[i] = 0.313}
	if (CodonSeqRef[i] == "GTT"){CodonCAI[i] = 0.578}
	if (CodonSeqRef[i] == "GTC"){CodonCAI[i] = 0.422}	

	#Ile
	if (CodonSeqRef[i] == "ATT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "ATC"){CodonCAI[i] = 0.914}
	if (CodonSeqRef[i] == "ATA"){CodonCAI[i] = 0.148}

	#Leu
	if (CodonSeqRef[i] == "TTA"){CodonCAI[i] = 0.319}
	if (CodonSeqRef[i] == "TTG"){CodonCAI[i] = 0.338}
	if (CodonSeqRef[i] == "CTT"){CodonCAI[i] = 0.232}
	if (CodonSeqRef[i] == "CTC"){CodonCAI[i] = 0.256}
	if (CodonSeqRef[i] == "CTA"){CodonCAI[i] = 0.063}
	if (CodonSeqRef[i] == "CTG"){CodonCAI[i] = 1.000}

	#Phe
	if (CodonSeqRef[i] == "TTT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "TTC"){CodonCAI[i] = 0.652}

	#Tyr
	if (CodonSeqRef[i] == "TAT"){CodonCAI[i] = 1.000}
	if (CodonSeqRef[i] == "TAC"){CodonCAI[i] = 0.457}

	#Trp
	if (CodonSeqRef[i] == "TGG"){CodonCAI[i] = 1.000}
	i = i + 1;
}

GeneCAI = (prod(CodonCAI))^(1/length(CodonCAI));

CodonCAI[1+length(CodonCAI)] = GeneCAI;

#Concatenante optimized codon vector & output into work folder------------------------------------------------------ 
write(CodonCAI, file = "CodonCAIwt.txt", sep = "\n")
#-------------------------------------------------------------------------------------------------------------------

CodonCAI

CAI_GM = CodonCAI[1];
Count = 1
Codon = 0
i = 1;
CAIAverage = 0;
for(i in 2:length(CodonCAI)){
	CAI_GM = CAI_GM*CodonCAI[i];
	Codon = Codon + 1;
	if(Codon == 9){
		Codon = 0;
		CAIAverage[Count] = ((CAI_GM)^(1/9));
		CAI_GM = CodonCAI[i];
		Count = Count + 1;
	}
	i = i+1;
}
CAIAverage
#Concatenante optimized codon vector & output into work folder------------------------------------------------------ 
write(CAIAverage, file = "CAI_GM.txt", sep = "\n")
#-------------------------------------------------------------------------------------------------------------------

CAIAverage
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

	#Lys Check
	if (CodonSeqSt[i] == "AAG"){
		b = i - 1;
		if(CodonSeqSt[b] == "AAA" | CodonSeqSt[b] == "AAG"){
			CodonSeqXX[i] = "AAG"
			CodonSeqXX[b] = "AAG" 
		}
	}

	#Glu Check
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
	#Ser AG Check
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
write(OpSeq, file = "PD1_01222019.txt", sep = "")
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
write(OpSeq, file = "PD1_5P_01232019.txt", sep = "")
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
write(OpSeq, file = "PD1_RNA3p_5p_09032018.txt", sep = "")
#------------------------------------------------------------------------




#Concatenante optimized codon vector & output into work folder----------- 
OpSeq = paste(CodonRNA5P,collapse="")
RNA5pSeq = paste(CodonRNA5P,collapse="")
RNA3pSeq = paste(CodonRNA5P,collapse="")
write(OpSeq, file = "Codon_RNA3p09032018.txt", sep = "")
#------------------------------------------------------------------------


#RNAFOLD Analysis-------------------------------------------------------
#Loop Conditions--------------------------------------------------------
Open = 0
Close = 0
CodonWindow = 20
looplength = (length(CodonSeqSt)/CodonWindow)
for(i in 1:looplength){
	CTSeq = paste(CodonSeqSt[10:20],collapse="")
	CTRNASeq = paste(CodonRNA5P[1:16],collapse="")
	RNAfold.path = '"C:/Program Files (x86)/ViennaRNA Package/RNAfold.exe"'
	Test = LncFinder::run_RNAfold(CTSeq, RNAfold.path=RNAfold.path, parallel.cores = 2)
	TestVector <- s2c(paste(Test, collapse = " "));
	MFE = TestVector[(length(TestVector)-6):(length(TestVector)- 2)]
	if(MFE[1] == "\""){MFE = MFE[2:5]}
		RNAMFE = as.numeric(paste(MFE,collapse=""))
}

CTSeq = paste(CodonSeqSt[10:20],collapse="")
CTRNASeq = paste(CodonRNA5P[1:16],collapse="")
RNAfold.path = '"C:/Program Files (x86)/ViennaRNA Package/RNAfold.exe"'
Test = LncFinder::run_RNAfold(CTSeq, RNAfold.path=RNAfold.path, parallel.cores = 2)
TestVector <- s2c(paste(Test, collapse = " "));
MFE = TestVector[(length(TestVector)-6):(length(TestVector)- 2)]
if(MFE[1] == "\""){MFE = MFE[2:5]}
RNAMFE = as.numeric(paste(MFE,collapse=""))
#------------------------------------------------------------------------




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

