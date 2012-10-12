Geometric Analysis of Structural Variants (GASV)

GASV is software for identifying and structural variants (SVs) from paired-end sequencing data. 

GASV
Version: 2.0
Version Date: May 30, 2012

GASVPro
Version: 1.2
Version Date: October 12, 2012

Contact: Benjamin Raphael at gasv@cs.brown.edu

If you use GASV in your research, please cite: 
S. Sindi, E. Helman, A. Bashir, B.J. Raphael. A Geometric Approach for
Classification and Comparison of Structural Variants. Joint: 17th
Annual International Conference on Intelligent Systems in Molecular
Biology and 8th Annual International European Conference on
Computational Biology (ISMB/ECCB 09).  Bioinformatics 2009
25(12):i222-i230. 

If you use GASVPro in your research, please cite:
S. Sindi, S. Onal, L. Peng, H.T. Wu, B.J Raphael. An Integrative Probabilistic 
Model for Identification of Structural Variation in Sequencing Data. 
Genome Biology 2012 27;13(3). 

SUMMARY ==============================================

This README file lists and annotates the files in the
svn repository hosted at

<http://code.google.com/p/gasv>

A read-only copy can be checked out 

% svn checkout http://gasv.googlecode.com/svn/trunk/ gasv-read-only

CONTENTS =============================================

(i) Installation Script

	To install from the command line, simply type:
	
	% ./install
	
	Installation creates the following:
	
	bin/GASV.jar  
	bin/BAMToGASV.jar 
	bin/BAMToGASV_AMBIG.jar
	bin/GASVPro-HQ.sh 
	bin/GASVPro.sh
	bin/GASVPro-CC		
	bin/GASVPro-graph
	bin/GASVPro-mcmc

	Note that GASV and BAMToGASV are maintained in Java. 
	Executable jar files are built using ant, available 
	at:
	
	<http://ant.apache.org/>

(ii) Documentation (doc/ subdirectory):

* GASV_UserGuide.pdf
	The complete GASV manual including a "Quick Start" 
	introduction from an example BAM file "Example.bam" 
	to a set of predicted structural variants. 

* GASV Release Notes; RELEASE_NOTES.txt
	Documentation of features, options developed with each GASV 
	Release.

* Software License Information; LICENSE.txt

(iii) Software (in src/)

* Source code for GASV.jar is in src/gasv

* Source code for BAMToGASV.jar and BAMToGASV_AMBIG.jar 
  is in src/bamtogasv

* Source code for GASVPro is in src/gasvPro

* Running ./install builds executable GASV.jar and BAMToGASV.jar in bin/ 

(iv) Other Dependencies

* jarfiles/ subdirectory:
	Contains jar files from Picard:
	- FixMateInformation.jar	
	- sam-1.68.jar 
		
	BAMToGASV.jar re-packages these jar files.	

	For more information visit:
       <http://picard.sourceforge.net/>	
	
* scripts/ subdirectory:
     - sortPR.bash 
	    If you create your own GASV input files, without using
		BAMToGASV.jar, the paired-read files must be sorted with
		sortPR.bash prior to running GASV.jar.

     - GASVPruneClusters.pl
	   Code used by GASVPro.sh and GASVPro-HQ.sh for pruning redundant
           predictions from set of clusters.

