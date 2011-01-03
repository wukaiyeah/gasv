Geometric Analysis of Structural Variants (GASV)

Software for identifying and comparing structural variants by computing 
intersections of breakpoint regions.

If you use GASV in your research, please cite: 
S. Sindi, E. Helman, A. Bashir, B.J. Raphael. A Geometric Approach for
Classification and Comparison of Structural Variants. Joint: 17th
Annual International Conference on Intelligent Systems in Molecular
Biology and 8th Annual International European Conference on
Computational Biology (ISMB/ECCB 09).  Bioinformatics 2009
25(12):i222-i230. 

This product includes software developed by the Solution Engineering,
Inc. (http://www.seisw.com/). 

Contact: Ben Raphael at gasv@cs.brown.edu

Version: 1.5.1
Version Date: Jan 3 2011

SUMMARY ==============================================

This README file lists and annotates the files distributed
with GASV Release 1.5.1, GASV_RELEASE_1.5.1.tar.gz.

CONTENTS =============================================

(i) Documentation (under doc/ subdirectory):

* Quickstart Guide; GASV_Quickstart_Guide.txt
	A quick introduction to GASV going from the 
	example BAM file "Example.bam" to a set of predicted
	structural variants. The Example.bam is available
	as a separate download.

* GASV Manual; MANUAL.txt	
	A complete description of GASV features and options.

* BAM Preprocessor Documentation; BAM_preprocessor.txt
	Detailed description of BAM_preprocessor.pl, used
	to convert from the BAM format to GASV input format.

* GASV Release Notes; RELEASE_NOTES.txt
	Documentation of features, options developed with each GASV 
	Release.

* GASV FAQ; FAQ.txt
	A list of "Frequently Asked Questions" regarding GASV.

* GPL license info; COPYING.txt
	A copy of the GPL license.

(ii) GASV Software (in src/java/gasv)
* Source code in src/java/gasv
* Run make_jar_cmd (unix-like OS only) to build executable gasv.jar in lib/ 

(iii) Preprocessors to convert from BAM/MAQ to GASV input files:
	- BAM Preprocessor:
		Run install.sh for full functionality (LminLmaxProcessor).
		Executable: bin/BAM_preprocessor.pl, generate_GASV.pl,
		Other source code: 	src/cpp/LminLmaxProcessor_Source/ 
					bin/generate_GASV.pl
	- MAQ Preprocessor:
		Source code in src/java/maqpreprocessor/
		Run make_jar_cmd (unix-like OS only) to build executable maqpreprocessor.jar in lib/ 

(iv) Sample Files (under sample/ subdirectory):

* Sample ESP Files: Cancer, sample 

* Sample CGH file: sampleCGH 

* Sample Batch Input Files: sample.in, Cancer.in, sampleCGH.in



