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

Version: 2.0
Version Date: May 30, 2012

SUMMARY ==============================================

This README file lists and annotates the files in the
svn repository hosted at

<http://code.google.com/p/gasv>

A read-only copy can be checked out 

% svn checkout http://gasv.googlecode.com/svn/trunk/ gasv-read-only

CONTENTS =============================================

(i) Installation Script (build.xml)

	The GASV software package is maintained in Java. 
	Executable jar files can be built usint ant, available 
	at:
	
	<http://ant.apache.org/>
	
	To install from the command line, simply type:
	
	% ant
	
	Installation creates the following:
	
	bin/GASV.jar  
	bin/BAMToGASV.jar 
	
(ii) Documentation (doc/ subdirectory):

* GASV_UserGuide.pdf
	The complete GASV manual including a "Quick Start" 
	introduction from an example BAM file "Example.bam" 
	to a set of predicted structural variants. 

* GASV Release Notes; RELEASE_NOTES.txt
	Documentation of features, options developed with each GASV 
	Release.

* Software License Information; LICENSE.txt

(iii) GASV Software (in src/)

* Source code for GASV.jar is in src/gasv

* Source code for BAMToGASV.jar is in src/bamtogasv

* Run ant to build executable GASV.jar and BAMToGASV.jar in bin/ 

(iv) Other Dependencies

* jarfiles/ subdirectory:
	Contains:
	- FixMateInformation.jar	
	- sam-1.68.jar 
	
	BAMToGASV.jar re-packages these two jar files from Picard.
	For more information visit:
	
    <http://picard.sourceforge.net/>	
	
* scripts/ subdirectory:
     - sortPR.bash 
	    If you create your own GASV input files, without using
		BAMToGASV.jar, the paired-read files must be sorted with
		sortPR.bash prior to running GASV.jar.



