/**
 * Copyright 2010 Benjamin Raphael, Suzanne Sindi, Hsin-Ta Wu, Anna Ritz, Luke Peng
 *
 *  This file is part of gasv.
 * 
 *  gasv is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  gasv is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with gasv.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package maqpreprocessor;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;



/**
 * Takes a MAP file and an Lmin and an Lmax partitions the file into valid and invalid subsets.
 * 
 *
 * @author Anna Ritz
 * @date 1/14/09
 *
 */
public class ExtractInvalidPairs {

	private String prefix;
	private BufferedWriter rawreads_invs,rawreads_dels,rawreads_divs,rawreads_trans;

	public ExtractInvalidPairs(String p,boolean method) throws IOException {
		prefix = p;
		if (method) {
		rawreads_invs = new BufferedWriter(new FileWriter(p+"_inversions.reads"));
		rawreads_dels = new BufferedWriter(new FileWriter(p+"_deletions.reads"));
		rawreads_divs = new BufferedWriter(new FileWriter(p+"_divergent.reads"));
		} else {
			rawreads_trans = new BufferedWriter(new FileWriter(p+"_translocations.reads"));
		}
		
	}

	/**
	 * main method.  There are two ways to call this method.
	 * A: Generate invalid files for inversions, convergent and divergent pairs
	 * B: Generate invalid files for pairs whose ends are on Chromosomes I and J.
	 * 
	 * A PARAMS:
	 * @param args - vector of 
	 * (0) true for Method A, false for Method B
	 * (1)an *.txt file produced by
	 * % maq mapview in.map > out.txt
	 * (2)prefix for other output files
	 * (3) Lmin
	 * (4) Lmax
	 * (5) min inversion length
	 * (6) min deletion length
	 * (7) whether we print the diffchrom file.
	 * 
	 * B PARAMS:
	 * @param args - vector of 
	 * (0) true for Method A, false for Method B
	 * (1)an *.txt file produced by
	 * % maq mapview in.map > out.txt
	 * (2) an *.txt file produced by
	 * % maq mapview in.map > out.txt
	 * (3)prefix for other output files
	 */
	public static void main(String[] args) {
		long start = System.currentTimeMillis();
		// put input file names into an arraylist
		try {
			boolean method = new Boolean(args[0]).booleanValue();
			if (method) { // Method A
				System.out.println("Extract inversions and deletions");
				ExtractInvalidPairs eip = new ExtractInvalidPairs(args[2],method);

				boolean printDiffChrom = false;
				if (args[7].equals("1"))
					printDiffChrom = true;

				// use Lmin and Lmax to generate .invalid files.
				eip.printInvalids(args[1],new Integer(args[3]).intValue(),
						new Integer(args[4]).intValue(),
						new Integer(args[5]).intValue(),
						new Integer(args[6]).intValue(),printDiffChrom);
				eip.rawreads_invs.close();
				eip.rawreads_dels.close();
				eip.rawreads_divs.close();
			} else if (!new Boolean(args[0]).booleanValue()) { // Method B
				ExtractInvalidPairs eip = new ExtractInvalidPairs(args[3],method);

				// read diff chrom reads to generate diffChrom.invalid file
				System.out.println("Printing translocation invalid file");
				eip.printDiffChromInvalids(args[1],args[2]);
				eip.rawreads_trans.close();
			} else
				throw new Exception("Incorrect number of parameters.");

			// catch any exception: FileNotFoundException, IOException
		} catch (Exception e) {
			e.printStackTrace();
		}

		long end = System.currentTimeMillis();
		System.out.println("time taken: " + (end-start));
	}


	/**
	 * Takes a file and outputs five different files:
	 * (2) file of reads flagged as inversions > 10KB
	 * (3) file of reads flagged as divergent (RF) pairs
	 * (4) file of reads flagged as convergent (FR) pairs smaller than
	 * Lmin or greater than Lmax
	 * (5) file of reads flagged as deletions > 5KB
	 * 
	 * @param filename - data
	 * @param Lmin - min valid read length
	 * @param Lmax - max valid read length
	 * @param printDiffChrom 
	 * 
	 * @throws IOException
	 */
	public void printInvalids(String filename,int Lmin,int Lmax,int minInvLen, int minDelLen, boolean printDiffChrom) throws IOException {
		MAPCursor cursor = new MAPCursor(filename);

		BufferedWriter inversions = new BufferedWriter(new FileWriter(prefix+"_inversions.invalid"));
		BufferedWriter deletions = new BufferedWriter(new FileWriter(prefix+"_deletions.invalid"));
		BufferedWriter divergent = new BufferedWriter(new FileWriter(prefix+"_divergent.invalid"));

		BufferedWriter diffchrom = null;
		if (printDiffChrom) 
			diffchrom = new BufferedWriter(new FileWriter(filename+".translocations"));
		
		HashMap<String,MaqRead> invPairs  = new HashMap<String,MaqRead>();
		HashMap<String,MaqRead> delPairs  = new HashMap<String,MaqRead>();
		HashMap<String,MaqRead> divPairs  = new HashMap<String,MaqRead>();
		
		int counter =0,invcounter =0,divcounter=0,delcounter = 0;
		MaqRead read;
		while(cursor.hasNext()) {
			cursor.advance();
			
			if (printDiffChrom && cursor.diffChrom())
				diffchrom.write(cursor.line+"\n");
			
			// DEBUG
			//cursor.printLine();
			
			counter++;
			if (counter % 100000 == 0) {
				System.out.println("\tLine " + counter +": "+ invcounter + " inversion fragments, " +
						delcounter + " deletion fragments, " + 
						invPairs.size()+ " unpaired inversions, " + delPairs.size() + " unpaired deletions.");
			}

			// ensure that the read is unique (non-zero mapping quality) and
			// that the quality is high-scoring (above 30)
			if (cursor.mappingquality >= 30) {
				
				//System.out.println("DEBUG: passing mapping quality");
				if (cursor.convergent() && cursor.pairLen >= minDelLen) {
					
					// insertions/deletions on the same chromosome
					if (delPairs.containsKey(cursor.readname)) {
						read = delPairs.get(cursor.readname);
						
						if (cursor.strand=='+') 
							writeline(deletions,cursor.readname,cursor.genOutput(),read.outputstring);
						else
							writeline(deletions,cursor.readname,read.outputstring,cursor.genOutput());
						rawreads_dels.write(cursor.line+"\n");
						rawreads_dels.write(read.line+"\n");
						delcounter++;
						delPairs.remove(cursor.readname);
					} else
						delPairs.put(cursor.readname,new MaqRead(cursor.genOutput(),cursor.startposition,cursor.line));
				} else if (cursor.inversion() && cursor.pairLen >= minInvLen) {
					
					// inversions on the same chromosome
					if (invPairs.containsKey(cursor.readname)) {
						read = invPairs.get(cursor.readname);
						if (read.startposition < cursor.startposition) 
							writeline(inversions,cursor.readname,read.outputstring,cursor.genOutput());
						else
							writeline(inversions,cursor.readname,cursor.genOutput(),read.outputstring);
						rawreads_invs.write(cursor.line+"\n");
						rawreads_invs.write(read.line+"\n");
						invcounter++;
						invPairs.remove(cursor.readname);
					} else
						invPairs.put(cursor.readname,new MaqRead(cursor.genOutput(),cursor.startposition,cursor.line));
				} else if (cursor.divergent()) {
					
					// more complicated, RF oriented pairs on the same chromosome
					if (divPairs.containsKey(cursor.readname)) {
						read = divPairs.get(cursor.readname);
						if (read.startposition < cursor.startposition) 
							writeline(divergent,cursor.readname,read.outputstring,cursor.genOutput());
						else
							writeline(divergent,cursor.readname,cursor.genOutput(),read.outputstring);
						rawreads_divs.write(cursor.line+"\n");
						rawreads_divs.write(read.line+"\n");
						divcounter++;
						divPairs.remove(cursor.readname);
					} else
						divPairs.put(cursor.readname,new MaqRead(cursor.genOutput(),cursor.startposition,cursor.line));
				}
			} 
		}

		System.out.println("\tLine " + counter +": "+ invcounter + " inversion fragments, " +
				delcounter + " deletion fragments, " + 
				invPairs.size()+ " unpaired inversions, " + delPairs.size() + " unpaired deletions.");
		// close writers
		inversions.close();
		divergent.close();
		deletions.close();
		if (printDiffChrom)
			diffchrom.close();
	}

	public void writeline(BufferedWriter writer,String name, String partA,String partB) throws IOException {
		writer.write(name+"\t"+partA+"\t"+partB+"\n");
	}
	/**
	 * Takes two files and reports the pairs whose ends span the two chromosomes.
	 * @throws IOException 
	 */
	public void printDiffChromInvalids(String fname1, String fname2) throws IOException {

		HashMap<String,MaqRead> pairs  = new HashMap<String,MaqRead>();
		int counter = 0;

		// take first file, report all the pairs with ends on different chromosomes.
		MAPCursor cursor = new MAPCursor(fname1);
		while(cursor.hasNext()) {
			cursor.advance();

			counter++;
			if (counter % 1000000 == 0) {
				System.out.println("\tFile 1: Line " + counter +": " 
						+ pairs.size() + " reads with pairs on different chromosomes");
				// testing
				//break;
			}

			if (cursor.mappingquality >= 30 && cursor.diffChrom()) {
				pairs.put(cursor.readname,new MaqRead(cursor.genOutput(),cursor.startposition,cursor.line));
			}
		}
		cursor.close();

		BufferedWriter writer = new BufferedWriter(new FileWriter(prefix+"_translocations.invalid"));

		// take second file, write all completed pairs to file.
		counter = 0;
		cursor = new MAPCursor(fname2);
		MaqRead read;
		while(cursor.hasNext()) {
			cursor.advance();

			counter++;
			if (counter % 1000000 == 0) {
				System.out.println("\tFile 2: Line " + counter +": " 
						+ pairs.size() + " reads with pairs on different chromosomes");
			}

			if (cursor.mappingquality >= 30 && cursor.diffChrom() && pairs.containsKey(cursor.readname)) {
				read = pairs.get(cursor.readname);
				writeline(writer,cursor.readname,read.outputstring,cursor.genOutput());
				rawreads_trans.write(cursor.line+"\n");
				rawreads_trans.write(read.line+"\n");
				pairs.remove(cursor.readname);
			}
		}
		cursor.close();

		writer.close();
	}
}
