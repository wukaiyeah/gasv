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
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.InputMismatchException;
import java.util.Scanner;

/**
 * MAPCursor is an iterator that reads the output produced by a .map
 * file line-by-line, storing the relevant information.
 *  **FROM maq.pdf
The output from mapview contains the following columns (from maq.pdf):
(1) Read Name
(2) Chromosome
(3) Position
(4) Strand
(5) Insert size from the outer coordinates of a pair
(6) paired flag
(7) mapping quality
(8) single-end mapping quality
(9) alternative mapping quality
(10) number of mismatches of the best hit
(11) sum of the qualities of mismatched bases of the best hit
(12) number of 0-mismatch hits of the first 24bp
(13) number of 1-mismatch hits of the first 24 bp
(14) length of the read
(15) read sequence
(16) read sequence's quality

Column 6 provides the paired flag:
- First 4 bits denote orientation: 1=FF,2=FR,4=RF,8=RR
- Higher bits: 16=pair end requirement is met, 32=reads mapped to different
chromosomes, 64=one of the two reads cannot be mapped.
- Note that a correct pair flag is always 18.

 * @author aritz
 *
 */
public class MAPCursor {
	// scanner to read file.
	private Scanner scanner,linescanner;
	public ArrayList<Integer> flags;
	private int counter;
		
	// variables to store each time cursor advances:
	public String line;
	public String readname;
	public int chromosome;
	public int startposition;
	public int endposition;
	public char strand;
	public int flag;
	public int readlen;
	public int pairLen;
	public int mappingquality;
	public int readnum;
	
	// Constructor.  
	public MAPCursor(String filename) throws FileNotFoundException {
		scanner = new Scanner(new File(filename));
		flags = new ArrayList<Integer>();
		nullVariables();
		counter = 0;
	}
	
	// Constructor.  
	public MAPCursor(File file) throws FileNotFoundException {
		scanner = new Scanner(file);
		flags = new ArrayList<Integer>();
		nullVariables();
		counter = 0;
	}
	
	// sets all the relevant variables to null.
	private void nullVariables() {
		readname = null;
		chromosome = -1;
		startposition = -1;
		endposition = -1;
		strand = '0';
		flag = -1;
		readlen = -1;
		pairLen = -1;
	}
	
	public boolean hasNext() {
		if (scanner.hasNext())
			return true;
		return false;
			
	}
	
	public void close() {
		nullVariables();
		scanner.close();
	}
	
	/*
	 *  "advances" the cursor by loading the next line.  
	 *  Returns false if scanner doesn't have another line
	 */
	public boolean advance() throws FileNotFoundException {
		if (!hasNext()) {
			// set variables to null
			nullVariables();
			return false;
		}
		counter++; 
		line = scanner.nextLine();
		linescanner = new Scanner(line);
		try{
		 
		// replace variables with next line
		readname = linescanner.next();
		
		String chr = linescanner.next();
		//handle case where chr is something like chr17
		if (chr.length() > 2) {
			chr = chr.substring(3, chr.length());
		}
		if (chr.equalsIgnoreCase("X"))
			chromosome = 23;
		else if (chr.equalsIgnoreCase("Y"))
			chromosome = 24;
		else if (chr.equalsIgnoreCase("MT"))
			chromosome = 25;
		else
		chromosome = new Integer(chr).intValue();
		
		startposition = linescanner.nextInt();
		strand = linescanner.next().charAt(0);
		if (strand == 'p' || strand == 'P') {
			strand = '+';
		} else if (strand == 'm' || strand == 'M') {
			strand = '-';
		}
		pairLen = Math.abs(linescanner.nextInt()); 
		flag = linescanner.nextInt();

		linescanner.next(); // mapping quality
		linescanner.next(); // single-end mapping quality
		mappingquality = linescanner.nextInt(); // alternative mapping quality
		linescanner.next(); // # of mismatches of the best hit
		linescanner.next(); // sum of the qualities of mismatched bases of best hit
		linescanner.next(); // # of 0-mismatch hits in first 24bp
		linescanner.next(); // # of 1-mismatch hits in first 24bp
		readlen = linescanner.nextInt(); // read length
		endposition = startposition+readlen;
				
		// remove the last two characters of readname.
		readnum = new Integer(""+readname.charAt(readname.length()-1)).intValue();
		readname = readname.substring(0,readname.length()-2);

		// if flag hasn't been seen before, add it.
		if (!flags.contains(flag))
			flags.add(flag);

		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in line "+ counter + " before '"+linescanner.nextLine()+"'.  Skipping.");
			linescanner.close();
			return false;
		}
		//printLine();
		linescanner.close();
		return true;
	}

	/*
	 * Checks to see if the flag indicates the the paired read is in 
	 * forward/reverse orientation and on the same chromosome
	 */
	public boolean convergent() {
		if (flag==18 ||  flag == 2)
			return true;
		else return false;
	}
	
	/*
	 * checks to see if the flag indicates that the paired read is in FF or RR
	 * orientation and on the same chromosome.
	 */
	public boolean inversion() {
		if (flag==1 || flag == 8)
			return true;
		else return false;
	}
	
	public boolean diffChrom() {
		if (flag == 32) 
			return true;
		else return false;
	}
	
	public boolean divergent() {
		if (flag == 4)
			return true;
		else return false;
	}
	
	/*
	 * Checks to see if the name passed in is paired with readname.
	 * As far as I can tell, the names can be split into two parts:
	 * fragmentname/1 and fragmentname/2, where fragmentname is the same.
	 * 
	 * Only keep fragmentname, and compare this.
	 */
	public boolean isPair(String name) {
		return readname.equals(name);
	}
	
	public String genOutput() {
		return chromosome+"\t"+startposition+"\t"+endposition+"\t"+strand;
	}
	
	public void printLine() {
		System.out.println("name="+readname+",chr=" + chromosome+",startpos=" +startposition+",endpos=" +endposition+",strand="+strand+",flag=" + flag + ",readlen=" + readlen + ",pairlen=" + pairLen + ",qual= " + mappingquality + ",num=" + readnum);
	}
	
	

}
