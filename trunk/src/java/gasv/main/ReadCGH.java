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
package gasv.main;
import gasv.common.Out;
import gasv.common.Constants;
import java.io.IOException;
import java.util.ArrayList;
import java.io.RandomAccessFile;


public class ReadCGH extends ReadFile {
	
	/* ReadCGH reads in file of CGH breakpoint region in form of one pair of genomic locations.
	 *  input file should be of format:
	 *  chr pos1 pos2
	 *  
	 * Rectangles are created from CGH breakpoint locations as follows:
	 *  Large rectangles are created using each given CGH genomic region and each entire chromosome.
	 *  Those CGH pair that was used to create a rectangle that intersected with ESP data, is then paired with every other CGH pair
	 *  	to create many small candidate rectangles.
	 *  These rectangles are then put through the program once more.
	 */

	private String nextLine_;

	//private java.io.FileReader f_;
	//private java.io.BufferedReader b_;
	private boolean diffChrPairReached_;
	private int curLeftChr_;
	private int curRightChr_;
	private int endWindowPos_ = 0;
	private RandomAccessFile cghFile_;

	public ReadCGH(String file, ArrayList<BreakRegion> breakRegions) throws IOException{
		super(file, breakRegions);

		//open in read-only mode
		cghFile_ = new RandomAccessFile(file, "r");
		//f_ = new java.io.FileReader(file_);
		//b_ = new java.io.BufferedReader(f_);
		curLeftChr_ = 0;
		curRightChr_ = 0;
		diffChrPairReached_ = false; 
		nextLine_ = null;
	}

	//ID probe1start probe1end probe2start probe2end coord1 coord2 coord3 coord4 chr
	/*private CGH parsePairedCGHFromString(String nextLine, int targetChr, 
			boolean useWindowSize, boolean matchChromosomes) 
			throws NumberFormatException, Exception {
		String[] line = nextLine.split("\\s+");

		
		//if operation fails, will throw NumberFormatException and skip row
		int chr = chrFormatToNumber(line[9]);
		double coord1 = makeSigned(line[5], "+");
		double coord2 = makeSigned(line[6], "+");
		double coord3 = makeSigned(line[7], "+");
		double coord4 = makeSigned(line[8], "+");

		// if breakRegionsArray is null, meaning that only want matching chr #'s,
		// ignore unless the current line's chromosomes matches the target chromosomes
		if (matchChromosomes && (chr != targetChr)) {
			Out.print3("ReadCGH: chr"+chr+ " doesn't match target chr"+targetChr + " so returning null");
			diffChrPairReached_ = true;
			//Out.print2("ReadESP: different chromosome encountered!");
			return null;

		}


		if (useWindowSize && (coord1 > endWindowPos_)) {
					//+ " and endWindowPos_=" + endWindowPos_);
			return null;
		} else {
			return new PairedCGH(line[0], line[1], line[2], line[3], line[4],
					chr, coord1, coord2, coord3, coord4);
		}
	}
	*/

	private CGH parseCGHFromString(String nextLine, int targetChr, 
			boolean useWindowSize, boolean matchChromosomes) 
			throws NumberFormatException, Exception {
		String[] line = nextLine.split("\\s+");
		//if (line.length == Constants.NUM_COLS_IN_PAIR_CGH_FILE) {
			//return this.parsePairedCGHFromString(nextLine, targetChr, useWindowSize, matchChromosomes);
		//}

		if (line.length != Constants.NUM_COLS_IN_SING_CGH_FILE) {
			System.out.println("Found a row with incorrect"
					+ " number of items. Expect "
					+ Constants.NUM_COLS_IN_SING_CGH_FILE
					+ " in an ESP file, found: " 
					+ line.length + ", so skipping row.");
			throw new Exception();
		}

		//if operation fails, will throw NumberFormatException and skip row
		int chr = chrFormatToNumber(line[2]);

		// if breakRegionsArray is null, meaning that only want matching chr #'s,
		// ignore unless the current line's chromosomes matches the target chromosomes
		if (matchChromosomes && (chr != targetChr)) {
			Out.print3("ReadCGH: chr"+chr+ " doesn't match target chr"+targetChr + " so returning null");
			diffChrPairReached_ = true;
			//Out.print2("ReadESP: different chromosome encountered!");
			return null;

		}

		double x = makeSigned(line[3],"+");
		double y  = makeSigned(line[4],"+");


		if (useWindowSize && (x > endWindowPos_)) {
					//+ " and endWindowPos_=" + endWindowPos_);
			return null;
		} else {
			return new CGH(line[0], line[1], chr, x, y);
		}
	}	

	public boolean readNextBreakRegions(int leftChr, int rightChr, 
			ArrayList<BreakRegion> br, int endWindowPos, 
			boolean startNewChr) throws IOException {
		return false;
	}
	/*
	public boolean readNextBreakRegions(int leftChr, int rightChr, 
			ArrayList<BreakRegion> br, int endWindowPos, 
			boolean startNewChr) throws IOException {
		
		if (leftChr != rightChr) {
			Out.print("ERROR: readNextBreakRegions for CGH data can only be on one chr, but "
					+ " leftChr: " + leftChr + " and rightChr: " + rightChr);
			throw new IOException();	
		}

		//if starting on a new chr, then need to clear the following flag
		if (startNewChr) {
			diffChrPairReached_ = false;
		}
		endWindowPos_ = endWindowPos;
		// only need to do this the first time through. After that, the previously read
		// nextLine_ will need to be parsed again to see if it matches this next chr pair or it will
		// need to be added to this new window
		if (nextLine_ == null) {
			nextLine_ = b_.readLine();
		}
		//will break out of this loop when first ESP outside of the current window is encountered
		while(nextLine_ != null){

			//Skip comment lines (those beginning with "#")
			if (nextLine_.startsWith("#")) {
				nextLine_ = b_.readLine();
				continue;
			}

			CGH c = null;

			try {
				c = this.parseCGHFromString(nextLine_, leftChr, true, true);
			} catch (NumberFormatException ex) {
				Out.print("Couldn't parse chromosome identifier: " 
						+ " for line: " + nextLine_ 
						+ " so skipping current line");
				nextLine_ = b_.readLine();
				continue;
			} catch (Exception ex) {
				//encountered a row with incorrect number of columns
				//Already printed error message so just skip to next line
				nextLine_ = b_.readLine();
				continue;
			}

			// If c is null we either reached a read on different chromosomes or a read that is beyond
			// the current window, but we haven't reached the end of the file yet, so return false 
			if (c == null) {
				return false;
			}


			// if no 2D array provided, put it in the ArrayList 
			br.add(c);

			nextLine_ = b_.readLine();
		}

		// if we reach the end of the file, close streams and return true
		b_.close();
		f_.close();
		return true;
	}
	*/

	public boolean getDiffChrPairReached() {
		return diffChrPairReached_;
	}

	/**
	  * Returns a 2 element int array. 
	  * First element is -1 if no CGH data matching this chromosome is found.
	  * If data is found, the first element will be the starting position of the first line of matching data.
	  * and the second element will be the end position (a.k.a. start pos of first line of non-matching data)
	  * Even if no data is found, the second element will still be the end position...
	  * unless end of file reached without hitting a different chromosome, then 2nd element will be -1.
	  */
	public long[] readBreakRegions(int chr, ArrayList<BreakRegion> cgh, long startPos) throws IOException {

		long[] ret = new long[2];
		ret[0] = -1;
		ret[1] = -1;
		cghFile_.seek(startPos);
		long curPos = startPos;
		
		String nextLine = cghFile_.readLine();
		//will break out of this loop when first CGH with non-matching chromosome is encountered 
		while(nextLine != null){

			//Skip comment lines (those beginning with "#")
			if (nextLine.startsWith("#")) {
				//update position to beginning of next line
				curPos = cghFile_.getFilePointer();
				nextLine = cghFile_.readLine();
				continue;
			}
			//Out.print2("ReadCGH.java: looking at line: " + nextLine);

			CGH c = null;

			try {
				//false so that we read in all CGH in the chromosome, not limited by win size
				// and true so that we match only CGH w/same chromosome
				c = this.parseCGHFromString(nextLine, chr, false, true);
			} catch (NumberFormatException ex) {
				Out.print("Couldn't parse chromosome identifier: " 
						+ " for line: " + nextLine 
						+ " so skipping current line");
				//update position to beginning of next line
				curPos = cghFile_.getFilePointer();
				nextLine = cghFile_.readLine();
				continue;
			} catch (Exception ex) {
				//encountered a row with incorrect number of columns
				//Already printed error message so just skip to next line
				//update position to beginning of next line
				curPos = cghFile_.getFilePointer();
				nextLine = cghFile_.readLine();
				continue;
			}

			// If c is null we reached a read on different chromosome so return the 
			// value of curPos which is the beginning of the line we just read
			if (c == null) {
				ret[1] = curPos; 
				return ret;
			}

			//otherwise we got a valid piece of CGH data
			//will need to return the beginning of the first such line
			if (ret[0] == -1) {
				ret[0] = curPos;
			}

			cgh.add(c);

			//update position to beginning of next line
			curPos = cghFile_.getFilePointer();
			nextLine = cghFile_.readLine();
		}

		// if we reach the end of the file, return the ret values.
		// To indicate end of file reached, the second return value should be left
		// as -1.
		return ret;
	}

		/*
	public void readBreakRegions(int leftChr, int rightChr) throws IOException{
		String nextLine = b.readLine();
		int count = 0;
		while(nextLine != null){
			String[] line = nextLine.split("\\s+");
			
			int chr = chrFormatToNumber(line[0]);
			if(clones[0][chr-1]==null)
				clones[0][chr-1] = new ArrayList<BreakRegion>();
			CGH c = new CGH("cgh_"+ count,chr,Integer.parseInt(line[1]),Integer.parseInt(line[2]));
			clones[0][chr-1].add(c);
			
			count++;
			nextLine = b.readLine();
		}
		sortByFirstPos(clones);
	}
		*/

	//CURRENTLY UNSUPPORTED!
	public void readBreakRegions(ArrayList<BreakRegion>[][] breakRegionsArray) throws IOException {
		readBreakRegions(-1, -1, breakRegionsArray);
	}

	//CURRENTLY UNSUPPORTED!
	public void readBreakRegions(int targetLeftChr, int targetRightChr) throws IOException{
		readBreakRegions(targetLeftChr, targetRightChr, null);
	}
	//CURRENTLY UNSUPPORTED!
	private void readBreakRegions(int targetLeftChr, int targetRightChr,
			ArrayList<BreakRegion>[][] breakRegionsArray) throws IOException{
		
	}
	
}

