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
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;


/**
 *  * Takes a MAP file, generates a distribution of the mapped clone lengths, and
 *  outputs the .001% .999% quantiles as the Lmin and Lmax boundaries for valid pairs.
 *  This needs to be done before we can run the chromosomes, since Lmin and Lmax must
 *  be determined once for the entire experiment.
 *   
 * @author aritz
 * @date 1/27/09
 */
public class DetermineValidPairBounds {
private HashMap<Integer,Integer> lengths;
private int numReadsPerFile;
private double quantile;

public DetermineValidPairBounds(int n,double q) {
		lengths = new HashMap<Integer,Integer>();
		numReadsPerFile = n;
		quantile = q;
	}
	
/**
 * Main method.
 * @param args: 
 * (1) the directory to look in
 * (2) regular expression for files to be used (or a single file)
 * (3) either 
 *   (a) the number of correctly-oriented pairs to take for each file OR
 *   (b) 'all' to use all correctly oriented pairs.
 * (4) the output log filename
 * @throws IOException 
 */
	public static void main(String[] args) throws IOException {
		
		  File thisdir = new File (".");
		  System.out.println(thisdir.getCanonicalPath()+File.separator+args[0]);
			 
		  File[] files;
		  if (new File(thisdir.getCanonicalPath()+File.separator+args[0]).exists()) {
			  files = new File(thisdir.getCanonicalPath()+File.separator+args[0]).listFiles(new DirFilter(args[1]));
		  }
		  else
			  files = new File(args[0]).listFiles(new DirFilter(args[1]));
		  
		System.out.println("There are " + files.length + " files that match '" + 
				args[1] + "' in directory " + args[0]);
		
		System.out.println("Writing information to "+args[3]);
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(args[3])));
		
		// write command line arguments & date executed
		String date =  (new SimpleDateFormat("MM-dd-yyyy")).format(Calendar.getInstance().getTime());
		writer.write("# executed on " + date + "\n");
		writer.write("# looking in directory " +args[0]+"\n");
		writer.write("# considering files satisfying regex '"+args[1]+"'\n");
		writer.write("# taking "+args[2]+" correctly-oriented pairs in each file\n");
		
		// for now, quantile values are hard-coded in.
		double q = .001;
		writer.write("# setting Lmin and Lmax to be "+q+" and "+ (1-q)+", respectively\n");
		
		int n = -1;
		if (args[2].equalsIgnoreCase("all")) 
			n = Integer.MAX_VALUE;
		else
			n = new Integer(args[2]).intValue();
		
		DetermineValidPairBounds dvpb = new DetermineValidPairBounds(n,q);
		
		// for each file, take top n correctly oriented pairs.
		for (int i=0;i<files.length;i++) {
			dvpb.addValidFragments(files[i]);
		}
		
		// take the .0001% and .999% quantiles.
		int[] lbounds = dvpb.getBounds();
		writer.write("\n# Lmin and Lmax\n");
		writer.write("Lmin="+lbounds[0]+"\n");
		writer.write("Lmax="+lbounds[1]+"\n");

		// print histogram of fragment lengths.
		writer.write("\n# Histogram of correctly oriented fragment lengths\n");
		writer.write("Fragment Length\tCount\n");
		Object[] keys = dvpb.lengths.keySet().toArray();
		Arrays.sort(keys);
		int len;
		for (int i=0;i<keys.length;i++) {
			len = (Integer)keys[i];
			writer.write(len+"\t"+dvpb.lengths.get(len)+"\n");
		}
		
		// close log file.
		writer.close();
	}


	public void addValidFragments(File file) throws IOException {
		System.out.println("reading file " + file);
		
		MAPCursor cursor = new MAPCursor(file);

		// assume this is done separately for each chromosome
		int counter = 0,validCounter=0;
		boolean a;
		while (cursor.hasNext()) {
			a = cursor.advance();
			if (!a) {
				break;
			}
		
			// only consider reads with high-scoring unique mapping locations 
			// (mapping quality cannot equal 0). Only consider reads 
			// with the correct paired orientation (FR) and ends are 
			// on the same chromosome.  If the strand is positive, then 
			// add it to the counts in the HashMap lengths.
			if (cursor.mappingquality >= 30 && 
					cursor.convergent() &&
					cursor.strand == '+') {
				validCounter++;
				if (lengths.get(cursor.pairLen) == null)
					lengths.put(cursor.pairLen,1);
				else
					lengths.put(cursor.pairLen,lengths.get(cursor.pairLen)+1);
			} 

			// if we have seen a certain number of valid reads, quit
			// TODO: determine whether this is sound.
			if (validCounter == numReadsPerFile)
				break;
		}
		cursor.close();
	
	}
	
	public int[] getBounds() {
		int[] lbounds = new int[2];		
		int numclones = 0;
		Iterator<Integer> vals = lengths.values().iterator();
		while (vals.hasNext()) 
			numclones+=vals.next();
		
		// expand hashmap to an integer
		Iterator<Integer> keys = lengths.keySet().iterator();
		int ind = 0;
		int len;
		int[] expanded = new int[numclones];
		while(keys.hasNext()) {
			len = (Integer)keys.next();
			for (int j=0;j<lengths.get(len);j++) {
				expanded[ind] = len;
				ind++;
			}
		}
		
		// get lmin 
		int lowQ = (int)Math.ceil(numclones*quantile);
		lbounds[0] = expanded[lowQ];
		
		// get lmax
		int highQ = (int)Math.floor(numclones*(1-quantile));
		lbounds[1] = expanded[highQ];
		
		System.out.println("Lmin = " + lbounds[0] + " Lmax = " + lbounds[1]);
		
		return lbounds;
	}
	
}
