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
import java.util.HashMap;
import java.util.Scanner;


/**
 * Removes duplicates in a file.  Assumes the file is in the new input format.
 * @author aritz
 *
 */
public class RmDup {

	public static void main(String[] args) throws IOException {
		String inputfile = args[0];
		String outputfile = args[1];
		String dupFile = args[2];
		Scanner scan = new Scanner(new File(inputfile));
		BufferedWriter writer = new BufferedWriter(new FileWriter(outputfile));
		BufferedWriter dups = new BufferedWriter(new FileWriter(dupFile));
		
		HashMap<String,String> seen = new HashMap<String,String>();
		
		String name,chrA,startposA,endposA,chrB,startposB,endposB;
		char signA,signB;
		String line,origline;
		Scanner linescan;
		while (scan.hasNext()) {
			
			origline = scan.nextLine();
			try{
			linescan = new Scanner(origline);
			name = linescan.next(); // name
			chrA = linescan.next();
			startposA = linescan.next();
			endposA = linescan.next();
			signA = linescan.next().charAt(0);
			chrB = linescan.next();
			startposB = linescan.next();
			endposB = linescan.next();
			signB = linescan.next().charAt(0);
			if (signA == '+' && signB == '+') 
				line = chrA+startposA+signA+chrB+startposB+signB;
			else if (signA == '+' && signB == '-') 
				line = chrA+startposA+signA+chrB+endposB+signB;
			else if (signA == '-' && signB == '+') 
				line = chrA+endposA+signA+chrB+startposB+signB;
			else if (signA == '-' && signB == '-')
				line = chrA+endposA+signA+chrB+endposB+signB;
			else {
				System.out.println("Error: signA and signB are other than +/-");
				line = "";
			}
			
			if (seen.containsKey(line)) {
				dups.write(seen.get(line)+"\t"+name+"\n");
			} else {
				writer.write(origline+"\n");
				seen.put(line,name);
			}

			} catch (Exception e) {
				System.out.println("skipping line " + origline);
			}
		}
		writer.close();
		dups.close();
	}
}
