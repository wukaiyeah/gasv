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
import gasv.geom.*;
import java.io.IOException;
import java.util.ArrayList;


//public class ReadSinglePred extends ReadFile {
	
	/* ReadCGH reads in file of CGH breakpoint region in form of one pair of genomic locations.
	 *  input file should be of format:
	 *  name chr pos1 pos2
	 *  
	 * Rectangles are created from CGH breakpoint locations as follows:
	 *  Large rectangles are created using each given CGH genomic region and each entire chromosome.
	 *  Those CGH pair that was used to create a rectangle that intersected with ESP data, is then paired with every other CGH pair
	 *  	to create many small candidate rectangles.
	 *  These rectangles are then put through the program once more.
	 */
/*
	public ReadSinglePred(String file, ArrayList<BreakRegion>[][] breakRegions) throws IOException{
		super(file, breakRegions);
	}
	public void readBreakRegions() throws IOException{
		String nextLine = b.readLine();
		int count = 0;
		while(nextLine != null){
			String[] line = nextLine.split("\\s+");
			
			int chr = chrFormatToNumber(line[1]);
			if(clones[0][chr-1]==null)
				clones[0][chr-1] = new ArrayList<BreakRegion>();
			CGH c = new CGH(line[0],chr,Integer.parseInt(line[2]),Integer.parseInt(line[3]));
			clones[0][chr-1].add(c);
			
			count++;
			nextLine = b.readLine();
		}
		sortByFirstPos(clones);
	}
*/
	// Make large rectangles from each clone and each chromosome region.
	/*
	public void makeSingletons(){
		System.out.println("Creating large rectangles");
		CGH[] chromReg = new CGH[numChrom];
		for(int i = 0;i < numChrom; i++){
			if(clones[0][i]!=null){
				ArrayList<BreakRegion> chr = clones[0][i];
				CGH firstCGH = (CGH) chr.get(0);
				int first = firstCGH.getX();
				CGH lastCGH = (CGH) chr.get(chr.size()-1);
				int last = lastCGH.getY();
				chromReg[i] = new CGH(""+i,i,first,last);
			}
		}
		for(int i = 0; i < numChrom; i++){ //for each bin containing a list of clones
			if(clones[0][i]!=null){
					ArrayList<BreakRegion> cloneList = clones[0][i];
					for(int j = 0; j < numChrom; j++){ // for each chromosome region (chromReg)
						ArrayList<Cluster> rectList = new ArrayList<Cluster>();
						for(int k = 0; k < cloneList.size(); k++) { //for each Clone in list
							CGH c = (CGH) cloneList.get(k);
							PolyDefault t = new Rectangle(c.getX(),c.getY(),chromReg[i].getX(),chromReg[i].getY());
							ArrayList<String> name = new ArrayList<String>();
							name.add(c.getName());
							name.add("chr" + (i+1));
							rectList.add(new Cluster(name,t));
						}
						shapes[i][j] = rectList;
					}	
			}
		}
	}
	*/

//}
