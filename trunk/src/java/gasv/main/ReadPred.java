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


//public class ReadPred extends ReadFile{
	
	/* ReadPred reads in file of breakpoint regions in form of two pairs of genomic locations, creating a rectangle in the plane.
	 *  input file should be of format:
	 *  clone_index chr1 brkpt1_pos1 brkpt1_pos2 chr2 brkpt2_pos1 brkpt2_pos2
	 */
/*
	public ReadPred(String file, ArrayList<BreakRegion>[][] breakRegions) throws IOException{
		super(file, breakRegions);
	}
	
	public void readBreakRegions() throws IOException{
		String nextLine = b.readLine();
		while(nextLine != null){
			String[] line = nextLine.split("\\s+");
			
			int leftChr = chrFormatToNumber(line[1]);
			int rightChr = chrFormatToNumber(line[4]);
			if(clones[leftChr-1][rightChr-1]==null)
				clones[leftChr-1][rightChr-1] = new ArrayList<BreakRegion>();
			PredictedRect c = new PredictedRect(line[0],leftChr,rightChr,Integer.parseInt(line[2]),Integer.parseInt(line[3]),
					Integer.parseInt(line[5]),Integer.parseInt(line[6]));
			clones[leftChr-1][rightChr-1].add(c);
			
			nextLine = b.readLine();
		}
		sortByFirstPos(clones);
	}
*/
	/*
	public void makeSingletons(){
		System.out.println("Creating rectangles");
		for(int i = 0; i < numChrom; i++){
			for(int j = 0; j < numChrom; j++){
				if(clones[i][j]!=null){
					ArrayList<BreakRegion> cloneList = clones[i][j];
					ArrayList<Cluster> rectList = new ArrayList<Cluster>();
					for(int k = 0; k < cloneList.size(); k++) {
						PredictedRect c = (PredictedRect) cloneList.get(k);
						PolyDefault t = new Rectangle(c.getX(),c.getX2(),c.getY(),c.getY2());
						ArrayList<String> name = new ArrayList<String>();
						name.add(c.getName());
						rectList.add(new Cluster(name,t));
					}
					shapes[i][j] = rectList;
				}
			}
		}
	}
	*/
//}
