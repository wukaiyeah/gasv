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
import java.awt.*;
import java.util.*;
import java.io.*;
import java.text.*;

public class ConvertFileToPolygons{

	private static final String ESP_HEAD_STD = "#Cluster_ID:\tLeftChr:\tLeftBreakPoint:\tRightChr:\tRightBreakPoint:\tNum PES:\tLocalization:\tType:\n";
	private static final String ESP_HEAD_READS = "#Cluster_ID:\tLeftChr:\tLeftBreakPoint:\tRightChr:\tRightBreakPoint:\tNum PES:\tLocalization:\tType:\tList of PES:\n";
	private static final String ESP_HEAD_REGIONS = "#Cluster_ID:\tNum PES:\tLocalization:\tType:\tList of PES:\t LeftChr:\tRightChr:\tBoundary Points:\n";
	public static void convertFileToPolys(String inputFile, int lmin, int lmax, boolean useBatch, boolean useFast) 
			throws java.io.IOException, NullPointerException{

		String inFileNameMinusPath = inputFile;
		int idxOfSlash = inputFile.lastIndexOf("/");
		if (idxOfSlash != -1) {
			inFileNameMinusPath = inputFile.substring(idxOfSlash + 1);
		}

		String fileName = GASVMain.OUTPUT_DIR + inFileNameMinusPath + ".noclusters"; 
		File polyFile = new File(fileName);
		if (polyFile.exists()) {
			polyFile.delete();
		}
		BufferedWriter polyWriter = new BufferedWriter(new FileWriter(polyFile));
		
		if (GASVMain.USE_HEADER) {
			//String headerLine = "#Clone Name: \tChrM:\tChrN:\tBoundary Points:\n";
			//String headerLine = "#Cluster_ID:\tNum PES:\tLocalization:\tList of PES:\t ChrM:\tChrN:\tBoundary Points:\n";
			String headerLine = ESP_HEAD_STD;
			if (GASVMain.OUT_MODE == GASVMain.GASV_OUTPUT_MODE.REGIONS) {
				headerLine = ESP_HEAD_REGIONS;
			} else if (GASVMain.OUT_MODE == GASVMain.GASV_OUTPUT_MODE.READS) {
				headerLine = ESP_HEAD_READS;
			} else if (GASVMain.OUT_MODE == GASVMain.GASV_OUTPUT_MODE.STANDARD) {
				headerLine = ESP_HEAD_STD;
			}
			polyWriter.write(headerLine);
		}
		

		ReadInput readInput = new ReadInput(inputFile);
		ArrayList<BreakRegion>[][] breakRegions = null; 
		boolean fileDone = false;

		//in fast mode, read PES from all chromosomes into memory at once
		if (useFast) {
			breakRegions = new ArrayList[GASVMain.NUM_CHROM][GASVMain.NUM_CHROM];
			if (useBatch) {
				readInput.readFiles(breakRegions);
			} else {
				readInput.readSingleFile(lmin, lmax, breakRegions);
			}
			fileDone = true;
		}

		//Printing at most one decimal point.
		DecimalFormat df = new DecimalFormat("#.#");

		int count = 0;
		String cid = "c";
		String tab = "\t";
		String num_pes = "1\t";
		boolean saveMemory = GASVMain.SAVE_MEMORY;
		for(int i = 0; i<GASVMain.NUM_CHROM; i++){
			for(int j = i; j<GASVMain.NUM_CHROM; j++){
				int chrx = i+1;
				int chry = j+1;
				Out.print1("ConvertFileToPolygons: processing chr " + chrx + ", chr" + chry);
				ArrayList<BreakRegion> c = null;
				do {
				if (useFast) {
					c = breakRegions[i][j];
				} else {
					if (c == null) {
						c = new ArrayList<BreakRegion>();
					}
					if (useBatch) {
						//c = readInput.readFiles(chrx, chry);
						fileDone = readInput.readWindowFromFiles(chrx, 
										chry, c);
					} else {
						//c = readInput.readSingleFile(chrx, chry, lmin, lmax);
						fileDone = readInput.readWindowFromSingleFile(chrx, 
										chry, lmin, lmax, c);
					}
				}
				while ((c != null) && c.size() > 0) {
					BreakRegion br = c.remove(0);
					if (br instanceof Clone) {
						Clone clone = (Clone) br; 
						++count;
						
						Cluster curCluster = new Cluster(clone.getChrX(), clone.getChrY());
						curCluster.addClone(clone);
						curCluster.setIsCliqueAndMakePoly(true);
						ClusterESP.writeSingleCluster(cid+count, curCluster, false, polyWriter);
						/*
						//ID
						polyWriter.write(cid + count + tab);
						// Num PES
						polyWriter.write(num_pes);
						//Localization
						polyWriter.write(df.format(Math.sqrt(clone.getPoly().getArea())) + tab);
						//List of PES
						polyWriter.write(clone.getName() + tab);
						//ChrM and ChrN
						polyWriter.write(chrx + tab + chry + tab);
						
						//Boundary Points
						polyWriter.write(ClusterESP.printPoly((PolyDefault)clone.getPoly()));
						polyWriter.write("\n");
						*/
					} else {
						Out.print("ConvertFileToPolygons: Error! Unrecognized "
								+ "BreakRegion type: " + br);
					}

				}//end while
				} while (saveMemory && !fileDone 
						&& !readInput.getDiffChrPairReached() );

				if (saveMemory && fileDone) {
					polyWriter.flush();
					polyWriter.close();
					return;
				}
				
			}//end for
		}//end for
		
		if (!fileDone) {
			Out.print("WARNING: Finished looping through all possible Chr combinations, but didn't make it"
				       + " through all input file(s)! The inputs file(s) are either not sorted "
				       + "properly or have additional data with chromosome number > " 
				       + GASVMain.NUM_CHROM + "!");
		}
	
		polyWriter.flush();
		polyWriter.close();
	}
}
