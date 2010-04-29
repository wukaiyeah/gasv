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
import java.io.IOException;
import java.util.ArrayList;

/**
 * Wrapper for MAQprocessor, formerly ESPfilter.  
 * @author aritz
 *
 */
public class MAQProcessorMain {

	private enum MAQ_MODE {
		GET_VALIDS,
		GET_INVALIDS
	}
	
	// defaults
	private MAQ_MODE mode = null;
	private boolean rmdup = false;
	private String prefix = "out";
	private String numlines = "all";
	private String logfile = "out.log";
	private String minInvLen = "10000";
	private String minDelLen = "5000";
	private String printDiffChrom = "0";

	
	public static void printUsage() {
		System.out.println("\n########################################");
		System.out.println("General Usage: java -jar maqprocessor.jar <mode> [<optional args>] <args>");
		System.out.println("\n########################################");
		System.out.println("\n#### If <mode> is --valid or -v: ####\n");
		System.out.println("This mode determines the valid fragment bounds Lmin and Lmax and prints information to a log file.");
		System.out.println("Usage: java -jar maqprocessor.jar --valid [<optional args>] <InputDir> <Pattern>");
		System.out.println("\nParameters:");
		System.out.println("<InputDir> \tDirectory where extracted maq files are");
		System.out.println("<Pattern> \tRegular Expression (or single file name) to process");
		System.out.println("\nOptional Arguments:");
		System.out.println("-n, --numlines [N,all] \tFind up to N valid fragments in each file or 'all' for all lines. Default is 'all'");
		System.out.println("-l. --logfile <Logfile> \tname of log file to output to.  Default is 'out.log'");
		System.out.println("\nOutput:");
		System.out.println("<Logfile> \tprogram information, Lmin and Lmax, and a histogram of convergent pairs.");
		System.out.println("\n########################################");

		System.out.println("\n#### If <mode> is --invalid or -i: ####\n");
		System.out.println("This mode takes one or two chromosomes and outputs the invalid fragments.\n");
		System.out.println("Single Chromosome Usage : java -jar maqprocessor.jar --invalid [<optional arguments>] <InputFile> <Lmin> <Lmax>");
		System.out.println("For a single chromosome A, reports all divergent fragments, convergent but incorrectly-sized fragments, inversions, and deletions.");
		System.out.println("\nParameters:");
		System.out.println("<InputFile> \tfile for ChrA produced by mapview");
		System.out.println("<Lmin> \tLmin determined by output file in --valid mode.");
		System.out.println("<Lmax> \tLmax determined by output file in --valid mode,");
		System.out.println("\nOptional Arguments:");
		System.out.println("-p, --prefix <Prefix> \toutput file prefix, which can include output directory info and chromosome name.  Default is 'out'");
		System.out.println("-r, --rmdup \t remove duplicate fragments.  See below for description and output files.");
		System.out.println("--mininvlen N \tminimum distance N between inversion breakpoints. Default is 10000");
		System.out.println("--mindellen N \tminimum distance N between deletion breakpoints. Default is 5000");
		System.out.println("-d, --diffchrom \tprint all reads flagged as spanning chromosomes to <InputFile>.diffchrom.");
		System.out.println("\nOutput:");
		System.out.println("<Prefix>_inversions.invalid \tPES file of inversions >= 10KB");
		System.out.println("<Prefix>_convergent.invalid \tPES file of deletions & insertions");
		System.out.println("<Prefix>_divergent.invalid \tPES file of -/+ oriented fragments");
		System.out.println("<Prefix>_deletions.invalid \tPES file of deletions >= 5KB\n");
		System.out.println("\nTwo-Chromosome Usage : java -jar maqprocessor.jar --invalid [<optional arguments>] <InputFileA> <InputFileB>");
		System.out.println("For chromsomes A and B, reports all translocations from A to B (A is assumed to be less than B).");
		System.out.println("\nParameters:");
		System.out.println("<InputFileA> \tfile for ChrA produced by mapview");
		System.out.println("<InputFileB> \tfile for ChrB produced by mapview");
		System.out.println("\nOptional Arguments:");
		System.out.println("-p, --prefix <Prefix> \toutput file prefix, which can include output directory info and chromosome name.  Default is 'out'");
		System.out.println("-r, --rmdup \t remove duplicate fragments.  See below for description and output files.");
		System.out.println("\nOutput:");
		System.out.println("<Prefix>_translocations.invalid \tPES file of translocations for ChrA and ChrB\n");
		System.out.println("\nIf --rmdup is specified, duplicate fragments are removed, where duplicate fragments are fragments mapped to identical coordinates in identical orientations.  "+
				"This is similar to Maq's rmdup command, but is extended for invalid pairs.  There are two additional types of output files.");
		System.out.println("<*.invalid.rmdup> \t '*.invalid' file with duplicate reads removed.");
		System.out.println("<*.invalid.dup> \t list of fragments that were kept and duplicate fragments that were removed.");
		System.out.println("\nLastly, <Prefix>.reads contain the subset of Maq reads that are considered invalid.");
			return;
	}
	
	public static void main(String[] args) throws IOException, CloneNotSupportedException, NullPointerException{
		
		// remove first element of args
		MAQProcessorMain maq = new MAQProcessorMain();
		String[] newargs = maq.parseOptions(args);
		String[] argsToPass;
		if (maq.mode == MAQ_MODE.GET_VALIDS && newargs.length == 2) {
			argsToPass = new String[4];
			argsToPass[0] = newargs[0];
			argsToPass[1] = newargs[1];
			argsToPass[2] = maq.numlines;
			argsToPass[3] = maq.logfile;
			DetermineValidPairBounds.main(argsToPass);
		} else if (maq.mode == MAQ_MODE.GET_INVALIDS && newargs.length  == 3) {
			argsToPass = new String[8];
			argsToPass[0] = "true";
			argsToPass[1] = newargs[0];
			argsToPass[2] = maq.prefix;
			argsToPass[3] = newargs[1];
			argsToPass[4] = newargs[2];
			argsToPass[5] = maq.minInvLen;
			argsToPass[6] = maq.minDelLen;
			argsToPass[7] = maq.printDiffChrom;
			ExtractInvalidPairs.main(argsToPass);
		}	else if (maq.mode == MAQ_MODE.GET_INVALIDS && newargs.length  == 2) {
			argsToPass = new String[4];
			argsToPass[0] = "false";
			argsToPass[1] = newargs[0];
			argsToPass[2] = newargs[1];
			argsToPass[3] = maq.prefix;
			ExtractInvalidPairs.main(argsToPass);
		} else if (maq.mode == null) {
			printUsage();
			return;
		} else {
			System.out.println("Incorrect number of arguments:");
			printArgs(newargs,maq);
			return;
		}
		
		if (maq.mode == MAQ_MODE.GET_INVALIDS && maq.rmdup) {
			String[] names;
			if (argsToPass[0].equals("true")) {
				names = new String[3];
				names[0] = maq.prefix+"_divergent.invalid";
				names[1] = maq.prefix+"_deletions.invalid";
				names[2] = maq.prefix+"_inversions.invalid";
			} else {
				names = new String[1];
				names[0] = maq.prefix+"_translocations.invalid";
			}
			argsToPass = new String[3];
			for(int i=0;i<names.length;i++) {
				System.out.println("\tremoving duplicates for file " + names[i]);
				argsToPass[0] = names[i];
				argsToPass[1] = names[i]+".rmdup";
				argsToPass[2] = names[i]+".dup";
				RmDup.main(argsToPass);
			}
		}
	}
	
	private static void printArgs(String[] newargs,MAQProcessorMain m) {
		System.out.println("MAQ_MODE mode = " + m.mode);
		System.out.println("boolean rmdup = "+m.rmdup);
		System.out.println("String prefix = "+m.prefix);
		System.out.println("String numlines = "+m.numlines);
		System.out.println("String logfile = "+m.logfile);
		System.out.println("String minInvLen = "+m.minInvLen);
		System.out.println("String minDelLen = "+m.minDelLen);
		System.out.println("String diffchrom = " + m.printDiffChrom);
		System.out.println("arguments to pass to program:");
		for(int i=0;i<newargs.length;i++) {
			System.out.println("arg " + i +":"+newargs[i]);
		}
		
	}

	// sets arguments from default, returns arguments to be passed to the 
	// various methods.
	private String[] parseOptions(String[] args) throws IllegalArgumentException {
		ArrayList<String> arguments = new ArrayList<String>();
		//skip the first since that is the <mode> option of either -f or -c
		try {
		for (int i=0; i<args.length; i++) {
			if (args[i].equalsIgnoreCase("--valid") || args[i].equalsIgnoreCase("-v")) {
				mode = MAQ_MODE.GET_VALIDS;
			} else if (args[i].equalsIgnoreCase("--invalid") || args[i].equalsIgnoreCase("-i")) {
				mode = MAQ_MODE.GET_INVALIDS;
			} else if (args[i].equalsIgnoreCase("--rmdup") || args[i].equalsIgnoreCase("-r")) {
				rmdup = true;
			} else if (args[i].equalsIgnoreCase("--prefix") || args[i].equalsIgnoreCase("-p")) {
				i++;
				if (args[i].startsWith("-")) { throw new IllegalArgumentException(args[i] + " is not a valid prefix name."); }
				prefix = args[i];
			} else if (args[i].equalsIgnoreCase("--logfile") || args[i].equalsIgnoreCase("-l")) {
				i++;
				if (args[i].startsWith("-")) { throw new IllegalArgumentException(args[i] + " is not a valid logfile name."); }
				logfile = args[i];
			}  else if (args[i].equalsIgnoreCase("--numlines") || args[i].equalsIgnoreCase("-n")) {
				i++;
				if (!args[i].equalsIgnoreCase("all")) {
					int tmp = new Integer(args[i]).intValue();
					numlines = args[i];
				}
			} else if (args[i].equalsIgnoreCase("--mininvlen")) {
				i++;
				int tmp = new Integer(args[i]).intValue();
				minInvLen = args[i];
			} else if (args[i].equalsIgnoreCase("--mindellen")) {
				i++;
				int tmp = new Integer(args[i]).intValue();
				minDelLen = args[i];
			} else if (args[i].equalsIgnoreCase("--diffchrom") || args[i].equalsIgnoreCase("-d"))  {
				printDiffChrom = "1";
			} else {
				if (args[i].startsWith("-")) {
					throw new IllegalArgumentException(args[i] + " is not an identified option.");
				}
				arguments.add(args[i]);
			}
		}
		} catch (Exception e) { throw new IllegalArgumentException(e.getMessage()); }
		
		String[] newargs = new String[arguments.size()];
		for(int i=0;i<arguments.size();i++) {
			newargs[i] = arguments.get(i);
		}	
		return newargs;
	}
}
