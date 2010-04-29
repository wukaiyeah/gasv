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
import java.io.FilenameFilter;
import java.util.regex.Pattern;

/**
 * Takes a regular expression (can be a single file name) and determines whether
 * a given file name in a given directory satisfies the expression.
 * @author aritz
 * @date 1/17/09
 */
public class DirFilter implements FilenameFilter {
	private String regex;

	public DirFilter(String r) {
		regex = r;
	}

	/*
	 * Method for interface FilenameFilter.
	 * @see java.io.FilenameFilter#accept(java.io.File, java.lang.String)
	 */
	public boolean accept(File dir, String fname) {
		if (Pattern.matches(regex,fname))
			return true;
		return false;
	}
}
