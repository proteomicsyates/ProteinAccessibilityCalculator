package edu.scripps.yates.pdb.util;

import java.io.File;
import java.util.Arrays;

import edu.scripps.yates.dbindex.util.PeptideFilterByMaxOccurrencies;
import edu.scripps.yates.utilities.fasta.dbindex.PeptideFilter;

public class FastaDigestionConfiguration {
	private final File fasta;
	private final char[] enzymeArray;
	private final int numMisscleavages;
	private final boolean semiCleavage;
	private final String peptideFilterString;
	private final boolean ignorePeptidesNotFoundInDB;

	public FastaDigestionConfiguration(File fasta, char[] enzymeArray, int numMisscleavages, boolean semiCleavage,
			String peptideFilterString, boolean ignorePeptidesNotFoundInDB) {
		super();
		this.fasta = fasta;
		this.enzymeArray = enzymeArray;
		this.numMisscleavages = numMisscleavages;
		this.semiCleavage = semiCleavage;
		this.peptideFilterString = peptideFilterString;
		this.ignorePeptidesNotFoundInDB = ignorePeptidesNotFoundInDB;
	}

	/**
	 * @return the fasta
	 */
	public File getFasta() {
		return fasta;
	}

	/**
	 * @return the enzymeArray
	 */
	public char[] getEnzymeArray() {
		return enzymeArray;
	}

	/**
	 * @return the numMisscleavages
	 */
	public int getNumMisscleavages() {
		return numMisscleavages;
	}

	/**
	 * @return the semiCleavage
	 */
	public boolean isSemiCleavage() {
		return semiCleavage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "[enzymeArray=" + Arrays.toString(enzymeArray) + ", numMisscleavages=" + numMisscleavages
				+ ", semiCleavage=" + semiCleavage + ", peptideFilterString=" + peptideFilterString
				+ ", ignorePeptidesNotFoundInDB=" + ignorePeptidesNotFoundInDB + "]";
	}

	public PeptideFilter getPeptideFilter() {
		String aa = String.valueOf(peptideFilterString.charAt(0));
		int numMax = Integer.valueOf(String.valueOf(peptideFilterString.charAt(1)));
		return new PeptideFilterByMaxOccurrencies(aa, numMax);
	}

	public boolean isIgnorePeptidesNotFoundInDB() {
		return ignorePeptidesNotFoundInDB;
	}

}
