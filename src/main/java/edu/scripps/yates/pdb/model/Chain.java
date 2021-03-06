package edu.scripps.yates.pdb.model;

import org.apache.log4j.Logger;

public class Chain {
	private final static Logger log = Logger.getLogger(Chain.class);
	private final String identifier;
	private int start;
	private final int end;
	private final String pdbID;
	private Float resolution;

	/**
	 *
	 * @param text
	 *            like 'A=328-363'
	 */
	public Chain(String pdbID, String text, Float resolution) {
		final String[] split = text.split("=");
		identifier = split[0].trim();
		if (split[1].contains("-") && !"-".equals(split[1])) {
			final String[] split2 = split[1].split("-");
			start = Integer.valueOf(split2[0].trim());
			end = Integer.valueOf(split2[1].trim());
		} else {
			start = -1;
			end = -1;
		}

		this.resolution = resolution;
		this.pdbID = pdbID.trim();
	}

	/**
	 * @return the identifier
	 */
	public String getIdentifier() {
		return identifier;
	}

	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}

	public boolean includesPosition(int position) {
		if (position >= start && position <= end) {
			return true;
		}
		return false;
	}

	/**
	 * @return the pdbID
	 */
	public String getPdbID() {
		return pdbID;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return pdbID + " [ chain:" + identifier + ", start=" + start + ", end=" + end + "]";
	}

	public Float getResolution() {
		return resolution;
	}

	public void setResolution(float resolution) {
		this.resolution = resolution;
	}
}
