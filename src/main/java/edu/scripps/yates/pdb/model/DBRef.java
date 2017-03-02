package edu.scripps.yates.pdb.model;

public class DBRef {
	private final String pdbID;
	private final String chainID;
	private final String uniprotID;

	public DBRef(String dbRefLine) {
		final String[] split = dbRefLine.split("\\s+");
		pdbID = split[1];
		chainID = split[2];
		if (split.length >= 7) {
			// some lines are like:
			// DBREF2 4KZX i 283837872 1 1863
			// in that case, we dont have uniprotID
			uniprotID = split[6];
		} else {
			uniprotID = null;
		}
	}

	/**
	 * @return the pdbID
	 */
	public String getPdbID() {
		return pdbID;
	}

	/**
	 * @return the chainID
	 */
	public String getChainID() {
		return chainID;
	}

	/**
	 * @return the uniprotID
	 */
	public String getUniprotID() {
		return uniprotID;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "DBRef [pdbID=" + pdbID + ", chainID=" + chainID + ", uniprotID=" + uniprotID + "]";
	}
}
