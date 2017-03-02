package edu.scripps.yates.pdb.model;

import edu.scripps.yates.pdb.read.PDBUtil;

public class Atom {

	private final int atomNumber;
	private final AtomType atomType;
	private final String aa;
	private int positionInPDB;
	private final String chainID;
	private final static String sep = "\t";

	public Atom(int atomNumber, AtomType atomType, String aa, int position, String chainID) {
		super();
		this.atomNumber = atomNumber;
		this.atomType = atomType;
		this.aa = aa;
		positionInPDB = position;
		this.chainID = chainID;
	}

	/**
	 * Constructor of {@link Atom} from a PDB formatted line
	 *
	 * @param atomLine
	 * @throws IllegalArgumentException
	 */
	public Atom(String atomLine) throws IllegalArgumentException {
		if (atomLine.contains("1964  CA  GLY C1461")) {
			System.out.println(atomLine);
		}
		final String[] split = atomLine.split("\\s+");
		if (atomLine.contains("184A")) {
			System.out.println(atomLine);
		}
		atomNumber = Integer.valueOf(split[1]);
		if (split[2].length() > 4) {
			atomType = AtomType.getByName(split[2].substring(0, split[2].length() - 3));
			aa = PDBUtil.parseAA(split[2].substring(split[2].length() - 3));
			chainID = split[3];
			positionInPDB = Integer.valueOf(split[4]);
		} else {
			atomType = AtomType.getByName(split[2]);
			if (split[3].length() == 1) {
				aa = split[3];
			} else {
				aa = PDBUtil.parseAA(split[3].substring(split[3].length() - 3));
			}
			chainID = split[4];
			positionInPDB = Integer.valueOf(split[5]);
		}
	}

	/**
	 * @return the atomNumber
	 */
	public int getAtomNumber() {
		return atomNumber;
	}

	/**
	 * @return the atomType
	 */
	public AtomType getAtomType() {
		return atomType;
	}

	/**
	 * @return the aa
	 */
	public String getAa() {
		return aa;
	}

	/**
	 * @return the position
	 */
	public int getPositionInPDB() {
		return positionInPDB;
	}

	/**
	 * @return the chainID
	 */
	public String getChainID() {
		return chainID;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(atomNumber).append(sep).append(atomType).append(sep).append(chainID);
		return sb.toString();
	}

	public static String getToStringHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("atom_number").append(sep).append("atom_type").append(sep).append("chain_id");
		return sb.toString();
	}

	public void setPositionInPDB(int position2) {
		positionInPDB = position2;

	}

	/**
	 * Parse an string formated as the function toString()
	 *
	 * @param string
	 * @return
	 */
	public static Atom getFromString(String string, String aa, int position) {
		try {
			final String[] split = string.split(sep);
			int atomNumber = Integer.valueOf(split[0]);
			AtomType atomType = AtomType.getByName(split[1]);
			String chainID = split[2];
			Atom atom = new Atom(atomNumber, atomType, aa, position, chainID);
			return atom;
		} catch (Exception e) {
			return null;
		}
	}

}
