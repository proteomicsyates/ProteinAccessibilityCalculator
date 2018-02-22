package edu.scripps.yates.pdb.model;

public class Atom extends Atom3D {

	public Atom(int atomNumber, AtomType atomType, String aa, int position, String chainID) {
		super(atomNumber, atomType, aa, position, chainID, null);
	}

	public Atom(String atomLine) throws IllegalArgumentException {
		super(atomLine, false);
	}
}
