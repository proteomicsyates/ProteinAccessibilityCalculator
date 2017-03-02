package edu.scripps.yates.pdb.model;

public enum AtomType {
	N, CA, C, O, CB, CG, CG1, CG2, CD, CE, CE1, CE2, CE3, //
	CD1, CD2, CZ, CZ2, CZ3, NE, NE1, NE2, NZ, ND1, ND2, CH2, //
	NH1, NH2, OD, OD1, OD2, OE, OE1, OE2, OG, OG1, OH, SG, SD, //
	OXT, H, HA, HB1, HB2, HB3, HZ, HB, HG2, HG3, HD2, HD3, HE, //
	HH11, HH12, HH21, HH22, HG12, HG13, HG21, HG22, HG23, HD11, //
	HD12, HD13, HD21, HD22, HD23, HG11, HA2, HA3, HG, HE1, HE2, //
	HE3, HD1, HZ1, HZ2, HZ3, HH2, HE21, HE22, HG1, HH, NH1A, //
	NH2A, HB2A, OE1A, NE2A, HB3A, HG2A, HG3A, HE21A, HE22A, OE1B, //
	NE2B, HB2B, HB3B, OE2B, OE2A, HZ3B, HZ2B, HZ1B, HE2B, HE3B, //
	HD2B, HD3B, HG2B, HG3B, HD2A, HD3A, HE2A, HZ1A, HZ2A, HZ3A, C4, C3, C2, N1, C6, C5, C8, C1, O2, O3, O4, OP2, OP1, P, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, O5, O6, OP3;

	public static AtomType getByName(String string) {
		if (string.endsWith("'")) {
			return AtomType.valueOf(string.substring(0, string.length() - 1));
		} else {
			return AtomType.valueOf(string);
		}
	}
}
