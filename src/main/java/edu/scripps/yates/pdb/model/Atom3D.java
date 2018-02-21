package edu.scripps.yates.pdb.model;

import javax.vecmath.Point3d;

import edu.scripps.yates.pdb.read.PDBUtil;

/**
 * According to the PDB format,
 * ({@link http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM}
 * Atom lines in PDB format are:<br>
 * COLUMNS DATA TYPE FIELD DEFINITION<br>
 * -------------------------------------------------------------------------------------<br>
 * 1 - 6 Record name "ATOM "<br>
 * 7 - 11 Integer serial Atom serial number.<br>
 * 13 - 16 Atom name Atom name.<br>
 * 17 Character altLoc Alternate location indicator.<br>
 * 18 - 20 Residue name resName Residue name.<br>
 * 22 Character chainID Chain identifier.<br>
 * 23 - 26 Integer resSeq Residue sequence number.<br>
 * 27 AChar iCode Code for insertion of residues.<br>
 * 31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.<br>
 * 39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.<br>
 * 47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.<br>
 * 55 - 60 Real(6.2) occupancy Occupancy.<br>
 * 61 - 66 Real(6.2) tempFactor Temperature factor.<br>
 * 77 - 78 LString(2) element Element symbol, right-justified.<br>
 * 79 - 80 LString(2) charge Charge on the atom.<br>
 * 
 * @author Salva
 *
 */
public class Atom3D {

	protected final int atomNumber;
	protected final AtomType atomType;
	protected final String aa;
	protected int positionInPDB;
	protected final String chainID;
	private final static String sep = "\t";
	private final Point3d coordinates;

	public Atom3D(int atomNumber, AtomType atomType, String aa, int position, String chainID, double x, double y,
			double z) {
		this(atomNumber, atomType, aa, position, chainID, new Point3d(x, y, z));
	}

	public Atom3D(int atomNumber, AtomType atomType, String aa, int position, String chainID, Point3d coordinates) {
		super();
		this.atomNumber = atomNumber;
		this.atomType = atomType;
		this.aa = aa;
		positionInPDB = position;
		this.chainID = chainID;
		this.coordinates = coordinates;
	}

	/**
	 * Constructor of {@link Atom3D} from a PDB formatted line
	 *
	 * @param atomLine
	 * @throws IllegalArgumentException
	 */
	public Atom3D(String atomLine, boolean parseCoordinates) throws IllegalArgumentException {

		boolean useNewParsing = true;
		if (useNewParsing) {
			atomNumber = Integer.valueOf(atomLine.substring(6, 11).trim());
			atomType = AtomType.getByName(atomLine.substring(12, 16).trim());
			String aaString = atomLine.substring(17, 20).trim();
			if (aaString.length() != 3) {
				throw new IllegalArgumentException(
						"Error reading atom for protein. This atom may not belong to a protein");
			}
			aa = PDBUtil.parseAA(aaString);
			chainID = atomLine.substring(21, 22);
			positionInPDB = Integer.valueOf(atomLine.substring(22, 26).trim());
			if (parseCoordinates) {
				double x = Double.valueOf(atomLine.substring(30, 38).trim());
				double y = Double.valueOf(atomLine.substring(38, 46).trim());
				double z = Double.valueOf(atomLine.substring(46, 54).trim());
				this.coordinates = new Point3d(x, y, z);
			} else {
				this.coordinates = null;
			}
		} else {
			final String[] split = atomLine.split("\\s+");

			atomNumber = Integer.valueOf(split[1]);
			if (split[2].length() > 4) {
				atomType = AtomType.getByName(split[2].substring(0, split[2].length() - 3));
				aa = PDBUtil.parseAA(split[2].substring(split[2].length() - 3));
				chainID = split[3];
				positionInPDB = Integer.valueOf(split[4]);
				if (parseCoordinates) {
					double x = Double.valueOf(split[5]);
					double y = Double.valueOf(split[6]);
					double z = Double.valueOf(split[7]);
					this.coordinates = new Point3d(x, y, z);
				} else {
					this.coordinates = null;
				}
			} else {
				atomType = AtomType.getByName(split[2]);
				if (split[3].length() == 1) {
					aa = split[3];
				} else {
					aa = PDBUtil.parseAA(split[3].substring(split[3].length() - 3));
				}
				chainID = split[4];
				positionInPDB = Integer.valueOf(split[5]);
				if (parseCoordinates) {
					double x = Double.valueOf(split[6]);
					double y = Double.valueOf(split[7]);
					double z = Double.valueOf(split[8]);
					this.coordinates = new Point3d(x, y, z);
				} else {
					this.coordinates = null;
				}
			}
		}
	}

	public double distance(Atom3D atom) {
		return getCoordinates().distance(atom.getCoordinates());
	}

	public Point3d getCoordinates() {
		return coordinates;
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
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return toString(sep);
	}

	public String toString(String separator) {
		StringBuilder sb = new StringBuilder();
		sb.append(atomNumber).append(separator).append(atomType).append(separator).append(chainID).append(separator)
				.append(positionInPDB).append(separator).append(aa).append(separator)
				.append(getCoordinatesString(separator));
		return sb.toString();
	}

	private String getCoordinatesString(String separator) {
		StringBuilder sb = new StringBuilder();
		if (coordinates != null) {
			sb.append(coordinates.x + separator + coordinates.y + separator + coordinates.z);

		} else {
			sb.append("-" + separator + "-" + separator + "-");
		}
		return sb.toString();
	}

	public static String getToStringHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("atom_number").append(sep).append("atom_type").append(sep).append("chain_id").append(sep)
				.append("position_in_PDB").append(sep).append("AA").append(sep).append("x").append(sep).append("y")
				.append(sep).append("z");
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
	public static Atom3D getFromString(String string) {
		try {
			final String[] split = string.split(sep);
			int atomNumber = Integer.valueOf(split[0]);
			AtomType atomType = AtomType.getByName(split[1]);
			String chainID = split[2];
			int positionInPDB = Integer.valueOf(split[3]);
			String aa = split[4];
			Point3d coordinates = null;
			if (split.length == 6) {
				coordinates = new Point3d(Double.valueOf(split[5]), Double.valueOf(split[6]), Double.valueOf(split[7]));
			}
			Atom3D atom = new Atom3D(atomNumber, atomType, aa, positionInPDB, chainID, coordinates);
			return atom;
		} catch (Exception e) {
			return null;
		}
	}

	public static int numElementsInPrint() {
		return getToStringHeaders().split(sep).length;
	}

}
