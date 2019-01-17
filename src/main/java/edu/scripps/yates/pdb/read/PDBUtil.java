package edu.scripps.yates.pdb.read;

import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.pdb.model.DBRef;
import edu.scripps.yates.utilities.annotations.uniprot.xml.DbReferenceType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.PropertyType;
import gnu.trove.set.hash.THashSet;

public class PDBUtil {
	private final static Logger log = Logger.getLogger(PDBUtil.class);

	public static String parseAA(String threeLetterAA) {
		if (threeLetterAA.length() != 3) {
			return " ";
		}
		switch (threeLetterAA) {
		case "ALA":
			return "A";
		case "HIS":
			return "H";
		case "ARG":
			return "R";
		case "PHE":
			return "F";
		case "CYS":
			return "C";
		case "GLY":
			return "G";
		case "GLN":
			return "Q";
		case "GLU":
			return "E";
		case "ASP":
			return "D";
		case "LYS":
			return "K";
		case "LEU":
			return "L";
		case "MET":
			return "M";
		case "ASN":
			return "N";
		case "SER":
			return "S";
		case "TYR":
			return "Y";
		case "THR":
			return "T";
		case "ILE":
			return "I";
		case "TRP":
			return "W";
		case "PRO":
			return "P";
		case "VAL":
			return "V";
		case "CGU":
			return "E";
		case "UNK":// unknown
			return "X";
		default:
			log.warn(threeLetterAA + " is not a recognizable AA");
			return " ";
		}
	}

	public static String getPropertyValueFromDbReferenceType(DbReferenceType dbReferenceType, String propertyName) {
		if (dbReferenceType != null) {
			for (PropertyType prop : dbReferenceType.getProperty()) {
				if (prop.getType().equals(propertyName)) {
					return prop.getValue();
				}
			}
		}
		return null;
	}

	public static DBRef getDBRef(PDBParser parser, String chainID) {
		if (parser != null) {
			Set<String> individualChainIDs = new THashSet<String>();
			if (chainID.contains("/")) {
				String[] split = chainID.split("/");
				for (int i = 0; i < split.length; i++) {
					individualChainIDs.add(split[i]);
				}
			} else {
				individualChainIDs.add(chainID);
			}
			final List<DBRef> dbRefs = parser.getDBRefs();
			for (DBRef dbRef : dbRefs) {
				for (String individualChainID : individualChainIDs) {
					if (dbRef.getChainID().equals(individualChainID)) {
						return dbRef;
					}
				}
			}
		}
		return null;
	}

}
