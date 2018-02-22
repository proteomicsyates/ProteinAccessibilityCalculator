package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.jmol.api.JmolViewer;

import edu.scripps.yates.pdb.JMolScript;
import edu.scripps.yates.pdb.distance.model.Distance;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.DBRef;
import edu.scripps.yates.pdb.util.JMolCommandsUtil;
import gnu.trove.map.hash.TIntObjectHashMap;

public class PDBParser {
	private final static Logger log = Logger.getLogger(PDBParser.class);
	private final static String COMPND = "COMPND";
	private final static String DBREF = "DBREF";
	private final static String SEQRES = "SEQRES";
	private final static String EXPDTA = "EXPDTA";
	private final static String ATOM = "ATOM";
	private final static String CA_ATOMS_ONLY = "CA ATOMS ONLY";
	private final static String MDLTYP = "MDLTYP";
	private static final int MIN_PEP_LENGTH = 6;
	private static final String MUTATION = "MUTATION: YES";
	private final String filePath;
	private static JmolViewer viewer;
	private ArrayList<Atom3D> atomList;
	private static boolean initialized = false;
	private TIntObjectHashMap<List<Atom3D>> atomsByPosition;
	private boolean opened = false;
	private final String pdbID;
	private String selectedChainID;
	private String experimentalMethod;
	private Boolean mutation;
	private Map<String, List<String>> linesMap = new HashMap<String, List<String>>();
	private ArrayList<DBRef> dbRefs;
	private final boolean parseCoordinates;
	private final Map<String, Map<AtomType, List<Atom3D>>> atomsByAminoacidAndType = new HashMap<String, Map<AtomType, List<Atom3D>>>();
	private static double minDistanceBetweenDifferentChains = Double.MAX_VALUE;

	public PDBParser(String filePath, String pdbID, boolean parseCoordinates) throws IOException {
		this.filePath = filePath;
		this.pdbID = pdbID;
		this.parseCoordinates = parseCoordinates;
	}

	private void init(boolean forceOpen) {
		if (!initialized) {
			// set JMol log level to ERROR
			org.jmol.util.Logger.setLogLevel(org.jmol.util.Logger.LEVEL_ERROR);
			viewer = JmolViewer.allocateViewer(null, null);
			viewer.setScreenDimension(600, 600);

			initialized = true;
			opened = false;
		}
		if (!opened || forceOpen) {
			opened = true;
			executeCommands(new JMolScript("load " + filePath));
			// final String openFileError = viewer.openFile(filePath);

			// if (openFileError != null) {
			// log.error(openFileError);
			// throw new IllegalArgumentException(openFileError);
			// }
		}
	}

	private JmolViewer getViewer() {
		init(false);
		return viewer;
	}

	public Boolean getMutation() {
		if (mutation == null) {
			final List<String> lines = getLinesContaining(MUTATION);
			if (lines != null && !lines.isEmpty()) {
				mutation = true;
			} else {
				mutation = false;
			}
		}
		return mutation;
	}

	public Atom3D getAtom(String chainID, char aa, AtomType atomType, int positionInPDB) {

		final List<Atom3D> atomsByAAPosition = getAtomsByAAPosition(positionInPDB);
		if (atomsByAAPosition == null) {
			return null;
		}
		for (Atom3D atom : atomsByAAPosition) {
			if (atom.getPositionInPDB() == positionInPDB && atom.getChainID().equals(chainID)
					&& atom.getAa().equals(String.valueOf(aa)) && atom.getAtomType() == atomType) {
				return atom;
			}

		}
		return null;
	}

	public Double getSurfaceAccessibilityOfAtom(Atom3D atom, boolean removeOtherChains, boolean removeOtherMolecules) {

		// if (!atom.getChainID().equals(selectedChainID)) {
		init(true);
		executeCommands(JMolCommandsUtil.getSelectChainJMolScriptByAtom(atom, removeOtherChains, removeOtherMolecules));
		selectedChainID = atom.getChainID();
		// } else {
		// init(false);
		// }
		// executeCommands(JMolCommandsUtil.getSelectAtomScript(atom));
		String strOutput = executeCommands(JMolCommandsUtil.getCalculateSurfaceScript(atom));

		return parseAccessibilityOutput(strOutput);
	}

	public String executeCommands(JMolScript commandSet) {
		init(false);

		log.info("Executing JMol command set:" + commandSet.getCommandsToExecuteIndifferentLines());
		// getViewer().evalString(commandSet.getCommandsToExecute());
		String strOutput = (String) getViewer().scriptWaitStatus(commandSet.getCommandsToExecute(), null);
		if (!"".equals(strOutput)) {
			log.debug("Output of the command in JMol is " + strOutput);
		}
		return strOutput;
	}

	private Double parseAccessibilityOutput(String strOutput) {
		List<String> lines = new ArrayList<String>();
		if (strOutput.contains("\n")) {
			final String[] split = strOutput.split("\n");
			for (String string : split) {
				lines.add(string);
			}
		} else {
			lines.add(strOutput);
		}
		for (String line : lines) {
			if (line.startsWith("isosurfaceArea")) {
				line = line.replace(" ", "");
				if (line.contains("=")) {
					String s = line.split("=")[1];
					if (s.startsWith("[")) {
						s = s.substring(1);
					}
					if (s.endsWith("]")) {
						s = s.substring(0, s.length() - 1);
					}
					if (s.contains(",")) {
						final String[] numberArray = s.split(",");
						double total = 0.0;
						for (String string : numberArray) {
							total += Double.valueOf(string);
						}
						return total;
					} else {
						return Double.valueOf(s);
					}
				}
			}
		}
		return null;
	}

	private List<Atom3D> getAtoms() {
		if (atomList == null || atomList.isEmpty()) {
			atomList = new ArrayList<Atom3D>();

			// if (!getLinesContaining(CA_ATOMS_ONLY).isEmpty()) {
			// return atomList;
			// }
			Set<String> chainsReaded = new HashSet<String>();
			log.info("Reading PDB file...");
			List<String> lines = getLinesStarting(ATOM);
			log.info("Parsing " + lines.size() + " ATOM lines");
			for (String string : lines) {
				try {

					Atom3D atom = new Atom3D(string, parseCoordinates);
					if (!chainsReaded.contains(atom.getChainID())) {
						log.info(atom.getChainID() + " chain ");
						chainsReaded.add(atom.getChainID());
					}
					if (!atomsByAminoacidAndType.containsKey(atom.getAa())) {
						atomsByAminoacidAndType.put(atom.getAa(), new HashMap<AtomType, List<Atom3D>>());
					}
					if (!atomsByAminoacidAndType.get(atom.getAa()).containsKey(atom.getAtomType())) {
						atomsByAminoacidAndType.get(atom.getAa()).put(atom.getAtomType(), new ArrayList<Atom3D>());
					}
					atomsByAminoacidAndType.get(atom.getAa()).get(atom.getAtomType()).add(atom);
					// if (atom.getAtomType() == AtomType.N) {
					// position++;
					// }
					// atom.setPositionInPDB(position);
					atomList.add(atom);
				} catch (IllegalArgumentException e) {
					log.debug(e);
				} catch (Exception e) {
					log.debug(e);
				}
			}
			log.info(atomList.size() + " atom list acquired");
		}
		return atomList;
	}

	private List<Atom3D> getAtomsByAAPosition(int aaPosition) {
		if (atomsByPosition == null) {
			atomsByPosition = new TIntObjectHashMap<List<Atom3D>>();
			final List<Atom3D> atoms = getAtoms();
			for (Atom3D atom : atoms) {
				final int position = atom.getPositionInPDB();
				if (atomsByPosition.containsKey(position)) {
					atomsByPosition.get(position).add(atom);
				} else {
					List<Atom3D> list = new ArrayList<Atom3D>();
					list.add(atom);
					atomsByPosition.put(position, list);
				}
			}
		}
		return atomsByPosition.get(aaPosition);
	}

	public List<Atom3D> getAtoms(String chainID, String aa, AtomType atomType)
			throws IOException, NotValidPDBException {
		List<Atom3D> ret = new ArrayList<Atom3D>();
		for (Atom3D atom : getAtoms()) {
			if (atom.getChainID().equals(chainID) && atom.getAa().equals(aa) && atom.getAtomType() == atomType) {
				ret.add(atom);
			}
		}
		return ret;

	}

	public List<String> getChainIDs() {

		List<String> ret = new ArrayList<String>();
		List<DBRef> dbRefs = getDBRefs();
		for (DBRef dbRef : dbRefs) {
			String chainID = dbRef.getChainID();
			if (!ret.contains(chainID)) {
				ret.add(chainID);
			}
		}
		Collections.sort(ret);
		return ret;
	}

	public List<DBRef> getDBRefs() {
		if (dbRefs == null) {
			final List<String> lines = getLinesContaining(DBREF);
			dbRefs = new ArrayList<DBRef>();
			for (String dbLine : lines) {
				dbRefs.add(new DBRef(dbLine));
			}
		}
		return dbRefs;

	}

	public String getExperimentalMethod() {
		if (experimentalMethod == null) {

			List<String> lines = getLinesStarting(EXPDTA);
			if (!lines.isEmpty()) {
				String string = lines.get(0);
				experimentalMethod = string.substring(string.indexOf(EXPDTA) + EXPDTA.length()).trim();
			}

		}
		return experimentalMethod;
	}

	public String getSequence(DBRef dbRef) {
		StringBuilder sb = new StringBuilder();
		int position = -1;
		List<Atom3D> atoms = getAtoms();

		for (Atom3D atom : atoms) {
			if (dbRef.getChainID().equals(atom.getChainID())) {
				if ("".equals(sb.toString())) {
					// if this is the first AA and the position is not 1, fill
					// the sequence with "X"s until that position
					if (atom.getPositionInPDB() > 1) {
						for (int i = 1; i < atom.getPositionInPDB(); i++) {
							sb.append("?");
						}
					}
				}
				if (position != atom.getPositionInPDB()) {
					if (sb.toString().length() + 1 < atom.getPositionInPDB()) {
						for (int i = sb.toString().length() + 1; i < atom.getPositionInPDB(); i++) {
							sb.append("?");
						}
					}
					sb.append(atom.getAa());
					position = atom.getPositionInPDB();
				}
			}
		}
		if ("".equals(sb.toString())) {
			log.info(dbRef);
		}
		return sb.toString();
	}

	private List<String> getLinesContaining(String containing) {

		// When filteredLines is closed, it closes underlying stream as well as
		// underlying file.
		if (linesMap.containsKey(containing)) {
			log.info(containing + " already readed");
			return linesMap.get(containing);
		}
		try {
			List<String> list = Files.lines(Paths.get(filePath)).filter(s -> s.contains(containing))
					.collect(Collectors.toList());
			linesMap.put(containing, list);
			return list;
		} catch (IOException e) {
			return Collections.emptyList();
		}
	}

	private List<String> getLinesStarting(String starting) {
		Path path = Paths.get(filePath);
		// When filteredLines is closed, it closes underlying stream as well as
		// underlying file.
		try {
			return Files.lines(path).filter(s -> s.startsWith(starting)).collect(Collectors.toList());
		} catch (IOException e) {
			return Collections.emptyList();
		}
	}

	/**
	 * @return the pdbID
	 */
	public String getPdbID() {
		return pdbID;
	}

	/**
	 * @return the filePath
	 */
	public String getFilePath() {
		return filePath;
	}

	/**
	 * Gets the folder in which the PDB file is located
	 *
	 * @return
	 */
	public File getFileFolder() {
		File folder = new File(filePath).getParentFile();
		return folder;
	}

	public boolean containsUniprotReference(String uniprotAcc) {

		List<DBRef> dbRefs = getDBRefs();
		for (DBRef dbRef : dbRefs) {
			if (dbRef.getUniprotID() != null && dbRef.getUniprotID().equals(uniprotAcc)) {
				return true;
			}
		}

		return false;
	}

	// public Collection<Distance> getDistancesOfAtom(double distanceThreshold,
	// String aa1, AtomType atomType1, String aa2,
	// AtomType atomType2) {
	// Set<Distance> ret = new HashSet<Distance>();
	// List<String> chainIDs = getChainIDs();
	// for (int i = 0; i < chainIDs.size(); i++) {
	// String chain1ID = chainIDs.get(i);
	// for (int j = i; j < chainIDs.size(); j++) {
	// String chain2ID = chainIDs.get(j);
	// executeCommands(JMolCommandsUtil.getMeasureOffScript());
	// JMolScript calculateDistanceScript =
	// JMolCommandsUtil.getCalculateDistanceScript(chain1ID, aa1,
	// atomType1, chain2ID, aa2, atomType2);
	// executeCommands(calculateDistanceScript);
	// String output = executeCommands(JMolCommandsUtil.getMeasureListScript());
	// ret.addAll(parseDistanceOutput(distanceThreshold, output));
	// }
	// }
	// return ret;
	// }
	public Collection<Distance> getDistancesOfAtom(double distanceThreshold, Atom3D atom1, String aa2,
			AtomType atomType2) {
		Set<Distance> ret = new HashSet<Distance>();

		List<Atom3D> atoms2 = getAtomsByAminacidAndType(aa2, atomType2);

		for (int j = 0; j < atoms2.size(); j++) {
			Atom3D atom2 = atoms2.get(j);
			double distance = atom1.distance(atom2);
			if (!atom2.getChainID().contentEquals(atom1.getChainID())) {
				if (minDistanceBetweenDifferentChains > distance) {
					minDistanceBetweenDifferentChains = distance;
					log.info("distance of atoms of different chains = " + distance + "\t" + atom1.getPositionInPDB()
							+ atom1.getChainID() + "\t" + atom2.getPositionInPDB() + "-" + atom2.getAa() + "-"
							+ atom2.getChainID());

				}
			}

			if (distance <= distanceThreshold) {
				ret.add(new Distance(atom1, atom2, distance));
			}
		}

		return ret;
	}

	private List<Atom3D> getAtomsByAminacidAndType(String aa, AtomType atomType) {
		getAtoms();
		if (atomsByAminoacidAndType.containsKey(aa)) {
			Map<AtomType, List<Atom3D>> atomsByType = atomsByAminoacidAndType.get(aa);
			if (atomsByType.containsKey(atomType)) {
				return atomsByType.get(atomType);
			}
		}
		return Collections.emptyList();
	}

	// private Set<Distance> parseDistanceOutput(double distanceThreshold,
	// String output) {
	// Set<Distance> ret = new HashSet<Distance>();
	// List<String> lines = new ArrayList<String>();
	//
	// if (output.contains("\n")) {
	// boolean first = true;
	// String[] split = output.split("\n");
	// for (String line : split) {
	// if (!first) {
	// lines.add(line);
	// } else {
	// first = false;
	// }
	// }
	// } else {
	// throw new IllegalArgumentException("Check format of output: " + output);
	// }
	//
	// for (String line : lines) {
	// if (line.contains("\t")) {
	// String[] cols = line.split("\t");
	// try {
	// double distance = Double.valueOf(cols[1]);
	// if (distance <= distanceThreshold) {
	//
	// Atom3D originalAtom = parseAtomFromDistanceReport(cols[3]);
	// Atom3D endAtom = parseAtomFromDistanceReport(cols[4]);
	// ret.add(new Distance(originalAtom, endAtom, distance));
	// }
	// } catch (NumberFormatException e) {
	// e.printStackTrace();
	// }
	// } else {
	//
	// throw new IllegalArgumentException("Check format of a line of the output:
	// " + line);
	//
	// }
	// }
	// return ret;
	// }

	// /**
	// * Parses a string like this into a {@link Atom3D}: <br>
	// * <b>[LYS]263:A.NZ #9006</b>
	// *
	// * @param string
	// * @return
	// */
	// private Atom3D parseAtomFromDistanceReport(String string) {
	// try {
	// String[] spaceSplitted = string.split(" ");
	// int atomNumber = Integer.valueOf(spaceSplitted[1].substring(1));
	// AtomType atomType = AtomType.valueOf(spaceSplitted[0].split("\\.")[1]);
	// String threeLetterAA = string.substring(string.indexOf("[") + 1,
	// string.indexOf("]"));
	// String aa = JMolCommandsUtil.convertToOneLetterAA(threeLetterAA);
	// int positionInPDB = Integer.valueOf(string.substring(string.indexOf("]")
	// + 1, string.indexOf(":")));
	// String chainID = string.substring(string.indexOf(":") + 1,
	// string.indexOf("."));
	// Atom3D atom = new Atom3D(atomNumber, atomType, aa, positionInPDB,
	// chainID);
	// return atom;
	// } catch (Exception e) {
	// e.printStackTrace();
	// throw new IllegalArgumentException("Error parsing atom from output of
	// distances: " + string);
	// }
	//
	// }
}
