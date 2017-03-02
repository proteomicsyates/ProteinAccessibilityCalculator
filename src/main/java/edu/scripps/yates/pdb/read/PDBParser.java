package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.jmol.api.JmolViewer;

import edu.scripps.yates.pdb.model.Atom;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Chain;
import edu.scripps.yates.pdb.model.DBRef;
import edu.scripps.yates.pdb.surface.JMolScript;
import edu.scripps.yates.pdb.surface.SiteSurfaceAccessibilityReport;
import edu.scripps.yates.pdb.util.JMolCommandsUtil;
import edu.scripps.yates.pdb.util.SurfaceAccebilityInputParameters;
import edu.scripps.yates.utilities.strings.StringUtils;

public class PDBParser {
	private final static Logger log = Logger.getLogger(PDBParser.class);
	private final static String COMPND = "COMPND";
	private final static String DBREF = "DBREF";
	private final static String SEQRES = "SEQRES";
	private final static String ATOM = "ATOM";
	private final static String CA_ATOMS_ONLY = "CA ATOMS ONLY";
	private final static String MDLTYP = "MDLTYP";
	private static final int MIN_PEP_LENGTH = 6;
	private final String filePath;
	private static JmolViewer viewer;
	private ArrayList<Atom> atomList;
	private static boolean initialized = false;
	private HashMap<Integer, List<Atom>> atomsByPosition;
	private boolean opened = false;
	private final String pdbID;
	private final Map<String, SiteSurfaceAccessibilityReport> reportsByPDBPositionAndChain = new HashMap<String, SiteSurfaceAccessibilityReport>();
	private String selectedChainID;

	public PDBParser(String filePath, String pdbID) throws IOException {
		this.filePath = filePath;
		this.pdbID = pdbID;

	}

	private void init(boolean forceOpen) {
		if (!initialized) {
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

	public Map<String, SiteSurfaceAccessibilityReport> getSiteSurfaceAccesibilityFromParameters(
			SurfaceAccebilityInputParameters inputParameters) throws IOException, NotValidPDBException {
		if (inputParameters.getUniprotACC() != null) {
			Map<String, SiteSurfaceAccessibilityReport> ret = new HashMap<String, SiteSurfaceAccessibilityReport>();
			final SiteSurfaceAccessibilityReport report = getSiteSurfaceAccesibilityMappedToUniprot(inputParameters,
					true);
			ret.put(report.getPositionInPDB() + report.getAtom().getChainID(), report);
			return ret;
		} else {
			return getSiteSurfaceAccesibilityForPDBModel(inputParameters);
		}
	}

	private Map<String, SiteSurfaceAccessibilityReport> getSiteSurfaceAccesibilityForPDBModel(
			SurfaceAccebilityInputParameters inputParameters) throws IOException, NotValidPDBException {

		return getSiteSurfaceAccesibilityForPDBModel(inputParameters, inputParameters.getChain() != null);

	}

	private Map<String, SiteSurfaceAccessibilityReport> getSiteSurfaceAccesibilityForPDBModel(
			SurfaceAccebilityInputParameters inputParameters, boolean removeOtherChains)
			throws IOException, NotValidPDBException {
		// get the proteinSequence in PDB
		Chain chain = inputParameters.getChain();
		List<DBRef> dbRefs = new ArrayList<DBRef>();
		if (chain != null) {
			DBRef dbRef = PDBUtil.getDBRef(this, chain.getIdentifier());
			dbRefs.add(dbRef);
		} else {
			dbRefs.addAll(getDBRefs());
		}
		for (DBRef dbRef : dbRefs) {
			if (dbRef != null && dbRef.getChainID() != null) {
				chain = new Chain(pdbID, dbRef.getChainID() + "=0-0", -1.0f);
				final String pdbProteinSeq = getSequence(dbRef);
				if ("".equals(pdbProteinSeq)) {
					log.warn("Protein sequence for chain " + dbRef.getChainID() + " was not obtained");
					continue;
				}
				final List<Integer> aaPositionsInPDB = StringUtils.allPositionsOf(pdbProteinSeq,
						inputParameters.getAa());
				for (Integer positionInPDB : aaPositionsInPDB) {

					final String key = positionInPDB + chain.getIdentifier();

					if (reportsByPDBPositionAndChain.containsKey(key)) {
						log.info("Surface  accesibility already calculated for site " + positionInPDB + " in chain "
								+ chain.getIdentifier() + " of PDB " + chain.getPdbID());
						continue;
					}
					// get the appropiate atom in the chain of PDB
					final AtomType atomType = inputParameters.getAtomType();
					final String aa = inputParameters.getAa();
					Atom atom = getAtom(chain.getIdentifier(), aa, atomType, positionInPDB);
					if (atom == null) {
						log.warn("Atom is not found in chainID: '" + chain.getIdentifier() + "' in aa: '" + aa
								+ "' atomType: '" + atomType + "' position in PDB: '" + positionInPDB + "'");
						continue;
					}
					Double accesibility = getSurfaceAccesibilityOfAtom(atom, removeOtherChains);
					if (accesibility != null) {
						SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(dbRef.getPdbID(),
								accesibility, aa, atom, atomType, positionInPDB, inputParameters.getUniprotACC(),
								inputParameters.getPositionInUniprotProtein(), chain.getResolution());
						reportsByPDBPositionAndChain.put(key, report);

					}
				}
			} else {
				throw new NotValidPDBException("Chain " + chain.getIdentifier() + " not found in PDBParser " + this);
			}
		}
		return reportsByPDBPositionAndChain;
	}

	private SiteSurfaceAccessibilityReport getSiteSurfaceAccesibilityMappedToUniprot(
			SurfaceAccebilityInputParameters inputParameters, boolean removeOtherChains)
			throws IOException, NotValidPDBException {
		// get the proteinSequence in PDB
		final Chain chain = inputParameters.getChain();
		DBRef dbRef = PDBUtil.getDBRef(this, chain.getIdentifier());
		if (dbRef != null && dbRef.getChainID() != null) {
			final String pdbProteinSeq = getSequence(dbRef);
			if ("".equals(pdbProteinSeq)) {
				return null;
			}
			// see in which position the peptide is in the protein
			final String uniprotPeptideSeq = inputParameters.getUniprotPeptideSeq();
			if (pdbProteinSeq.contains(uniprotPeptideSeq)) {

				final int peptideStartInPDB = pdbProteinSeq.lastIndexOf(uniprotPeptideSeq) + 1;
				final int positionInPDB = peptideStartInPDB + inputParameters.getPositionInPeptide();
				final String key = positionInPDB + chain.getIdentifier();
				if (chain.getPdbID().equals("2RAI") && positionInPDB == 126) {
					System.out.println("asdf");
				}
				if (reportsByPDBPositionAndChain.containsKey(key)) {
					log.info("Surface accesibility already calculated for site " + positionInPDB + " in chain "
							+ chain.getIdentifier() + " of PDB " + chain.getPdbID());
					return reportsByPDBPositionAndChain.get(key);
				}
				// get the appropiate atom in the chain of PDB
				final AtomType atomType = inputParameters.getAtomType();
				final String aa = inputParameters.getAa();
				Atom atom = getAtom(chain.getIdentifier(), aa, atomType, positionInPDB);
				if (atom == null) {
					throw new NotValidPDBException(
							"Atom is not found with chainID: '" + chain.getIdentifier() + "' in aa: '" + aa
									+ "' atomType: '" + atomType + "' position in PDB: '" + positionInPDB + "'");
				}
				Double accesibility = getSurfaceAccesibilityOfAtom(atom, removeOtherChains);
				if (accesibility != null) {
					SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(dbRef.getPdbID(),
							accesibility, aa, atom, atomType, positionInPDB, inputParameters.getUniprotACC(),
							inputParameters.getPositionInUniprotProtein(), chain.getResolution());
					reportsByPDBPositionAndChain.put(key, report);
					return report;

				}
			} else {
				throw new NotValidPDBException("Peptide " + uniprotPeptideSeq + " is not found in PDB entry "
						+ dbRef.getPdbID() + " - " + dbRef.getChainID() + " with protein sequence: " + pdbProteinSeq);
			}
		} else {
			throw new NotValidPDBException("Chain " + chain.getIdentifier() + " not found in PDBParser " + this);
		}
		return null;
	}

	private Atom getAtom(String chainID, String aa, AtomType atomType, int positionInPDB)
			throws IOException, NotValidPDBException {

		final List<Atom> atomsByAAPosition = getAtomsByAAPosition(positionInPDB);
		if (atomsByAAPosition == null) {
			return null;
		}
		for (Atom atom : atomsByAAPosition) {
			if (atom.getPositionInPDB() == positionInPDB && atom.getChainID().equals(chainID) && atom.getAa().equals(aa)
					&& atom.getAtomType() == atomType) {
				return atom;
			}

		}
		return null;
	}

	private Double getSurfaceAccesibilityOfAtom(Atom atom, boolean removeOtherChains) {

		if (!atom.getChainID().equals(selectedChainID)) {
			init(true);
			executeCommands(JMolCommandsUtil.getSelectChainJMolScriptByAtom(atom, removeOtherChains));
			selectedChainID = atom.getChainID();
		} else {
			init(false);
		}
		// executeCommands(JMolCommandsUtil.getSelectAtomScript(atom));
		String strOutput = executeCommands(JMolCommandsUtil.getCalculateSurfaceScript(atom));

		return parseAccessibilityOutput(strOutput);
	}

	public String executeCommands(JMolScript commandSet) {
		init(false);

		log.info("Executing JMol command set:");
		log.info(commandSet.getCommandsToExecuteIndifferentLines());
		// getViewer().evalString(commandSet.getCommandsToExecute());
		String strOutput = (String) getViewer().scriptWaitStatus(commandSet.getCommandsToExecute(), null);

		log.info("Output of the command in JMol is " + strOutput);
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

	private List<Atom> getAtoms() throws IOException, NotValidPDBException {
		if (atomList == null) {
			atomList = new ArrayList<Atom>();
			if (!getLinesContaining(CA_ATOMS_ONLY).isEmpty()) {
				throw new NotValidPDBException(CA_ATOMS_ONLY + " entry");
			}
			final List<String> lines = getLinesStarting(ATOM);
			for (String string : lines) {
				try {
					Atom atom = new Atom(string);
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
		}
		return atomList;
	}

	private List<Atom> getAtomsByAAPosition(int aaPosition) throws IOException, NotValidPDBException {
		if (atomsByPosition == null) {
			atomsByPosition = new HashMap<Integer, List<Atom>>();
			final List<Atom> atoms = getAtoms();
			for (Atom atom : atoms) {
				final int position = atom.getPositionInPDB();
				if (atomsByPosition.containsKey(position)) {
					atomsByPosition.get(position).add(atom);
				} else {
					List<Atom> list = new ArrayList<Atom>();
					list.add(atom);
					atomsByPosition.put(position, list);
				}
			}
		}
		return atomsByPosition.get(aaPosition);
	}

	public List<Atom> getAtoms(String chainID, String aa, AtomType atomType) throws IOException, NotValidPDBException {
		List<Atom> ret = new ArrayList<Atom>();
		for (Atom atom : getAtoms()) {
			if (atom.getChainID().equals(chainID) && atom.getAa().equals(aa) && atom.getAtomType() == atomType) {
				ret.add(atom);
			}
		}
		return ret;

	}

	public List<DBRef> getDBRefs() throws IOException {
		final List<String> lines = getLinesContaining(DBREF);
		List<DBRef> list = new ArrayList<DBRef>();
		for (String dbLine : lines) {
			list.add(new DBRef(dbLine));
		}
		return list;
	}

	public String getSequence(DBRef dbRef) throws IOException, NotValidPDBException {
		StringBuilder sb = new StringBuilder();
		int position = -1;
		for (Atom atom : getAtoms()) {
			if (dbRef.getChainID().equals(atom.getChainID())) {
				if ("".equals(sb.toString())) {
					// if this is the first AA and the position is not 1, fill
					// the sequence with "X"s until that position
					if (atom.getPositionInPDB() > 1) {
						for (int i = 1; i < atom.getPositionInPDB(); i++) {
							sb.append("X");
						}
					}
				}
				if (position != atom.getPositionInPDB()) {
					sb.append(atom.getAa());
					position = atom.getPositionInPDB();
				}
			}
		}
		return sb.toString();
	}

	private List<String> getLinesContaining(String containing) throws IOException {

		// When filteredLines is closed, it closes underlying stream as well as
		// underlying file.
		return Files.lines(Paths.get(filePath)).filter(s -> s.contains(containing)).collect(Collectors.toList());
	}

	private List<String> getLinesStarting(String starting) throws IOException {
		Path path = Paths.get(filePath);
		// When filteredLines is closed, it closes underlying stream as well as
		// underlying file.
		return Files.lines(path).filter(s -> s.startsWith(starting)).collect(Collectors.toList());
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
		try {
			List<DBRef> dbRefs = getDBRefs();
			for (DBRef dbRef : dbRefs) {
				if (dbRef.getUniprotID() != null && dbRef.getUniprotID().equals(uniprotAcc)) {
					return true;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return false;
	}
}
