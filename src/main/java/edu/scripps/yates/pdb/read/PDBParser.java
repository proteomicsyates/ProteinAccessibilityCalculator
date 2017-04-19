package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
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
	private final static String EXPDTA = "EXPDTA";
	private final static String ATOM = "ATOM";
	private final static String CA_ATOMS_ONLY = "CA ATOMS ONLY";
	private final static String MDLTYP = "MDLTYP";
	private static final int MIN_PEP_LENGTH = 6;
	private static final String MUTATION = "MUTATION: YES";
	private final String filePath;
	private static JmolViewer viewer;
	private ArrayList<Atom> atomList;
	private static boolean initialized = false;
	private HashMap<Integer, List<Atom>> atomsByPosition;
	private boolean opened = false;
	private final String pdbID;
	private static final Map<String, SiteSurfaceAccessibilityReport> reportsByPDBPositionAndChainAndRemovals = new HashMap<String, SiteSurfaceAccessibilityReport>();
	private String selectedChainID;
	private String experimentalMethod;
	private Boolean mutation;

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
			SurfaceAccebilityInputParameters inputParameters) {
		if (inputParameters.getUniprotACC() != null) {
			Map<String, SiteSurfaceAccessibilityReport> ret = new HashMap<String, SiteSurfaceAccessibilityReport>();
			final SiteSurfaceAccessibilityReport report = getSiteSurfaceAccesibilityMappedToUniprot(inputParameters);
			if (report != null) {
				ret.put(report.getPositionInPDB() + report.getAtom().getChainID(), report);
			}
			return ret;
		} else {
			return getSiteSurfaceAccesibilityForPDBModel(inputParameters);
		}
	}

	private Map<String, SiteSurfaceAccessibilityReport> getSiteSurfaceAccesibilityForPDBModel(
			SurfaceAccebilityInputParameters inputParameters) {
		// get the proteinSequence in PDB
		Chain chain = inputParameters.getChain();
		List<DBRef> dbRefs = new ArrayList<DBRef>();
		if (chain != null) {
			DBRef dbRef = PDBUtil.getDBRef(this, chain.getIdentifier());
			dbRefs.add(dbRef);
		} else {
			dbRefs.addAll(getDBRefs());
		}
		boolean mutation = getMutation();
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

					final String key = positionInPDB + chain.getIdentifier() + inputParameters.isRemoveOtherChains()
							+ inputParameters.isRemoveOtherMolecules();

					if (reportsByPDBPositionAndChainAndRemovals.containsKey(key)) {
						log.info("Surface  accesibility already calculated for site " + positionInPDB + " in chain "
								+ chain.getIdentifier() + " of PDB " + chain.getPdbID() + " with removingOtherChains="
								+ inputParameters.isRemoveOtherChains() + " and removingOtherMolcules="
								+ inputParameters.isRemoveOtherMolecules());
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
					Double accesibility = getSurfaceAccesibilityOfAtom(atom, inputParameters.isRemoveOtherChains(),
							inputParameters.isRemoveOtherMolecules());
					if (accesibility != null) {
						SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(dbRef.getPdbID(),
								accesibility, aa, atom, atomType, positionInPDB, inputParameters.getUniprotACC(),
								inputParameters.getPositionInUniprotProtein(), chain.getResolution(),
								inputParameters.isRemoveOtherChains(), inputParameters.isRemoveOtherMolecules(),
								mutation, getExperimentalMethod());
						reportsByPDBPositionAndChainAndRemovals.put(key, report);

					}
				}
			} else {
				log.debug("Chain " + chain.getIdentifier() + " not found in PDBParser " + this);
			}
		}
		return reportsByPDBPositionAndChainAndRemovals;
	}

	private SiteSurfaceAccessibilityReport getSiteSurfaceAccesibilityMappedToUniprot(
			SurfaceAccebilityInputParameters inputParameters) {
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
				final String key = inputParameters.getReportKey();
				if (reportsByPDBPositionAndChainAndRemovals.containsKey(key)) {
					log.info("Surface accesibility already calculated for site " + positionInPDB + " in chain "
							+ chain.getIdentifier() + " of PDB " + chain.getPdbID() + " with removingOtherChains="
							+ inputParameters.isRemoveOtherChains() + " and removingOtherMolcules="
							+ inputParameters.isRemoveOtherMolecules());
					return reportsByPDBPositionAndChainAndRemovals.get(key);
				}
				// get the appropiate atom in the chain of PDB
				final AtomType atomType = inputParameters.getAtomType();
				final String aa = inputParameters.getAa();
				Atom atom = getAtom(chain.getIdentifier(), aa, atomType, positionInPDB);
				if (atom == null) {
					log.debug("Atom is not found with chainID: '" + chain.getIdentifier() + "' in aa: '" + aa
							+ "' atomType: '" + atomType + "' position in PDB: '" + positionInPDB + "'");
					return null;
				}
				Double accesibility = getSurfaceAccesibilityOfAtom(atom, inputParameters.isRemoveOtherChains(),
						inputParameters.isRemoveOtherMolecules());
				if (accesibility != null) {
					String experimentalMethod2 = getExperimentalMethod();
					SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(dbRef.getPdbID(),
							accesibility, aa, atom, atomType, positionInPDB, inputParameters.getUniprotACC(),
							inputParameters.getPositionInUniprotProtein(), chain.getResolution(),
							inputParameters.isRemoveOtherChains(), inputParameters.isRemoveOtherMolecules(),
							getMutation(), experimentalMethod2);
					reportsByPDBPositionAndChainAndRemovals.put(key, report);
					return report;

				}
			} else {
				log.debug("Peptide " + uniprotPeptideSeq + " is not found in PDB entry " + dbRef.getPdbID() + " - "
						+ dbRef.getChainID() + " with protein sequence: " + pdbProteinSeq);
			}
		} else {
			log.debug("Chain " + chain.getIdentifier() + " not found in PDBParser " + this);
		}
		return null;
	}

	private Boolean getMutation() {
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

	private Atom getAtom(String chainID, String aa, AtomType atomType, int positionInPDB) {

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

	private Double getSurfaceAccesibilityOfAtom(Atom atom, boolean removeOtherChains, boolean removeOtherMolecules) {

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

	private List<Atom> getAtoms() {
		if (atomList == null) {
			atomList = new ArrayList<Atom>();

			if (!getLinesContaining(CA_ATOMS_ONLY).isEmpty()) {
				return atomList;
			}

			List<String> lines = getLinesStarting(ATOM);

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

	private List<Atom> getAtomsByAAPosition(int aaPosition) {
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

	public List<DBRef> getDBRefs() {
		final List<String> lines = getLinesContaining(DBREF);
		List<DBRef> list = new ArrayList<DBRef>();
		for (String dbLine : lines) {
			list.add(new DBRef(dbLine));
		}
		return list;
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

	private List<String> getLinesContaining(String containing) {

		// When filteredLines is closed, it closes underlying stream as well as
		// underlying file.
		try {
			return Files.lines(Paths.get(filePath)).filter(s -> s.contains(containing)).collect(Collectors.toList());
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
}
