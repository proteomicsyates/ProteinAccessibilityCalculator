package edu.scripps.yates.pdb;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.util.DigestionConfiguration;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Chain;
import edu.scripps.yates.pdb.model.DBRef;
import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.read.PDBParser;
import edu.scripps.yates.pdb.read.PDBParserManager;
import edu.scripps.yates.pdb.read.PDBUtil;
import edu.scripps.yates.pdb.util.InputParameters;
import edu.scripps.yates.pdb.util.JMolCommandsUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.DbReferenceType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import gnu.trove.set.hash.TIntHashSet;

public abstract class Calculator<R extends ProteinReport<T>, T extends JMolAtomReport> {
	protected final UniprotProteinLocalRetriever uplr;

	protected final Map<Character, List<AtomType>> atomTypeMap;
	protected final static Logger log = Logger.getLogger(Calculator.class);
	protected final PDBParserManager pdbParserManager;
	protected boolean createImages = false;
	protected String uniprotVersion = null;
	protected boolean removeOtherChains;
	protected boolean removeOtherMolecules;
	protected JMolReportManager<R, T> reportManager;
	protected static TDoubleArrayList numPDBStructuresList = new TDoubleArrayList();
	protected final Map<String, T> reportsByPDBPositionAndChainAndRemovals = new THashMap<String, T>();
	private DigestionConfiguration fastaDigestionConfiguration;
	private boolean oneModelPerProtein;

	public Calculator(Map<Character, List<AtomType>> atomTypeMap, boolean removeOtherChains,
			boolean removeOtherMolecules, File parentPDBFolderContainer) {
		this(null, atomTypeMap, removeOtherChains, removeOtherMolecules, false, parentPDBFolderContainer);

	}

	public Calculator(UniprotProteinLocalRetriever uplr, Map<Character, List<AtomType>> atomTypeMap,
			boolean removeOtherChains, boolean removeOtherMolecules, boolean oneModelPerProtein,
			File parentPDBFolderContainer) {
		this.uplr = uplr;

		this.atomTypeMap = atomTypeMap;
		this.removeOtherChains = removeOtherChains;
		this.removeOtherMolecules = removeOtherMolecules;
		this.oneModelPerProtein = oneModelPerProtein;
		pdbParserManager = new PDBParserManager(parentPDBFolderContainer);

	}

	public Calculator(UniprotProteinLocalRetriever uplr, Character aa, AtomType atomType, boolean removeOtherChains,
			boolean removeOtherMolecules, boolean oneModelPerProtein, File parentPDBFolderContainer) {
		this.uplr = uplr;

		this.atomTypeMap = new HashMap<Character, List<AtomType>>();
		final List<AtomType> list = new ArrayList<AtomType>();
		list.add(atomType);
		atomTypeMap.put(aa, list);
		this.removeOtherChains = removeOtherChains;
		this.removeOtherMolecules = removeOtherMolecules;
		this.oneModelPerProtein = oneModelPerProtein;
		pdbParserManager = new PDBParserManager(parentPDBFolderContainer);

	}

	public String getUniprotVersion() {
		return uniprotVersion;
	}

	public void setUniprotVersion(String uniprotVersion) {
		this.uniprotVersion = uniprotVersion;
	}

	/**
	 * @return the createImages
	 */
	public boolean isCreateImages() {
		return createImages;
	}

	/**
	 * @param createImages the createImages to set
	 */
	public void setCreateImages(boolean createImages) {
		this.createImages = createImages;
	}

	protected Map<String, R> getReportFromProteins(Collection<Protein> proteins) {
		final Set<String> uniprotAccs = new THashSet<String>();
		for (final Protein protein : proteins) {
			final Accession acc = FastaParser.getACC(protein.getAcc());
			if (acc != null && acc.getAccessionType() == AccessionType.UNIPROT) {
				uniprotAccs.add(acc.getAccession());
			}
		}
		final Map<String, R> ret = new THashMap<String, R>();
		final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(getUniprotVersion(), uniprotAccs);
		for (final Protein protein : proteins) {
			if (annotatedProteins.containsKey(protein.getAcc())) {
				final R surfaceAccesibilityFromProtein = getReportFromProtein(protein, null,
						annotatedProteins.get(protein.getAcc()));
				if (surfaceAccesibilityFromProtein != null) {
					ret.put(protein.getAcc(), surfaceAccesibilityFromProtein);
				}
			}
		}
		return ret;
	}

	public R getReportFromPDBModel(String pdbID) {
		return getReportFromPDBModel(pdbID, null);
	}

	public R getReportFromPDBModel(String pdbID, String chainID) {

		pdbParserManager.clearParsers();

		final R proteinReport = createProteinReportObject(pdbID, null);
		InputParameters inputParameters = null;
		final PDBParser parser = pdbParserManager.getPDBParserByPDBID(pdbID, isParseCoordinates());
		if (parser != null) {
			if (chainID != null) {
				inputParameters = new InputParameters(atomTypeMap, pdbID, chainID, removeOtherChains,
						removeOtherMolecules);
			} else {
				inputParameters = new InputParameters(atomTypeMap, pdbID, removeOtherChains, removeOtherMolecules);
			}

			// if not exception, exit loop here
			// break;

			log.info("Using PDB model " + pdbID);

			final Map<String, T> siteReports = getSiteReportFromParameters(parser, inputParameters);

			if (siteReports != null) {
				for (final String key : siteReports.keySet()) {
					proteinReport.addReport(siteReports.get(key));
				}

			}

		}
		if (!proteinReport.getReports().isEmpty()) {
			log.info(" accessibilities calculated for " + pdbID + " in " + proteinReport.getReports().size()
					+ " different sites and structure models");
		}
		if (createImages) {
			getJPGImages(proteinReport);
		}
		if (proteinReport.isEmpty()) {
			return null;
		}
		return proteinReport;

	}

	public abstract R createProteinReportObject(String acc, String proteinSequence);

	private void sortByUniprotPosition(List<InputParameters> list) {
		Collections.sort(list, new Comparator<InputParameters>() {

			@Override
			public int compare(InputParameters o1, InputParameters o2) {
				return Integer.compare(o1.getPositionInUniprotProtein(), o2.getPositionInUniprotProtein());
			}
		});
	}

	public List<File> getJPGImages(R proteinReport) {
		final List<File> files = new ArrayList<File>();

		// do this per each individual pdb
		final Map<String, Set<T>> reportsByPDBID = proteinReport.getReportsByPDBID();
		for (final String pdbID : reportsByPDBID.keySet()) {

			final JMolScript selectAndAddLabels = JMolCommandsUtil.getSelectAndAddLabels(proteinReport, pdbID);
			final PDBParser pdbParser = pdbParserManager.getPDBParserByPDBID(pdbID, isParseCoordinates());
			pdbParser.executeCommands(selectAndAddLabels);

			final File file = new File(pdbParser.getFileFolder().getAbsolutePath() + File.separator
					+ proteinReport.getUniprotACC() + "_" + pdbID + ".png");
			final JMolScript saveImageCommand = new JMolScript(
					"write image 600 600 PNG 2 \"" + file.getAbsolutePath().replace(File.separator, "/") + "\"");
			// String saveImageCommand = "x = write(\"PNGJ\");";
			pdbParser.executeCommands(saveImageCommand);

			files.add(file);
		}
		return files;

	}

	public static String getUniprotProteinSequence(Entry entry) {
		if (entry != null) {
			if (entry.getSequence() != null) {
				return entry.getSequence().getValue().replace("\n", "");
			}
			log.warn("UniprotKB entry " + entry.getAccession().get(0) + " has no protein sequence");

		}
		return null;
	}

	/**
	 *
	 * @param entry
	 * @param positionInUniprot
	 * @return a List of Chain. Sorted by resolution.
	 */
	protected List<Chain> getPDBChainListSortedByResolution(Entry entry, String pdbID, int positionInUniprot) {
		final List<Chain> ret = new ArrayList<Chain>();
		double numPDBStructures = 0;
		final List<DbReferenceType> dbReferences = entry.getDbReference();
		if (dbReferences != null) {
			for (final DbReferenceType dbReferenceType : dbReferences) {
				if (dbReferenceType.getType().equals("PDB")) {
					numPDBStructures++;
					// see if this crystal structure include the position we are
					// interested in
					String resolutionString = PDBUtil.getPropertyValueFromDbReferenceType(dbReferenceType,
							"resolution");
					Float resolution = null;
					if (resolutionString != null) {
						if (resolutionString.equals("2.70 A")) {
							log.info(PDBUtil.getPropertyValueFromDbReferenceType(dbReferenceType, "resolution"));
						}
						try {
							resolution = Float.valueOf(resolutionString);
						} catch (final NumberFormatException e) {
							if (resolutionString.endsWith(" A")) {
								resolutionString = resolutionString.split(" ")[0];
								try {
									resolution = Float.valueOf(resolutionString);
								} catch (final NumberFormatException e2) {

								}
							}
						}
					}
					final String chainsString = PDBUtil.getPropertyValueFromDbReferenceType(dbReferenceType, "chains");
					final Set<Chain> chains = getChainsFromChainString(chainsString, dbReferenceType, resolution);
					for (final Chain chain : chains) {
						// if pdbID is provided, only take that model
						if (pdbID != null && !chain.getPdbID().equals(pdbID)) {
							continue;
						}
						if (chain.includesPosition(positionInUniprot)) {
							ret.add(chain);
						}
					}

				}
			}
		}
		numPDBStructuresList.add(numPDBStructures);
		if (!ret.isEmpty()) {
			// sort by resolution
			final Comparator<Chain> comparator = new Comparator<Chain>() {
				@Override
				public int compare(Chain o1, Chain o2) {
					final Float res1 = o1.getResolution() != null ? o1.getResolution() : -1;
					final Float res2 = o2.getResolution() != null ? o2.getResolution() : -1;
					return Float.compare(res2, res1);
				}
			};
			Collections.sort(ret, comparator);
			// log.debug(ret.size() + " models for protein " +
			// entry.getAccession().get(0) + " site " + positionInUniprot);
			// log.debug("Best model for position " + positionInUniprot + " in "
			// + entry.getAccession().get(0) + " is "
			// + ret.get(0));
		}
		return ret;
	}

	public static String getStatistics() {
		final StringBuilder sb = new StringBuilder();
		sb.append(numPDBStructuresList.size() + " proteins analyzed\n");

		final double mean = Maths.mean(numPDBStructuresList);
		final double stddev = Maths.stddev(numPDBStructuresList);
		final double sum = numPDBStructuresList.sum();
		final double max = numPDBStructuresList.max();
		sb.append(sum + " total PDB structures matched\n" + mean + "(" + stddev
				+ ") structures matched per protein in average (stdev). Max=" + max);
		return sb.toString();
	}

	public static void clearStatistics() {
		numPDBStructuresList.clear();
	}

	private Set<Chain> getChainsFromChainString(String chainsString, DbReferenceType dbReferenceType,
			Float resolution) {
		if (chainsString == null) {
			return Collections.emptySet();
		}
		final Set<Chain> ret = new THashSet<Chain>();
		if (chainsString.contains(",")) {
			final String[] split = chainsString.split(",");
			for (final String individualChaingString : split) {
				ret.addAll(getChainsFromIndividualChainString(individualChaingString, dbReferenceType, resolution));
			}
		} else {
			ret.addAll(getChainsFromIndividualChainString(chainsString, dbReferenceType, resolution));
		}
		return ret;
	}

	private Set<Chain> getChainsFromIndividualChainString(String individualChaingString,
			DbReferenceType dbReferenceType, Float resolution) {
		final Set<Chain> ret = new THashSet<Chain>();
		if (individualChaingString.contains("=")) {
			final String[] split = individualChaingString.split("=");
			final String chainIdString = split[0];
			final String positionsString = split[1];
			if (chainIdString.contains("/")) {
				final String[] split2 = chainIdString.split("/");
				for (final String individualChainId : split2) {
					final String text = individualChainId + "=" + positionsString;
					ret.add(new Chain(dbReferenceType.getId(), text, resolution));
				}
			} else {
				final String text = chainIdString + "=" + positionsString;
				ret.add(new Chain(dbReferenceType.getId(), text, resolution));
			}
		}
		return ret;
	}

	/**
	 * @return the pdbParserManager
	 */
	public PDBParserManager getPdbParserManager() {
		return pdbParserManager;
	}

	/**
	 * @return the uplr
	 */
	public UniprotProteinLocalRetriever getUplr() {
		return uplr;
	}

	public boolean isRemoveOtherChains() {
		return removeOtherChains;
	}

	public boolean isRemoveOtherMolecules() {
		return removeOtherMolecules;
	}

	public void setManager(JMolReportManager<R, T> reportManager) {
		this.reportManager = reportManager;

	}

	public Map<String, T> getSiteReportFromParameters(PDBParser parser, InputParameters inputParameters) {
		if (inputParameters.getUniprotACC() != null) {
			final Map<String, T> ret = new THashMap<String, T>();
			final T report = getSiteReportMappedToUniprot(parser, inputParameters);
			if (report != null) {
				ret.put(report.getAtom().getPositionInPDB() + report.getAtom().getChainID(), report);
			}
			return ret;
		} else {
			return getSiteReportForPDBModel(parser, inputParameters);
		}
	}

	private Map<String, T> getSiteReportForPDBModel(PDBParser parser, InputParameters inputParameters) {
		// get the proteinSequence in PDB
		Chain chain = inputParameters.getChain();
		final List<DBRef> dbRefs = new ArrayList<DBRef>();
		if (chain != null) {
			final DBRef dbRef = PDBUtil.getDBRef(parser, chain.getIdentifier());
			dbRefs.add(dbRef);
		} else {
			dbRefs.addAll(parser.getDBRefs());
		}
		for (final DBRef dbRef : dbRefs) {
			if (dbRef != null && dbRef.getChainID() != null) {
				chain = new Chain(parser.getPdbID(), dbRef.getChainID() + "=0-0", -1.0f);
				final String pdbProteinSeq = parser.getSequence(dbRef);
				if ("".equals(pdbProteinSeq)) {
					log.warn("Protein sequence for chain " + dbRef.getChainID() + " was not obtained");
					continue;
				}
				for (final Character aa : inputParameters.getAtomTypeMap().keySet()) {
					final TIntArrayList aaPositionsInPDB = StringUtils.allPositionsOf(pdbProteinSeq, aa);
					for (final int positionInPDB : aaPositionsInPDB.toArray()) {

						final String key = positionInPDB + chain.getIdentifier() + inputParameters.isRemoveOtherChains()
								+ inputParameters.isRemoveOtherMolecules();

						if (reportsByPDBPositionAndChainAndRemovals.containsKey(key)) {
							log.info("Surface  accesibility already calculated for site " + positionInPDB + " in chain "
									+ chain.getIdentifier() + " of PDB " + chain.getPdbID()
									+ " with removingOtherChains=" + inputParameters.isRemoveOtherChains()
									+ " and removingOtherMolcules=" + inputParameters.isRemoveOtherMolecules());
							continue;
						}
						// get the appropriate atom in the chain of PDB

						final List<AtomType> atomTypeList = inputParameters.getAtomTypeMap().get(aa);
						for (final AtomType atomType : atomTypeList) {
							final Atom3D atom = parser.getAtom(chain.getIdentifier(), aa, atomType, positionInPDB);
							if (atom == null) {
								log.warn("Atom is not found in chainID: '" + chain.getIdentifier() + "' in aa: '" + aa
										+ "' atomType: '" + atomType + "' position in PDB: '" + positionInPDB + "'");
								continue;
							}
							final T report = calculateReport(parser, inputParameters, atom, positionInPDB,
									chain.getResolution());
							if (report != null) {
								reportsByPDBPositionAndChainAndRemovals.put(key, report);
							}
						}
					}

				}
			} else {
				log.debug("Chain " + chain.getIdentifier() + " not found in PDBParser " + this);
			}
		}
		return reportsByPDBPositionAndChainAndRemovals;
	}

	public abstract T calculateReport(PDBParser parser, InputParameters inputParameters, Atom3D atom, int positionInPDB,
			Float resolution);

	private T getSiteReportMappedToUniprot(PDBParser parser, InputParameters inputParameters) {
		// get the proteinSequence in PDB
		final Chain chain = inputParameters.getChain();
		final DBRef dbRef = PDBUtil.getDBRef(parser, chain.getIdentifier());
		if (dbRef != null && dbRef.getChainID() != null) {
			final String pdbProteinSeq = parser.getSequence(dbRef);
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
					log.debug("Surface accesibility already calculated for site " + positionInPDB + " in chain "
							+ chain.getIdentifier() + " of PDB " + chain.getPdbID() + " with removingOtherChains="
							+ inputParameters.isRemoveOtherChains() + " and removingOtherMolcules="
							+ inputParameters.isRemoveOtherMolecules());
					return reportsByPDBPositionAndChainAndRemovals.get(key);
				}
				// get the appropriate atom in the chain of PDB
				for (final Character aa : inputParameters.getAtomTypeMap().keySet()) {
					final List<AtomType> atomTypeList = inputParameters.getAtomTypeMap().get(aa);
					for (final AtomType atomType : atomTypeList) {

						final Atom3D atom = parser.getAtom(chain.getIdentifier(), aa, atomType, positionInPDB);
						if (atom == null) {
							log.debug("Atom is not found with chainID: '" + chain.getIdentifier() + "' in aa: '" + aa
									+ "' atomType: '" + atomType + "' position in PDB: '" + positionInPDB + "'");
							return null;
						}
						final T report = calculateReport(parser, inputParameters, atom, positionInPDB,
								chain.getResolution());
						if (report != null) {
							reportsByPDBPositionAndChainAndRemovals.put(key, report);
							return report;
						}
					}
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

	public R getReportFromProtein(String proteinAcc, String pdbID, Entry entry) {
		if (proteinAcc != null) {
			if (entry == null) {
				final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(uniprotVersion, proteinAcc);
				if (annotatedProtein.containsKey(proteinAcc)) {
					entry = annotatedProtein.get(proteinAcc);
				}
			}
			final String uniprotProteinSeq = getUniprotProteinSequence(entry);
			if (uniprotProteinSeq == null) {
				return null;
			}
			if (fastaDigestionConfiguration != null) {
				final List<String> peptides = fastaDigestionConfiguration.digestProtein(uniprotProteinSeq);
				return getReportFromProtein(proteinAcc, pdbID, peptides, entry);
			}
		} else {
			return getReportFromPDBModel(pdbID, null);
		}
		return null;
	}

	public R getReportFromProtein(Protein protein, String pdbID, Entry entry) {
		final Set<String> peptideSequences = new HashSet<String>();
		for (final Peptide peptide : protein.getPeptides()) {
			peptideSequences.add(peptide.getSequence());
		}
		return getReportFromProtein(protein.getAcc(), pdbID, peptideSequences, entry);
	}

	public R getReportFromProtein(String proteinAcc, String pdbID, Collection<String> peptides, Entry entry) {
		if (entry == null) {
			final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(uniprotVersion, proteinAcc);
			if (annotatedProtein.containsKey(proteinAcc)) {
				entry = annotatedProtein.get(proteinAcc);
			}
		}
		final String uniprotProteinSeq = getUniprotProteinSequence(entry);
		if (uniprotProteinSeq == null) {
			return null;
		}
		pdbParserManager.clearParsers();

		final R proteinReport = createProteinReportObject(proteinAcc, uniprotProteinSeq);
		final TIntHashSet positionsInUniprotProteinProcessed = new TIntHashSet();
		final Map<String, List<InputParameters>> parametersByPDBID = new THashMap<String, List<InputParameters>>();
		int numCalculationsToDo = 0;
		for (final String uniprotPeptideSequence : peptides) {

			final TIntArrayList peptidePositionInUniprots = StringUtils.allPositionsOf(uniprotProteinSeq,
					uniprotPeptideSequence);
			if (peptidePositionInUniprots.size() > 1) {
				log.debug(uniprotPeptideSequence + "_ is present " + peptidePositionInUniprots.size()
						+ " times in protein " + proteinAcc);
			}

			for (final int peptidePositionInUniprot : peptidePositionInUniprots.toArray()) {
				for (int positionInPeptide = 0; positionInPeptide < uniprotPeptideSequence
						.length(); positionInPeptide++) {
					// iterate over aas
					for (final Character aa : atomTypeMap.keySet()) {
						if (uniprotPeptideSequence.charAt(positionInPeptide) == aa) {
							final int positionInUniprotProtein = peptidePositionInUniprot + positionInPeptide;
							if (!positionsInUniprotProteinProcessed.contains(positionInUniprotProtein)) {
								positionsInUniprotProteinProcessed.add(positionInUniprotProtein);
								// if it is already done, don't do it again
								if (proteinReport.containsReportsForPosition(positionInUniprotProtein)) {
									continue;
								}
								// get a list of ranked chains from the best to
								// the
								// worst proteinPDB structure
								final List<Chain> chains = getPDBChainListSortedByResolution(entry, pdbID,
										positionInUniprotProtein);
								if (chains.isEmpty()) {
									log.debug("No PDB structures for position " + positionInUniprotProtein + " protein "
											+ proteinAcc);
								}
								for (final Chain chain : chains) {
									final PDBParser parser = pdbParserManager.getPDBParserByPDBID(chain.getPdbID(),
											isParseCoordinates());
									if (parser != null) {
										if (!parser.containsUniprotReference(proteinAcc)) {
											continue;
										}

										final InputParameters inputParameters = new InputParameters(atomTypeMap,
												positionInUniprotProtein, chain, uniprotPeptideSequence,
												positionInPeptide, proteinAcc, removeOtherChains, removeOtherMolecules);
										final String reportKey = inputParameters.getReportKey();
										final R proteinReportByProtein = reportManager
												.getReportByKey(reportKey);
										if (proteinReportByProtein != null) {
											for (final T report : proteinReportByProtein.getReports()) {
												proteinReport.addReport(report);
											}

										} else {
											numCalculationsToDo++;
											log.info(inputParameters);
											if (parametersByPDBID.containsKey(inputParameters.getPdbID())) {
												parametersByPDBID.get(inputParameters.getPdbID()).add(inputParameters);
											} else {
												final List<InputParameters> list = new ArrayList<InputParameters>();
												list.add(inputParameters);
												parametersByPDBID.put(inputParameters.getPdbID(), list);
											}
										}
										if (oneModelPerProtein) {
											break;
										}
									}
									// if not exception, exit loop here
									// break;
								}

							}
						}
					}
				}
			}
		}

		if (numCalculationsToDo > 0) {
			log.info(numCalculationsToDo + " calculations to do for " + proteinAcc);
		}
		for (final String pdbID2 : parametersByPDBID.keySet()) {
			final PDBParser parser = pdbParserManager.getPDBParserByPDBID(pdbID2, isParseCoordinates());
			log.info("Using PDB model " + pdbID2 + " for protein " + proteinAcc);
			if (parser != null) {
				if (!parser.containsUniprotReference(proteinAcc)) {
					continue;
				}
				final List<InputParameters> list = parametersByPDBID.get(pdbID2);
				sortByUniprotPosition(list);
				for (final InputParameters inputParameters : list) {
					log.debug("Using PDB model " + pdbID2 + " for protein " + proteinAcc + " position "
							+ inputParameters.getPositionInUniprotProtein());
					final Map<String, T> siteReports = getSiteReportFromParameters(parser, inputParameters);
					if (siteReports != null) {
						for (final String key : siteReports.keySet()) {
							final T report = siteReports.get(key);
							proteinReport.addReport(report);

						}

					}
				}

			}
		}

		if (!proteinReport.getReports().isEmpty() || numCalculationsToDo > 0) {
			log.info("Surface accessibilities calculated for " + proteinAcc + " in " + proteinReport.getReports().size()
					+ " different sites and structure models");
		}
		if (createImages) {
			getJPGImages(proteinReport);
		}
		if (proteinReport.isEmpty()) {
			return null;
		}
		this.reportManager.dumpToFile();
		return proteinReport;

	}

	public abstract boolean isParseCoordinates();

	public abstract CalculationType getCalculationType();

	public DigestionConfiguration getDigestionConfiguration() {
		return fastaDigestionConfiguration;
	}

	public void setDigestionConfiguration(DigestionConfiguration digestion) {
		this.fastaDigestionConfiguration = digestion;
	}
}
