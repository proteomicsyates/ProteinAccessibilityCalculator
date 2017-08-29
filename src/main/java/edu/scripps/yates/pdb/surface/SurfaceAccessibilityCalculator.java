package edu.scripps.yates.pdb.surface;

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
import edu.scripps.yates.annotations.uniprot.xml.DbReferenceType;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Chain;
import edu.scripps.yates.pdb.model.SurfacePeptide;
import edu.scripps.yates.pdb.model.SurfaceProtein;
import edu.scripps.yates.pdb.read.PDBParser;
import edu.scripps.yates.pdb.read.PDBParserManager;
import edu.scripps.yates.pdb.read.PDBUtil;
import edu.scripps.yates.pdb.util.JMolCommandsUtil;
import edu.scripps.yates.pdb.util.SurfaceAccebilityInputParameters;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.Pair;

public class SurfaceAccessibilityCalculator {
	private static final int MIN_LENTH = 6;
	private static final int MAX_LENTH = 100;
	private final UniprotProteinLocalRetriever uplr;
	private final String aa;
	private final AtomType atomType;
	private final static Logger log = Logger.getLogger(SurfaceAccessibilityCalculator.class);
	private final PDBParserManager pdbParserManager;
	private boolean createImages = false;
	private String uniprotVersion = null;
	private final boolean removeOtherChains;
	private final boolean removeOtherMolecules;
	private SurfaceAccessibilityManager manager;
	private final boolean testAllPositilities;
	private static List<Double> numPDBStructuresList = new ArrayList<Double>();

	public SurfaceAccessibilityCalculator(UniprotProteinLocalRetriever uplr, String aa, AtomType atomType,
			boolean removeOtherChains, boolean removeOtherMolecules, boolean testAllPositilities,
			File parentPDBFolderContainer) {
		this.uplr = uplr;
		this.aa = aa;
		this.atomType = atomType;
		this.removeOtherChains = removeOtherChains;
		this.removeOtherMolecules = removeOtherMolecules;
		this.testAllPositilities = testAllPositilities;
		pdbParserManager = new PDBParserManager(parentPDBFolderContainer);
	}

	public SurfaceAccessibilityProteinReport getSurfaceAccesibilityFromProtein(SurfaceProtein protein) {
		List<SurfaceProtein> list = new ArrayList<SurfaceProtein>();

		list.add(protein);
		return getSurfaceAccesibilityFromProteins(list).get(protein.getAcc());
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
	 * @param createImages
	 *            the createImages to set
	 */
	public void setCreateImages(boolean createImages) {
		this.createImages = createImages;
	}

	protected Map<String, SurfaceAccessibilityProteinReport> getSurfaceAccesibilityFromProteins(
			Collection<SurfaceProtein> proteins) {
		Set<String> uniprotAccs = new HashSet<String>();
		for (SurfaceProtein protein : proteins) {
			final Pair<String, String> acc = FastaParser.getACC(protein.getAcc());
			if (acc != null && acc.getSecondElement().equals("UNIPROT")) {
				if (acc.getFirstelement().equals("P09211")) {
					System.out.println("asdf");
				}
				uniprotAccs.add(acc.getFirstelement());
			}
		}
		Map<String, SurfaceAccessibilityProteinReport> ret = new HashMap<String, SurfaceAccessibilityProteinReport>();
		final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(getUniprotVersion(), uniprotAccs);
		for (SurfaceProtein protein : proteins) {
			if (annotatedProteins.containsKey(protein.getAcc())) {
				final SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein = getSurfaceAccesibilityFromProtein(
						protein, annotatedProteins.get(protein.getAcc()));
				if (surfaceAccesibilityFromProtein != null) {
					ret.put(protein.getAcc(), surfaceAccesibilityFromProtein);
				}
			}
		}
		return ret;
	}

	public SurfaceAccessibilityProteinReport getSurfaceAccesibilityFromPDBModel(String pdbID, String chainID) {
		List<SurfaceAccessibilityProteinReport> list = new ArrayList<SurfaceAccessibilityProteinReport>();
		if (this.testAllPositilities) {
			list.add(getSurfaceAccesibilityFromPDBModel(pdbID, chainID, true, true));
			list.add(getSurfaceAccesibilityFromPDBModel(pdbID, chainID, false, true));
			list.add(getSurfaceAccesibilityFromPDBModel(pdbID, chainID, false, false));
		} else {
			list.add(getSurfaceAccesibilityFromPDBModel(pdbID, chainID, true, true));
		}
		if (list.isEmpty()) {
			return null;
		}
		// merge to return just one
		SurfaceAccessibilityProteinReport ret = list.get(0);
		for (int i = 1; i < list.size(); i++) {
			for (SiteSurfaceAccessibilityReport report : list.get(i).getReports()) {
				ret.addSurfaceAccesibilityReport(report);
			}
		}
		return ret;
	}

	private SurfaceAccessibilityProteinReport getSurfaceAccesibilityFromPDBModel(String pdbID, String chainID,
			boolean removeOtherChains, boolean removeOtherMolecules) {
		pdbParserManager.clearParsers();

		SurfaceAccessibilityProteinReport proteinReport = new SurfaceAccessibilityProteinReport(pdbID, null);
		SurfaceAccebilityInputParameters inputParameters = null;
		PDBParser parser = pdbParserManager.getPDBParserByPDBID(pdbID);
		if (parser != null) {
			if (chainID != null) {
				inputParameters = new SurfaceAccebilityInputParameters(aa, atomType, pdbID, chainID, removeOtherChains,
						removeOtherMolecules);
			} else {
				inputParameters = new SurfaceAccebilityInputParameters(aa, atomType, pdbID, removeOtherChains,
						removeOtherMolecules);
			}
			SurfaceAccessibilityProteinReport proteinAccessibilityReportByProtein = manager
					.getProteinAccessibilityReportByProtein(inputParameters.getReportKey());
			if (proteinAccessibilityReportByProtein != null) {
				for (SiteSurfaceAccessibilityReport report : proteinAccessibilityReportByProtein.getReports()) {
					proteinReport.addSurfaceAccesibilityReport(report);
				}

			} else {
				// if not exception, exit loop here
				// break;

				log.info("Using PDB model " + pdbID);

				Map<String, SiteSurfaceAccessibilityReport> siteReports = parser
						.getSiteSurfaceAccesibilityFromParameters(inputParameters);

				if (siteReports != null) {
					for (String key : siteReports.keySet()) {
						proteinReport.addSurfaceAccesibilityReport(siteReports.get(key));
					}

				}
			}
		}
		if (!proteinReport.getReports().isEmpty()) {
			log.info("Surface accessibilities calculated for " + pdbID + " in " + proteinReport.getReports().size()
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

	private SurfaceAccessibilityProteinReport getSurfaceAccesibilityFromProtein(SurfaceProtein protein, Entry entry) {
		List<SurfaceAccessibilityProteinReport> list = new ArrayList<SurfaceAccessibilityProteinReport>();

		if (this.testAllPositilities) {
			SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein = getSurfaceAccesibilityFromProtein(
					protein, entry, true, true);
			if (surfaceAccesibilityFromProtein != null) {
				list.add(surfaceAccesibilityFromProtein);
			}
			SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein2 = getSurfaceAccesibilityFromProtein(
					protein, entry, false, true);
			if (surfaceAccesibilityFromProtein2 != null) {
				list.add(surfaceAccesibilityFromProtein2);
			}
			SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein3 = getSurfaceAccesibilityFromProtein(
					protein, entry, false, false);
			if (surfaceAccesibilityFromProtein3 != null) {
				list.add(surfaceAccesibilityFromProtein3);
			}
		} else {
			SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein = getSurfaceAccesibilityFromProtein(
					protein, entry, this.removeOtherChains, this.removeOtherMolecules);
			if (surfaceAccesibilityFromProtein != null) {
				list.add(surfaceAccesibilityFromProtein);
			}
		}
		if (list.isEmpty()) {
			return null;
		}
		// merge to return just one
		SurfaceAccessibilityProteinReport ret = list.get(0);
		for (int i = 1; i < list.size(); i++) {
			List<SiteSurfaceAccessibilityReport> reports = list.get(i).getReports();
			if (reports == null) {
				continue;
			}
			for (SiteSurfaceAccessibilityReport report : reports) {
				ret.addSurfaceAccesibilityReport(report);
			}
		}
		return ret;
	}

	private SurfaceAccessibilityProteinReport getSurfaceAccesibilityFromProtein(SurfaceProtein protein, Entry entry,
			boolean removeOtherChains, boolean removeOtherMolecules) {
		final String uniprotProteinSeq = getUniprotProteinSequence(entry);
		if (uniprotProteinSeq == null) {
			return null;
		}
		pdbParserManager.clearParsers();

		if (protein.getAcc().equals("P09211")) {
			System.out.println(removeOtherChains + "\t" + removeOtherMolecules);
		}
		final Set<SurfacePeptide> uniprotPeptides = protein.getPeptides();// getPeptides(uniprotProteinSeq);
		final String acc = protein.getAcc();
		SurfaceAccessibilityProteinReport proteinReport = new SurfaceAccessibilityProteinReport(acc, uniprotProteinSeq);
		Set<Integer> positionsInUniprotProteinProcessed = new HashSet<Integer>();
		Map<String, List<SurfaceAccebilityInputParameters>> parametersByPDBID = new HashMap<String, List<SurfaceAccebilityInputParameters>>();
		int numCalculationsToDo = 0;
		for (SurfacePeptide uniprotPeptide : uniprotPeptides) {
			String uniprotPeptideSequence = uniprotPeptide.getSequence();
			List<Integer> peptidePositionInUniprots = StringUtils.allPositionsOf(uniprotProteinSeq,
					uniprotPeptideSequence);
			if (peptidePositionInUniprots.size() > 1) {
				log.debug(uniprotPeptideSequence + "_ is present " + peptidePositionInUniprots.size()
						+ " times in protein " + acc);
			}
			for (Integer peptidePositionInUniprot : peptidePositionInUniprots) {
				for (int positionInPeptide = 0; positionInPeptide < uniprotPeptideSequence
						.length(); positionInPeptide++) {
					if (uniprotPeptideSequence.charAt(positionInPeptide) == aa.charAt(0)) {
						final int positionInUniprotProtein = peptidePositionInUniprot + positionInPeptide;
						if (!positionsInUniprotProteinProcessed.contains(positionInUniprotProtein)) {
							positionsInUniprotProteinProcessed.add(positionInUniprotProtein);
							// if it is already done, don't do it
							if (proteinReport.containsReportsForPosition(positionInUniprotProtein)) {
								// continue;
							}
							// get a list of ranked chains from the best to the
							// worst proteinPDB structure
							final List<Chain> chains = getPDBChainListSortedByResolution(entry,
									positionInUniprotProtein);
							if (chains.isEmpty()) {
								log.debug("No PDB structures for position " + positionInUniprotProtein + " protein "
										+ acc);
							}
							for (Chain chain : chains) {
								if (chain.getPdbID().equals("10GS")) {
									System.out.println(chain);
								}
								PDBParser parser = pdbParserManager.getPDBParserByPDBID(chain.getPdbID());
								if (parser != null) {
									if (!parser.containsUniprotReference(acc)) {
										continue;
									}

									SurfaceAccebilityInputParameters inputParameters = new SurfaceAccebilityInputParameters(
											aa, atomType, positionInUniprotProtein, chain, uniprotPeptideSequence,
											positionInPeptide, acc, removeOtherChains, removeOtherMolecules);
									SurfaceAccessibilityProteinReport proteinAccessibilityReportByProtein = manager
											.getProteinAccessibilityReportByProtein(inputParameters.getReportKey());
									if (proteinAccessibilityReportByProtein != null) {
										for (SiteSurfaceAccessibilityReport report : proteinAccessibilityReportByProtein
												.getReports()) {
											proteinReport.addSurfaceAccesibilityReport(report);
										}

									} else {
										numCalculationsToDo++;
										if (parametersByPDBID.containsKey(inputParameters.getPdbID())) {
											parametersByPDBID.get(inputParameters.getPdbID()).add(inputParameters);
										} else {
											List<SurfaceAccebilityInputParameters> list = new ArrayList<SurfaceAccebilityInputParameters>();
											list.add(inputParameters);
											parametersByPDBID.put(inputParameters.getPdbID(), list);
										}
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
		if (protein.getAcc().equals("P09211")) {
			System.out.println(removeOtherChains + "\t" + removeOtherMolecules);
		}
		if (numCalculationsToDo > 0) {
			log.info(numCalculationsToDo + " calculations to do for " + acc);
		}
		for (String pdbID : parametersByPDBID.keySet()) {
			if (protein.getAcc().equals("P09211") && pdbID.equals("10GS")) {
				log.info("adsf");
			}
			PDBParser parser = pdbParserManager.getPDBParserByPDBID(pdbID);
			log.info("Using PDB model " + pdbID + " for protein " + acc);
			if (parser != null) {
				if (!parser.containsUniprotReference(acc)) {
					continue;
				}
				List<SurfaceAccebilityInputParameters> list = parametersByPDBID.get(pdbID);
				sortByUniprotPosition(list);
				for (SurfaceAccebilityInputParameters surfaceAccebilityInputParameters : list) {
					log.debug("Using PDB model " + pdbID + " for protein " + acc + " position "
							+ surfaceAccebilityInputParameters.getPositionInUniprotProtein());
					Map<String, SiteSurfaceAccessibilityReport> siteReports = parser
							.getSiteSurfaceAccesibilityFromParameters(surfaceAccebilityInputParameters);
					if (siteReports != null) {
						for (String key : siteReports.keySet()) {
							proteinReport.addSurfaceAccesibilityReport(siteReports.get(key));
						}

					}
				}

			}
		}
		if (protein.getAcc().equals("P09211")) {
			System.out.println(removeOtherChains + "\t" + removeOtherMolecules);
		}
		if (!proteinReport.getReports().isEmpty() || numCalculationsToDo > 0) {
			log.info("Surface accessibilities calculated for " + acc + " in " + proteinReport.getReports().size()
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

	private void sortByUniprotPosition(List<SurfaceAccebilityInputParameters> list) {
		Collections.sort(list, new Comparator<SurfaceAccebilityInputParameters>() {

			@Override
			public int compare(SurfaceAccebilityInputParameters o1, SurfaceAccebilityInputParameters o2) {
				return Integer.compare(o1.getPositionInUniprotProtein(), o2.getPositionInUniprotProtein());
			}
		});
	}

	public List<File> getJPGImages(SurfaceAccessibilityProteinReport proteinReport) {
		List<File> files = new ArrayList<File>();

		// do this per each individual pdb
		final Map<String, Set<SiteSurfaceAccessibilityReport>> reportsByPDBID = proteinReport.getReportsByPDBID();
		for (String pdbID : reportsByPDBID.keySet()) {

			final JMolScript selectAndAddLabels = JMolCommandsUtil.getSelectAndAddLabels(proteinReport, pdbID);
			final PDBParser pdbParser = pdbParserManager.getPDBParserByPDBID(pdbID);
			pdbParser.executeCommands(selectAndAddLabels);

			File file = new File(pdbParser.getFileFolder().getAbsolutePath() + File.separator
					+ proteinReport.getUniprotACC() + "_" + pdbID + ".png");
			JMolScript saveImageCommand = new JMolScript(
					"write image 600 600 PNG 2 \"" + file.getAbsolutePath().replace(File.separator, "/") + "\"");
			// String saveImageCommand = "x = write(\"PNGJ\");";
			pdbParser.executeCommands(saveImageCommand);

			files.add(file);
		}
		return files;

	}

	public static String getUniprotProteinSequence(Entry entry) {
		if (entry.getSequence() != null) {
			return entry.getSequence().getValue().replace("\n", "");
		}
		log.warn("UniprotKB entry " + entry.getAccession().get(0) + " has no protein sequence");
		return null;
	}

	/**
	 *
	 * @param entry
	 * @param positionInUniprot
	 * @return a List of Chain. Sorted by resolution.
	 */
	private List<Chain> getPDBChainListSortedByResolution(Entry entry, int positionInUniprot) {
		List<Chain> ret = new ArrayList<Chain>();
		double numPDBStructures = 0;
		final List<DbReferenceType> dbReferences = entry.getDbReference();
		if (dbReferences != null) {
			for (DbReferenceType dbReferenceType : dbReferences) {
				if (dbReferenceType.getType().equals("PDB")) {
					numPDBStructures++;
					// see if this crystal structure include the position we are
					// interested in
					String resolutionString = PDBUtil.getPropertyValueFromDbReferenceType(dbReferenceType,
							"resolution");
					Float resolution = null;
					if (resolutionString != null) {
						resolution = Float.valueOf(resolutionString);
					}
					String chainsString = PDBUtil.getPropertyValueFromDbReferenceType(dbReferenceType, "chains");
					Set<Chain> chains = getChainsFromChainString(chainsString, dbReferenceType, resolution);
					for (Chain chain : chains) {
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
			Comparator<Chain> comparator = new Comparator<Chain>() {
				@Override
				public int compare(Chain o1, Chain o2) {
					Float res1 = o1.getResolution() != null ? o1.getResolution() : -1;
					Float res2 = o2.getResolution() != null ? o2.getResolution() : -1;
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
		StringBuilder sb = new StringBuilder();
		sb.append(numPDBStructuresList.size() + " proteins analyzed\n");
		double mean = Maths.mean(numPDBStructuresList.toArray(new Double[0]));
		double stddev = Maths.stddev(numPDBStructuresList.toArray(new Double[0]));
		double sum = Maths.sum(numPDBStructuresList.toArray(new Double[0]));
		double max = Maths.max(numPDBStructuresList.toArray(new Double[0]));
		sb.append(sum + " total PDB structures matched\n" + mean + "(" + stddev
				+ ") structures matched per protein in average (stdev). Max=" + max);
		return sb.toString();
	}

	public static void clearStatistics() {
		numPDBStructuresList.clear();
	}

	private Set<Chain> getChainsFromChainString(String chainsString, DbReferenceType dbReferenceType,
			Float resolution) {
		Set<Chain> ret = new HashSet<Chain>();
		if (chainsString.contains(",")) {
			final String[] split = chainsString.split(",");
			for (String individualChaingString : split) {
				ret.addAll(getChainsFromIndividualChainString(individualChaingString, dbReferenceType, resolution));
			}
		} else {
			ret.addAll(getChainsFromIndividualChainString(chainsString, dbReferenceType, resolution));
		}
		return ret;
	}

	private Set<Chain> getChainsFromIndividualChainString(String individualChaingString,
			DbReferenceType dbReferenceType, Float resolution) {
		Set<Chain> ret = new HashSet<Chain>();
		if (individualChaingString.contains("=")) {
			final String[] split = individualChaingString.split("=");
			String chainIdString = split[0];
			String positionsString = split[1];
			if (chainIdString.contains("/")) {
				final String[] split2 = chainIdString.split("/");
				for (String individualChainId : split2) {
					String text = individualChainId + "=" + positionsString;
					ret.add(new Chain(dbReferenceType.getId(), text, resolution));
				}
			} else {
				String text = chainIdString + "=" + positionsString;
				ret.add(new Chain(dbReferenceType.getId(), text, resolution));
			}
		}
		return ret;
	}

	// private List<com.compomics.util.protein.Protein> getPeptides(String
	// proteinSequence) {
	// com.compomics.util.protein.Protein fastaProtein = new
	// com.compomics.util.protein.Protein("header",
	// proteinSequence);
	// List<com.compomics.util.protein.Protein> list = new
	// ArrayList<com.compomics.util.protein.Protein>();
	// final com.compomics.util.protein.Protein[] cleave =
	// enzyme.cleave(fastaProtein, MIN_LENTH, MAX_LENTH);
	// for (com.compomics.util.protein.Protein pep : cleave) {
	// list.add(pep);
	// }
	// return list;
	//
	// }

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

	public void setManager(SurfaceAccessibilityManager surfaceAccessibilityManager) {
		this.manager = surfaceAccessibilityManager;

	}

}
