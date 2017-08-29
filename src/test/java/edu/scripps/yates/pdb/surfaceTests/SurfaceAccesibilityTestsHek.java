package edu.scripps.yates.pdb.surfaceTests;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.poi.ss.usermodel.Row;
import org.junit.Test;

import edu.scripps.yates.annotations.go.GORetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinRetrievalSettings;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.excel.ExcelReader;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.SurfacePeptide;
import edu.scripps.yates.pdb.model.SurfaceProtein;
import edu.scripps.yates.pdb.surface.SiteSurfaceAccessibilityReport;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityCalculator;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityManager;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityProteinReport;
import edu.scripps.yates.pdb.util.SurfaceAccesibilityReportWriter;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.jaxb.xpathquery.JAXBXPathQuery;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.factories.RatioEx;
import edu.scripps.yates.utilities.strings.StringUtils;
import psidev.psi.tools.ontology_manager.interfaces.OntologyTermI;

public class SurfaceAccesibilityTestsHek {
	private boolean removeOtherChains = true;
	private boolean removeOtherMolecules = true;

	@Test
	public void surfaceAccesibilityFromRemoteFiles_DNABinding_FromTable() {
		String GO = "GO:0003677";
		boolean justTheBestPerPSM = true;
		boolean filterByAnnotation = false;
		boolean testAllPosibilities = false;
		Double minRatio = -Double.MAX_VALUE;
		// Double minRatio = 0.0;
		// this gets the experimental data from a table
		OntologyTermI goTerm = new GORetriever(new File("z:\\share\\Salva\\data\\go")).getGOTermByID(GO);
		String goName = "";
		if (goTerm != null) {
			goName = goTerm.getPreferredName();
		}
		String bestPeptide = "";
		if (justTheBestPerPSM) {
			bestPeptide = "_BestCorrelatingStructure";
		}
		File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_" + goName
				+ "_ratio_gte=" + minRatio + bestPeptide + ".csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		File experimentalDataFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\experimental_proteins.csv");
		try {
			Set<String> psmIds = new HashSet<String>();
			Set<String> peptideSequences = new HashSet<String>();
			Set<String> peptideSequencesValid = new HashSet<String>();
			Set<String> uniquePositionsValid = new HashSet<String>();
			Set<String> proteinsWithNoPDBInformation = new HashSet<String>();
			Set<String> totalProteins = new HashSet<String>();
			Set<SurfacePeptide> surfacePeptides = new HashSet<SurfacePeptide>();
			List<String> accs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 0, true);
			List<String> seqs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 1, true);
			List<String> ratios = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 2, true);
			Map<String, SurfaceProtein> surfaceProteinMap = new HashMap<String, SurfaceProtein>();
			Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);

			for (int i = 0; i < accs.size(); i++) {

				Double ratio = Double.valueOf(ratios.get(i));
				if (ratio >= minRatio) {
					String acc = accs.get(i);
					if (acc.equals("P09211")) {
						System.out.println(acc);
					}
					Entry entry = annotatedProteins.get(acc);
					if (entry != null) {
						if (acc.equals("P09211")) {
							System.out.println(acc);
						}
						List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=" + GO, "$id");
						boolean containsAnnotation = !query.isEmpty();
						if (filterByAnnotation && !containsAnnotation) {
							continue;
						}
					}
					totalProteins.add(acc);
					SurfacePeptide peptide = new SurfacePeptide(seqs.get(i), ratio);
					surfacePeptides.add(peptide);
					if (surfaceProteinMap.containsKey(acc)) {
						surfaceProteinMap.get(acc).addPeptide(peptide);
					} else {
						SurfaceProtein protein = new SurfaceProtein(acc);
						protein.addPeptide(peptide);
						surfaceProteinMap.put(acc, protein);
					}
				}

			}
			System.out.println(surfaceProteinMap.size() + " proteins with experimental data");
			String aa = "K";
			AtomType atomType = AtomType.NZ;

			boolean printOnlyTheMostAccessibleSite = false;

			int psmsNoAA = 0;
			int psmsWithAANoRatios = 0;
			int discardedContainingMoreThanOneK = 0;
			int validPSMs = 0;

			FileWriter fw = null;

			fw = new FileWriter(outputFile);
			fw.write("\n\nremoveOtherChains=true, removeOtherMolecules=true\n\n");
			removeOtherChains = true;
			removeOtherMolecules = true;
			SurfaceAccessibilityCalculator calc = new SurfaceAccessibilityCalculator(uplr, aa, atomType,
					removeOtherChains, removeOtherMolecules, testAllPosibilities, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa, justTheBestPerPSM);
			System.out.println(calc.getStatistics());
			calc.clearStatistics();
			if (!testAllPosibilities) {

				fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=true\n\n");
				removeOtherChains = false;
				removeOtherMolecules = true;
				calc = new SurfaceAccessibilityCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
						testAllPosibilities, pdbFolder);
				run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa, justTheBestPerPSM);
				System.out.println(calc.getStatistics());
				calc.clearStatistics();
				fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=false\n\n");
				removeOtherChains = false;
				removeOtherMolecules = false;
				calc = new SurfaceAccessibilityCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
						testAllPosibilities, pdbFolder);
				run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa, justTheBestPerPSM);
				System.out.println(calc.getStatistics());
				calc.clearStatistics();
			}
			fw.write("\n\n\n\n\n");
			fw.write(psmIds.size() + " PSMs in total\n");
			fw.write(validPSMs + " PSMs valid (containing quantitative information)\n");
			fw.write(psmsNoAA + " PSMs discarded for not containing " + aa + "\n");
			fw.write(discardedContainingMoreThanOneK + " PSMs discarded having more than one " + aa + "\n");
			fw.write(psmsWithAANoRatios + " PSMs with no experimental quantitative information\n");
			fw.write(peptideSequences.size() + " peptides sequences in total\n");
			fw.write(peptideSequencesValid.size()
					+ " peptides sequences with some site with experimental and theoretical data\n");
			fw.write(uniquePositionsValid.size() + " different sites with experimental and theoretical data\n");
			fw.write(totalProteins.size() + " total different proteins\n");
			fw.write(proteinsWithNoPDBInformation.size() + " proteins with no PDB information\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	private void run(Collection<SurfacePeptide> surfacePeptides, Map<String, SurfaceProtein> surfaceProteinMap,
			SurfaceAccessibilityCalculator calc, FileWriter fw, Set<String> peptideSequences, String aa,
			boolean justTheBestPerPSM) throws IOException {
		printHeader(fw, null);
		SurfaceAccessibilityManager manager = new SurfaceAccessibilityManager(calc);
		manager.setCalculateIfNotPresent(true);
		final Map<String, SurfaceAccessibilityProteinReport> surfaceAccesibilityFromProteins = manager
				.getSurfaceAccesibilityFromProteins(surfaceProteinMap.values());
		System.out.println(
				surfaceAccesibilityFromProteins.size() + " reports for " + surfaceProteinMap.size() + " proteins");
		for (SurfacePeptide surfacePeptide : surfacePeptides) {
			final String peptideSequence = surfacePeptide.getSequence();

			peptideSequences.add(peptideSequence);

			if (!peptideSequence.contains(aa)) {
				continue;
			}
			if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
				// discard if contains more than one K
				continue;
			}

			Set<String> proteinAccs = getSurfaceProteinAccessions(surfacePeptide.getProteins());
			for (String acc : proteinAccs) {
				if (surfaceAccesibilityFromProteins.containsKey(acc)) {
					final SurfaceAccessibilityProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
							.get(acc);
					RatioEx ratio = new RatioEx(surfacePeptide.getRatio(), null, null, "ratio", AggregationLevel.PSM);
					SurfaceAccesibilityReportWriter.printReportForPsm(fw, ratio,
							String.valueOf(surfacePeptide.hashCode()), peptideSequence,
							surfaceAccesibilityProteinReport, aa, false, justTheBestPerPSM);
				}
			}

		}

	}

	private Set<String> getSurfaceProteinAccessions(Set<SurfaceProtein> proteins) {
		Set<String> proteinAccessionsFromAccessions = new HashSet<String>();

		for (SurfaceProtein protein : proteins) {
			proteinAccessionsFromAccessions.add(protein.getAcc());
		}
		return proteinAccessionsFromAccessions;
	}

	private void printHeader(FileWriter fw, List<String> ratioScoreNames) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("ACC").append("\t").append("PSMID").append("\t").append("SEQ").append("\t")
				.append("Position in peptide").append("\t").append(SiteSurfaceAccessibilityReport.getToStringHeaders())
				.append("\t").append("AREA_RATIO").append("\t");
		if (ratioScoreNames != null && !ratioScoreNames.isEmpty()) {
			for (String ratioScoreName : ratioScoreNames) {
				sb.append(ratioScoreName).append("\t");
			}
		}

		sb.append("RATIO").append("\t");
		if (ratioScoreNames != null) {
			for (String ratioScoreName : ratioScoreNames) {
				sb.append(ratioScoreName).append("\t");
			}
		}
		sb.append("\n");
		fw.write(sb.toString());
	}

	@Test
	public void statisticsOnRemovingOrNotElementsInPDB() {
		boolean filterByRemoval = false;
		boolean removeOtherChains3 = true;
		boolean removeOtherMolecules3 = true;
		String GO = "GO:0003677";
		boolean filterByAnnotation = false;
		Double minRatio = -Double.MAX_VALUE;
		// Double minRatio = 0.0;
		// this gets the experimental data from a table
		OntologyTermI goTerm = new GORetriever(new File("z:\\share\\Salva\\data\\go")).getGOTermByID(GO);
		String goName = "";
		if (goTerm != null) {
			goName = goTerm.getPreferredName();
		}
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		final String ACC = "ACC";
		final String POSITION_IN_UNIPROT = "position_in_uniprot";
		final String PSMID = "PSMID";
		final String PDB_ID = "pdb_ID";
		final String SURFACE_ACCESSIBILITY = "Surface_accessibility";
		final String CHAIN_ID = "chain_id";
		final String OTHER_CHAINS_REMOVED = "OtherChainsRemoved";
		final String OTHER_MOLECULES_REMOVED = "OtherMoleculesRemoved";
		final String AREA_RATIO = "AREA_RATIO";
		final String ATOM_NUMBER = "atom_number";
		File inputFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_DNA binding_ratio_gte=0.0.xlsx");
		ExcelReader reader;
		try {
			reader = new ExcelReader(inputFile, 0, 0);
			List<String> columnNames = reader.getColumnNames().get(0);
			Map<String, Integer> indexByCol = new HashMap<String, Integer>();
			int index = 0;
			for (String columnName : columnNames) {
				indexByCol.put(columnName, index++);
			}
			int rowIndex = 1;
			Row row = null;

			Map<String, List<Surface>> surfacesByPSM = new HashMap<String, List<Surface>>();
			Map<String, List<Surface>> surfacesByUniprotPosition = new HashMap<String, List<Surface>>();
			Map<String, List<Surface>> surfacesByPDBKey = new HashMap<String, List<Surface>>();
			while ((row = reader.getWorkbook().getSheetAt(0).getRow(rowIndex)) != null) {
				try {
					String acc = reader.getStringValue(0, rowIndex, indexByCol.get(ACC));
					if (filterByAnnotation) {
						Entry entry = uplr.getAnnotatedProtein(null, acc).get(acc);
						if (entry != null) {
							List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=" + GO, "$id");
							boolean containsAnnotation = !query.isEmpty();
							if (!containsAnnotation) {
								continue;
							}
						}
					}
					Double positionInUniprot = Double
							.valueOf(reader.getStringValue(0, rowIndex, indexByCol.get(POSITION_IN_UNIPROT)));
					String positionInUniprotKey = acc + "-" + positionInUniprot;
					String psmID = reader.getStringValue(0, rowIndex, indexByCol.get(PSMID));
					String pdbID = reader.getStringValue(0, rowIndex, indexByCol.get(PDB_ID));
					String chainID = reader.getStringValue(0, rowIndex, indexByCol.get(CHAIN_ID));
					int atomNumber = Double.valueOf(reader.getStringValue(0, rowIndex, indexByCol.get(ATOM_NUMBER)))
							.intValue();
					boolean removeOtherChains2 = Boolean
							.valueOf(reader.getStringValue(0, rowIndex, indexByCol.get(OTHER_CHAINS_REMOVED)));
					boolean removeOtherMolecules2 = Boolean
							.valueOf(reader.getStringValue(0, rowIndex, indexByCol.get(OTHER_MOLECULES_REMOVED)));
					if (filterByRemoval) {
						if (Boolean.compare(removeOtherChains2, removeOtherChains3) != 0) {
							continue;
						}
						if (Boolean.compare(removeOtherMolecules2, removeOtherMolecules3) != 0) {
							continue;
						}
					}
					double accessibility = Double
							.valueOf(reader.getNumberValue(0, rowIndex, indexByCol.get(SURFACE_ACCESSIBILITY)));
					double experimentalRatio = Double
							.valueOf(reader.getNumberValue(0, rowIndex, indexByCol.get(AREA_RATIO)));
					Surface surface = new Surface(pdbID, chainID, removeOtherChains2, removeOtherMolecules2,
							accessibility, experimentalRatio, psmID, atomNumber);
					String pdbKey = surface.getKey();
					if (surfacesByPDBKey.containsKey(pdbKey)) {
						surfacesByPDBKey.get(pdbKey).add(surface);
					} else {
						List<Surface> list = new ArrayList<Surface>();
						list.add(surface);
						surfacesByPDBKey.put(pdbKey, list);
					}
					if (surfacesByPSM.containsKey(psmID)) {
						surfacesByPSM.get(psmID).add(surface);
					} else {
						List<Surface> list = new ArrayList<Surface>();
						list.add(surface);
						surfacesByPSM.put(psmID, list);
					}
					if (surfacesByUniprotPosition.containsKey(positionInUniprotKey)) {
						surfacesByUniprotPosition.get(positionInUniprotKey).add(surface);
					} else {
						List<Surface> list = new ArrayList<Surface>();
						list.add(surface);
						surfacesByUniprotPosition.put(positionInUniprotKey, list);
					}
				} finally {
					rowIndex++;
				}
			}
			System.out.println("Remove other chains = " + removeOtherChains3);
			System.out.println("Remove other molecules = " + removeOtherMolecules3);
			System.out.println("By PSM:");
			printStatistics(surfacesByPSM);
			System.out.println("\nBy Position in protein:");
			printStatistics(surfacesByUniprotPosition);
			System.out.println("\nBy PDB key:");
			printStatistics(surfacesByPDBKey);
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}

	}

	private void printStatistics(Map<String, List<Surface>> map) {
		String pdbModelWithBiggerDifference = "";
		List<Double> surfacesSuperSTD = new ArrayList<Double>();
		List<Double> ratiosSuperSTD = new ArrayList<Double>();
		List<Integer> numberOfSurfacesPerGroup = new ArrayList<Integer>();
		Set<String> psmIDs = new HashSet<String>();
		List<Double> differences = new ArrayList<Double>();
		double superMaxDiff = 0;
		for (String key : map.keySet()) {
			double maxDiff = 0;
			Set<String> surfaceKeys = new HashSet<String>();
			List<Double> surfaces = new ArrayList<Double>();
			List<Double> ratios = new ArrayList<Double>();
			List<Surface> list = map.get(key);
			for (Surface surface : list) {
				surfaces.add(surface.getAccessibility());
				if (!psmIDs.contains(surface.getPSMID())) {
					ratios.add(surface.getExperimentalAccessibility());
				}
				String surfaceKey = surface.getKey();
				if (!surfaceKeys.contains(surfaceKey)) {
					surfaceKeys.add(surfaceKey);
				}
			}
			numberOfSurfacesPerGroup.add(surfaceKeys.size());
			Double[] surfacesArray = surfaces.toArray(new Double[0]);
			Double[] ratiosArray = ratios.toArray(new Double[0]);
			double surfacesStdev = Maths.stddev(surfacesArray);
			double ratiosStdev = Maths.stddev(ratiosArray);
			surfacesSuperSTD.add(surfacesStdev);
			ratiosSuperSTD.add(ratiosStdev);
			// difference
			if (surfaces.size() > 1) {
				double difference = Math.abs(surfaces.get(0) - surfaces.get(1));
				if (maxDiff < difference) {
					maxDiff = difference;

				}
				if (superMaxDiff < difference) {
					superMaxDiff = difference;
					pdbModelWithBiggerDifference = list.get(0).getPdbID();
				}
			}
			if (surfaces.size() > 2) {
				double difference = Math.abs(surfaces.get(0) - surfaces.get(2));
				if (maxDiff < difference) {
					maxDiff = difference;
				}
				if (superMaxDiff < difference) {
					superMaxDiff = difference;
					pdbModelWithBiggerDifference = list.get(0).getPdbID();
				}
				difference = Math.abs(surfaces.get(1) - surfaces.get(2));
				if (maxDiff < difference) {
					maxDiff = difference;
				}
				if (superMaxDiff < difference) {
					superMaxDiff = difference;
					pdbModelWithBiggerDifference = list.get(0).getPdbID();
				}
			}
			differences.add(maxDiff);
		}
		Integer[] numberOfSurfacesPerGroupArray = numberOfSurfacesPerGroup.toArray(new Integer[0]);
		Double[] surfacesSuperSTDArray = surfacesSuperSTD.toArray(new Double[0]);
		Double[] ratiosSuperSTDArray = ratiosSuperSTD.toArray(new Double[0]);
		Double[] differencesArray = differences.toArray(new Double[0]);
		System.out.println(map.size() + " groups");
		System.out.println("Mean number of models per group " + Maths.mean(numberOfSurfacesPerGroupArray));
		System.out
				.println("Mean standard deviation of model surface accessibity = " + Maths.mean(surfacesSuperSTDArray));
		System.out.println("Mean standard deviation of ratios = " + Maths.mean(ratiosSuperSTDArray));
		System.out.println("Maximum difference between surfaces in group = " + Maths.max(differencesArray));
		System.out.println("Average difference between surfaces in group = " + Maths.mean(differencesArray));
		System.out.println("PDB model with more differences = " + pdbModelWithBiggerDifference);
	}

}
