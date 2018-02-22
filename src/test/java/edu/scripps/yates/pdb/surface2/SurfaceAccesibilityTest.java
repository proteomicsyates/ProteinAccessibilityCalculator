package edu.scripps.yates.pdb.surface2;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

import edu.scripps.yates.annotations.go.GORetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinRetrievalSettings;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.dbindex.util.DigestionConfiguration;
import edu.scripps.yates.excel.proteindb.importcfg.adapter.ImportCfgFileReader;
import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.ProteinReportWriter;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.surface.SurfaceCalculator;
import edu.scripps.yates.pdb.surface.SurfaceProteinReport;
import edu.scripps.yates.pdb.surface.SurfaceProteinReportManager;
import edu.scripps.yates.pdb.surface.SurfaceReport;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.jaxb.xpathquery.JAXBXPathQuery;
import edu.scripps.yates.utilities.model.enums.AggregationLevel;
import edu.scripps.yates.utilities.model.factories.RatioEx;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Project;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;
import psidev.psi.tools.ontology_manager.interfaces.OntologyTermI;

public class SurfaceAccesibilityTest {

	private boolean removeOtherChains = true;
	private boolean removeOtherMolecules = true;

	@Test
	public void surfaceAccebilityTest1() {
		final Protein protein = new Protein("Q9UKV8");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		char aa = 'K';
		AtomType atomType = AtomType.NZ;

		SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
				true, pdbFolder);
		final SurfaceProteinReport proteinReport = calc.getReportFromProtein(protein, null, null);
		final TIntObjectHashMap<Set<SurfaceReport>> accesibilitiesByPositionInUniprotSeq = proteinReport
				.getReportsByPositionInUniprotSeq();
		List<Integer> positionList = new ArrayList<Integer>();
		for (int position : accesibilitiesByPositionInUniprotSeq.keys()) {
			positionList.add(position);
		}
		Collections.sort(positionList);
		for (Integer position : positionList) {
			final Set<SurfaceReport> surfaceAccesibilityReports = accesibilitiesByPositionInUniprotSeq.get(position);
			for (SurfaceReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
				System.out.println(surfaceAccessibilityReport);
			}

		}

	}

	@Test
	public void surfaceAccebilityTestWithManager() {

		// accs.add("Q9HD26-2");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		char aa = 'K';
		AtomType atomType = AtomType.NZ;

		SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
				true, pdbFolder);
		calc.setDigestionConfiguration(getDigestionConfiguration());
		SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
		final SurfaceProteinReport surfaceAccesibilityProteinReport = manager.getProteinReportByProtein("Q9UKV8", null);

		final TIntObjectHashMap<Set<SurfaceReport>> accesibilitiesByPositionInUniprotSeq = surfaceAccesibilityProteinReport
				.getReportsByPositionInUniprotSeq();
		List<Integer> positionList = new ArrayList<Integer>();
		for (int position : accesibilitiesByPositionInUniprotSeq.keys()) {
			positionList.add(position);
		}

		Collections.sort(positionList);
		for (Integer position : positionList) {
			final Set<SurfaceReport> surfaceAccesibilityReports = accesibilitiesByPositionInUniprotSeq.get(position);
			for (SurfaceReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
				System.out.println(surfaceAccessibilityReport);
			}

		}

	}

	private DigestionConfiguration getDigestionConfiguration() {
		char[] enzymeArray = { 'K' };
		int numMisscleavages = 1;
		boolean semiCleavage = false;
		String peptideFilterString = null;
		DigestionConfiguration ret = new DigestionConfiguration(enzymeArray, numMisscleavages, semiCleavage,
				peptideFilterString);
		return ret;
	}

	@Test
	public void surfaceAccesibilityFromRemoteFiles() {

		File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		File projectConfigFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		char aa = 'K';
		AtomType atomType = AtomType.NZ;
		boolean printOnlyTheMostAccessibleSite = false;

		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
					true, pdbFolder);
			Set<Protein> proteinAccs = getProteinAccsFromProject(project);

			List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			manager.setCalculateIfNotPresent(true);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(proteinAccs);
			Set<String> psmIds = new THashSet<String>();
			Set<String> peptideSequences = new THashSet<String>();
			Set<String> peptideSequencesValid = new THashSet<String>();
			Set<String> uniquePositionsValid = new THashSet<String>();
			Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (Condition condition : conditions) {
				final Set<PSM> psMs = condition.getPSMs();
				for (PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getPSMIdentifier())) {
						psmIds.add(psm.getPSMIdentifier());
						if (!peptideSequence.contains(String.valueOf(aa))) {
							psmsNoAA++;
							continue;
						}
						if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
							// discard if contains more than one K
							discardedContainingMoreThanOneK++;
							continue;
						}
						Ratio ratio = null;
						for (Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getPSMIdentifier(),
										peptideSequence, surfaceAccesibilityProteinReport, aa,
										printOnlyTheMostAccessibleSite);
								uniquePositionsValid
										.addAll(surfaceAccesibilityProteinReport.getUniquePositionsInProteinKeys());
							} else {
								proteinsWithNoPDBInformation.add(acc);
							}
						}
					}
				}
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

	@Test
	public void surfaceAccesibilityFromRemoteFiles_DNABinding() {

		File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_DNABinding.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		File projectConfigFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		char aa = 'K';
		AtomType atomType = AtomType.NZ;
		boolean printOnlyTheMostAccessibleSite = false;

		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
					true, pdbFolder);
			Set<Protein> proteinAccs = getProteinAccsFromProject(project);

			List<Protein> dnaBindingProteins = filterDNABindingProteins(proteinAccs, uplr);

			List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			manager.setCalculateIfNotPresent(true);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(dnaBindingProteins, null);
			Set<String> psmIds = new THashSet<String>();
			Set<String> peptideSequences = new THashSet<String>();
			Set<String> peptideSequencesValid = new THashSet<String>();
			Set<String> uniquePositionsValid = new THashSet<String>();
			Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (Condition condition : conditions) {
				final Set<PSM> psMs = condition.getPSMs();
				for (PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getPSMIdentifier())) {
						psmIds.add(psm.getPSMIdentifier());
						if (!peptideSequence.contains(String.valueOf(aa))) {
							psmsNoAA++;
							continue;
						}
						if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
							// discard if contains more than one K
							discardedContainingMoreThanOneK++;
							continue;
						}
						Ratio ratio = null;
						for (Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getPSMIdentifier(),
										peptideSequence, surfaceAccesibilityProteinReport, aa,
										printOnlyTheMostAccessibleSite);
								uniquePositionsValid
										.addAll(surfaceAccesibilityProteinReport.getUniquePositionsInProteinKeys());
							} else {
								proteinsWithNoPDBInformation.add(acc);
							}
						}
					}
				}
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

	@Test
	public void surfaceAccesibilityFromRemoteFiles_DNABinding_FromTable() {
		String GO = "GO:0003677";
		boolean filterByAnnotation = false;
		// Double minRatio = -Double.MAX_VALUE;
		Double minRatio = 0.0;
		// this gets the experimental data from a table
		OntologyTermI goTerm = new GORetriever(new File("z:\\share\\Salva\\data\\go")).getGOTermByID(GO);
		String goName = "";
		if (goTerm != null) {
			goName = goTerm.getPreferredName();
		}
		File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_" + goName
				+ "_ratio_gte=" + minRatio + ".csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		File experimentalDataFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\experimental_proteins.csv");
		try {
			Set<String> psmIds = new THashSet<String>();
			Set<String> peptideSequences = new THashSet<String>();
			Set<String> peptideSequencesValid = new THashSet<String>();
			Set<String> uniquePositionsValid = new THashSet<String>();
			Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			Set<String> totalProteins = new THashSet<String>();
			Set<Peptide> surfacePeptides = new THashSet<Peptide>();
			List<String> accs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 0, true);
			List<String> seqs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 1, true);
			List<String> ratios = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 2, true);
			Map<String, Protein> surfaceProteinMap = new THashMap<String, Protein>();
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
					Peptide peptide = new Peptide(seqs.get(i), ratio);
					surfacePeptides.add(peptide);
					if (surfaceProteinMap.containsKey(acc)) {
						surfaceProteinMap.get(acc).addPeptide(peptide);
					} else {
						Protein protein = new Protein(acc);
						protein.addPeptide(peptide);
						surfaceProteinMap.put(acc, protein);
					}
				}

			}
			System.out.println(surfaceProteinMap.size() + " proteins with experimental data");
			char aa = 'K';
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
			SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
					true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(calc.getStatistics());
			calc.clearStatistics();
			fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=true\n\n");
			removeOtherChains = false;
			removeOtherMolecules = true;
			calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(calc.getStatistics());
			calc.clearStatistics();
			fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=false\n\n");
			removeOtherChains = false;
			removeOtherMolecules = false;
			calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(calc.getStatistics());
			calc.clearStatistics();
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

	@Test
	public void filterProteinsByAnnotations() {
		File experimentalDataFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\accessible_experimental_proteins.csv");
		try {
			List<String> accs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 0, true);

			final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
			UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
			UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
			Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			for (String acc : accs) {
				Entry entry = annotatedProteins.get(acc);
				if (entry != null) {
					List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=GO:0003677", "$id");
					boolean containsAnnotation = !query.isEmpty();
					System.out.println(acc + "\t" + containsAnnotation);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void run(Collection<Peptide> surfacePeptides, Map<String, Protein> surfaceProteinMap,
			SurfaceCalculator calc, FileWriter fw, Set<String> peptideSequences, char aa) throws IOException {
		printHeader(fw, null);
		SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
		manager.setCalculateIfNotPresent(true);
		final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
				.getReportsFromProteins(surfaceProteinMap.values());
		System.out.println(
				surfaceAccesibilityFromProteins.size() + " reports for " + surfaceProteinMap.size() + " proteins");
		for (Peptide surfacePeptide : surfacePeptides) {
			final String peptideSequence = surfacePeptide.getSequence();

			peptideSequences.add(peptideSequence);

			if (!peptideSequence.contains(String.valueOf(aa))) {
				continue;
			}
			if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
				// discard if contains more than one K
				continue;
			}

			Set<String> proteinAccs = getSurfaceProteinAccessions(surfacePeptide.getProteins());
			for (String acc : proteinAccs) {
				if (surfaceAccesibilityFromProteins.containsKey(acc)) {
					final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
							.get(acc);
					RatioEx ratio = new RatioEx(surfacePeptide.getRatio(), null, null, "ratio", AggregationLevel.PSM);
					ProteinReportWriter.printReportForPsm(fw, ratio, String.valueOf(surfacePeptide.hashCode()),
							peptideSequence, surfaceAccesibilityProteinReport, aa, false);
				}
			}

		}

	}

	private List<Protein> filterDNABindingProteins(Set<Protein> proteinAccs, UniprotProteinLocalRetriever uplr) {

		List<Protein> ret = new ArrayList<Protein>();
		File dnaBindingProteins = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\DNABindingProteins.csv");
		try {
			FileWriter fw = null;

			Set<String> accs = new THashSet<String>();
			if (dnaBindingProteins.exists()) {
				accs.addAll(Files.lines(Paths.get(dnaBindingProteins.toURI())).collect(Collectors.toSet()));

			}

			fw = new FileWriter(dnaBindingProteins, true);

			ProgressCounter counter = new ProgressCounter(proteinAccs.size(), ProgressPrintingType.PERCENTAGE_STEPS, 1);
			for (Protein surfaceProtein : proteinAccs) {
				counter.increment();
				String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					System.out.println(printIfNecessary);
				}
				boolean accessibleRatio = false;
				for (Peptide peptide : surfaceProtein.getPeptides()) {
					if (peptide.getRatio() > 0) {
						accessibleRatio = true;
					}
				}
				if (!accessibleRatio) {
					continue;
				}

				if (accs.contains(surfaceProtein.getAcc())) {
					ret.add(surfaceProtein);
				} else {
					Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, surfaceProtein.getAcc());
					if (annotatedProtein.containsKey(surfaceProtein.getAcc())) {
						Entry entry = annotatedProtein.get(surfaceProtein.getAcc());
						List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=GO:0003677", "$id");
						if (!query.isEmpty()) {
							if (fw != null) {
								fw.write(surfaceProtein.getAcc() + "\n");
								fw.flush();
							}
							ret.add(surfaceProtein);
						}
					}
				}
			}
			System.out.println(ret.size() + " proteins that are annotated as DNA binding");

			if (fw != null) {
				fw.close();
			}
		} catch (IOException e) {

		}
		return ret;
	}

	@Test
	public void surfaceAccesibilityFromRemoteFilesMaxPerSite() {

		File outputFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\C_max_per_site_surface.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		File projectConfigFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		char aa = 'K';
		AtomType atomType = AtomType.NZ;
		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
					true, pdbFolder);
			Set<Protein> proteinAccs = getProteinAccsFromProject(project);
			List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			// everything is already calculated, so not go for it if not present
			manager.setCalculateIfNotPresent(false);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(proteinAccs);
			Set<String> psmIds = new THashSet<String>();
			Set<String> peptideSequences = new THashSet<String>();
			Set<String> peptideSequencesValid = new THashSet<String>();
			Set<String> uniquePositionsValid = new THashSet<String>();
			Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (Condition condition : conditions) {
				if (condition.getName().startsWith("HS")) {
					continue;
				}
				final Set<PSM> psMs = condition.getPSMs();
				for (PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getPSMIdentifier())) {
						psmIds.add(psm.getPSMIdentifier());
						if (!peptideSequence.contains(String.valueOf(aa))) {
							psmsNoAA++;
							continue;
						}
						if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
							// discard if contains more than one K
							discardedContainingMoreThanOneK++;
							continue;
						}
						Ratio ratio = null;
						for (Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getPSMIdentifier(),
										peptideSequence, surfaceAccesibilityProteinReport, aa, true);
								final List<String> uniquePositionsInProteinKeys = surfaceAccesibilityProteinReport
										.getUniquePositionsInProteinKeys();
								uniquePositionsValid.addAll(uniquePositionsInProteinKeys);
							} else {
								proteinsWithNoPDBInformation.add(acc);
							}
						}
					}
				}
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

	private void printHeader(FileWriter fw, List<String> ratioScoreNames) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("ACC").append("\t").append("PSMID").append("\t").append("SEQ").append("\t")
				.append("Position in peptide").append("\t").append(JMolAtomReport.getStaticHeaders()).append("\t")
				.append("AREA_RATIO").append("\t");
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

	private Set<String> getProteinAccessions(Set<edu.scripps.yates.utilities.proteomicsmodel.Protein> proteins) {
		Set<String> proteinAccessionsFromAccessions = new THashSet<String>();

		for (edu.scripps.yates.utilities.proteomicsmodel.Protein protein : proteins) {
			proteinAccessionsFromAccessions.add(protein.getAccession());
		}
		return proteinAccessionsFromAccessions;
	}

	private Set<String> getSurfaceProteinAccessions(Set<Protein> proteins) {
		Set<String> proteinAccessionsFromAccessions = new THashSet<String>();

		for (Protein protein : proteins) {
			proteinAccessionsFromAccessions.add(protein.getAcc());
		}
		return proteinAccessionsFromAccessions;
	}

	private Set<Protein> getProteinAccsFromProject(Project project) {
		Set<Protein> proteinAccessionsFromProject = new THashSet<Protein>();
		final Set<Condition> conditions = project.getConditions();
		Map<String, Peptide> surfacePeptideMap = new THashMap<String, Peptide>();
		for (Condition condition : conditions) {
			final Set<edu.scripps.yates.utilities.proteomicsmodel.Protein> proteins = condition.getProteins();
			for (edu.scripps.yates.utilities.proteomicsmodel.Protein protein : proteins) {
				Protein surfaceProtein = new Protein(FastaParser.getUniProtACC(protein.getAccession()));
				for (PSM psm : protein.getPSMs()) {
					if (psm.getRatios() != null && !psm.getRatios().isEmpty()) {
						Peptide surfacePeptide = null;
						if (surfacePeptideMap.containsKey(psm.getPSMIdentifier())) {
							surfacePeptide = surfacePeptideMap.get(psm.getPSMIdentifier());
						} else {
							surfacePeptide = new Peptide(psm.getSequence(),
									psm.getRatios().iterator().next().getValue());
						}
						surfaceProtein.getPeptides().add(surfacePeptide);
					}
				}
				proteinAccessionsFromProject.add(surfaceProtein);
			}
		}
		return proteinAccessionsFromProject;
	}

	private List<String> getRatioScoreNamesFromProject(Project project) {
		Set<String> ret = new THashSet<String>();
		final Set<Condition> conditions = project.getConditions();
		for (Condition condition : conditions) {
			final Set<PSM> psms = condition.getPSMs();
			for (PSM psm : psms) {
				final Set<Ratio> ratios = psm.getRatios();
				if (ratios != null) {
					for (Ratio ratio : ratios) {
						if (ratio.getAssociatedConfidenceScore() != null) {
							ret.add(ratio.getAssociatedConfidenceScore().getScoreName());
						}
					}
				}
			}
		}
		List<String> ratioScoreNames = new ArrayList<String>();
		ratioScoreNames.addAll(ret);
		Collections.sort(ratioScoreNames);
		return ratioScoreNames;
	}

	@Test
	public void getSurfaceAccebilityFromPDBSites() {
		UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		char aa = 'K';
		AtomType atomType = AtomType.NZ;
		String protein = "5A5T";
		String chain = null;
		SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
				true, pdbFolder);
		final SurfaceProteinReport surfaceAccesibilityProteinReport = calc.getReportFromPDBModel(protein, chain);
		final List<SurfaceReport> reports = surfaceAccesibilityProteinReport.getReports();
		boolean first = true;
		for (SurfaceReport report : reports) {
			if (first) {
				System.out.println(report.getToStringHeaders());
				first = false;
			}
			System.out.println(report);
		}

	}
}
