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
import java.util.stream.Stream;

import org.junit.Test;

import edu.scripps.yates.annotations.go.GORetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinRetrievalSettings;
import edu.scripps.yates.dbindex.util.DigestionConfiguration;
import edu.scripps.yates.excel.proteindb.importcfg.adapter.ImportCfgFileReader;
import edu.scripps.yates.pdb.Calculator;
import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.ProteinReportWriter;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.surface.SurfaceCalculator;
import edu.scripps.yates.pdb.surface.SurfaceProteinReport;
import edu.scripps.yates.pdb.surface.SurfaceProteinReportManager;
import edu.scripps.yates.pdb.surface.SurfaceReport;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.jaxb.xpathquery.JAXBXPathQuery;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Project;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AggregationLevel;
import edu.scripps.yates.utilities.proteomicsmodel.factories.RatioEx;
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
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;

		final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
				removeOtherMolecules, true, pdbFolder);
		final SurfaceProteinReport proteinReport = calc.getReportFromProtein(protein, null, null);
		final TIntObjectHashMap<Set<SurfaceReport>> accesibilitiesByPositionInUniprotSeq = proteinReport
				.getReportsByPositionInUniprotSeq();
		final List<Integer> positionList = new ArrayList<Integer>();
		for (final int position : accesibilitiesByPositionInUniprotSeq.keys()) {
			positionList.add(position);
		}
		Collections.sort(positionList);
		for (final Integer position : positionList) {
			final Set<SurfaceReport> surfaceAccesibilityReports = accesibilitiesByPositionInUniprotSeq.get(position);
			for (final SurfaceReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
				System.out.println(surfaceAccessibilityReport);
			}

		}

	}

	@Test
	public void surfaceAccebilityTestWithManager() {

		// accs.add("Q9HD26-2");
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;

		final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
				removeOtherMolecules, true, pdbFolder);
		calc.setDigestionConfiguration(getDigestionConfiguration());
		final SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
		final SurfaceProteinReport surfaceAccesibilityProteinReport = manager.getProteinReportByProtein("Q9UKV8", null);

		final TIntObjectHashMap<Set<SurfaceReport>> accesibilitiesByPositionInUniprotSeq = surfaceAccesibilityProteinReport
				.getReportsByPositionInUniprotSeq();
		final List<Integer> positionList = new ArrayList<Integer>();
		for (final int position : accesibilitiesByPositionInUniprotSeq.keys()) {
			positionList.add(position);
		}

		Collections.sort(positionList);
		for (final Integer position : positionList) {
			final Set<SurfaceReport> surfaceAccesibilityReports = accesibilitiesByPositionInUniprotSeq.get(position);
			for (final SurfaceReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
				System.out.println(surfaceAccessibilityReport);
			}

		}

	}

	private DigestionConfiguration getDigestionConfiguration() {
		final char[] enzymeArray = { 'K' };
		final int numMisscleavages = 1;
		final boolean semiCleavage = false;
		final String peptideFilterString = null;
		final DigestionConfiguration ret = new DigestionConfiguration(enzymeArray, numMisscleavages, semiCleavage,
				peptideFilterString);
		return ret;
	}

	@Test
	public void surfaceAccesibilityFromRemoteFiles() {

		final File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		final File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		final File projectConfigFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		final ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;
		final boolean printOnlyTheMostAccessibleSite = false;

		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
					removeOtherMolecules, true, pdbFolder);
			final Set<Protein> proteinAccs = getProteinAccsFromProject(project);

			final List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			final SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			manager.setCalculateIfNotPresent(true);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(proteinAccs);
			final Set<String> psmIds = new THashSet<String>();
			final Set<String> peptideSequences = new THashSet<String>();
			final Set<String> peptideSequencesValid = new THashSet<String>();
			final Set<String> uniquePositionsValid = new THashSet<String>();
			final Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			final Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (final Condition condition : conditions) {
				final Set<PSM> psMs = condition.getPSMs();
				for (final PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getIdentifier())) {
						psmIds.add(psm.getIdentifier());
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
						for (final Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						final Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (final String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getIdentifier(), peptideSequence,
										surfaceAccesibilityProteinReport, aa, printOnlyTheMostAccessibleSite);
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
		} catch (final IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	@Test
	public void surfaceAccesibilityFromRemoteFiles_DNABinding() {

		final File outputFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_DNABinding.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		final File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		final File projectConfigFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		final ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;
		final boolean printOnlyTheMostAccessibleSite = false;

		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
					removeOtherMolecules, true, pdbFolder);
			final Set<Protein> proteinAccs = getProteinAccsFromProject(project);

			final List<Protein> dnaBindingProteins = filterDNABindingProteins(proteinAccs, uplr);

			final List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			final SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			manager.setCalculateIfNotPresent(true);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(dnaBindingProteins, null);
			final Set<String> psmIds = new THashSet<String>();
			final Set<String> peptideSequences = new THashSet<String>();
			final Set<String> peptideSequencesValid = new THashSet<String>();
			final Set<String> uniquePositionsValid = new THashSet<String>();
			final Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			final Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (final Condition condition : conditions) {
				final Set<PSM> psMs = condition.getPSMs();
				for (final PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getIdentifier())) {
						psmIds.add(psm.getIdentifier());
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
						for (final Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						final Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (final String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getIdentifier(), peptideSequence,
										surfaceAccesibilityProteinReport, aa, printOnlyTheMostAccessibleSite);
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
		} catch (final IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	@Test
	public void surfaceAccesibilityFromRemoteFiles_DNABinding_FromTable() {
		final String GO = "GO:0003677";
		final boolean filterByAnnotation = false;
		// Double minRatio = -Double.MAX_VALUE;
		final Double minRatio = 0.0;
		// this gets the experimental data from a table
		final OntologyTermI goTerm = new GORetriever(new File("z:\\share\\Salva\\data\\go")).getGOTermByID(GO);
		String goName = "";
		if (goTerm != null) {
			goName = goTerm.getPreferredName();
		}
		final File outputFile = new File("z:\\share\\Salva\\data\\PINT projects\\molecular painting\\surface_" + goName
				+ "_ratio_gte=" + minRatio + ".csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		final File experimentalDataFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\experimental_proteins.csv");
		try {
			final Set<String> psmIds = new THashSet<String>();
			final Set<String> peptideSequences = new THashSet<String>();
			final Set<String> peptideSequencesValid = new THashSet<String>();
			final Set<String> uniquePositionsValid = new THashSet<String>();
			final Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			final Set<String> totalProteins = new THashSet<String>();
			final Set<Peptide> surfacePeptides = new THashSet<Peptide>();
			final List<String> accs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 0, true);
			final List<String> seqs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 1, true);
			final List<String> ratios = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 2, true);
			final Map<String, Protein> surfaceProteinMap = new THashMap<String, Protein>();
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);

			for (int i = 0; i < accs.size(); i++) {

				final Double ratio = Double.valueOf(ratios.get(i));
				if (ratio >= minRatio) {
					final String acc = accs.get(i);
					if (acc.equals("P09211")) {
						System.out.println(acc);
					}
					final Entry entry = annotatedProteins.get(acc);
					if (entry != null) {
						if (acc.equals("P09211")) {
							System.out.println(acc);
						}
						final List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=" + GO, "$id");
						final boolean containsAnnotation = !query.isEmpty();
						if (filterByAnnotation && !containsAnnotation) {
							continue;
						}
					}
					totalProteins.add(acc);
					final Peptide peptide = new Peptide(seqs.get(i), ratio);
					surfacePeptides.add(peptide);
					if (surfaceProteinMap.containsKey(acc)) {
						surfaceProteinMap.get(acc).addPeptide(peptide);
					} else {
						final Protein protein = new Protein(acc);
						protein.addPeptide(peptide);
						surfaceProteinMap.put(acc, protein);
					}
				}

			}
			System.out.println(surfaceProteinMap.size() + " proteins with experimental data");
			final char aa = 'K';
			final AtomType atomType = AtomType.NZ;

			final boolean printOnlyTheMostAccessibleSite = false;

			final int psmsNoAA = 0;
			final int psmsWithAANoRatios = 0;
			final int discardedContainingMoreThanOneK = 0;
			final int validPSMs = 0;

			FileWriter fw = null;

			fw = new FileWriter(outputFile);
			fw.write("\n\nremoveOtherChains=true, removeOtherMolecules=true\n\n");
			removeOtherChains = true;
			removeOtherMolecules = true;
			SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules,
					true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(Calculator.getStatistics());
			Calculator.clearStatistics();
			fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=true\n\n");
			removeOtherChains = false;
			removeOtherMolecules = true;
			calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(Calculator.getStatistics());
			Calculator.clearStatistics();
			fw.write("\n\nremoveOtherChains=false, removeOtherMolecules=false\n\n");
			removeOtherChains = false;
			removeOtherMolecules = false;
			calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, true, pdbFolder);
			run(surfacePeptides, surfaceProteinMap, calc, fw, peptideSequences, aa);
			System.out.println(Calculator.getStatistics());
			Calculator.clearStatistics();
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
		} catch (final IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	@Test
	public void filterProteinsByAnnotations() {
		final File experimentalDataFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\accessible_experimental_proteins.csv");
		try {
			final List<String> accs = FileUtils.readColumnFromTextFile(experimentalDataFile, ",", 0, true);

			final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
			UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			for (final String acc : accs) {
				final Entry entry = annotatedProteins.get(acc);
				if (entry != null) {
					final List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=GO:0003677", "$id");
					final boolean containsAnnotation = !query.isEmpty();
					System.out.println(acc + "\t" + containsAnnotation);
				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private void run(Collection<Peptide> surfacePeptides, Map<String, Protein> surfaceProteinMap,
			SurfaceCalculator calc, FileWriter fw, Set<String> peptideSequences, char aa) throws IOException {
		printHeader(fw, null);
		final SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
		manager.setCalculateIfNotPresent(true);
		final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
				.getReportsFromProteins(surfaceProteinMap.values());
		System.out.println(
				surfaceAccesibilityFromProteins.size() + " reports for " + surfaceProteinMap.size() + " proteins");
		for (final Peptide surfacePeptide : surfacePeptides) {
			final String peptideSequence = surfacePeptide.getSequence();

			peptideSequences.add(peptideSequence);

			if (!peptideSequence.contains(String.valueOf(aa))) {
				continue;
			}
			if (StringUtils.allPositionsOf(peptideSequence, aa).size() > 1) {
				// discard if contains more than one K
				continue;
			}

			final Set<String> proteinAccs = getSurfaceProteinAccessions(surfacePeptide.getProteins());
			for (final String acc : proteinAccs) {
				if (surfaceAccesibilityFromProteins.containsKey(acc)) {
					final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
							.get(acc);
					final Ratio ratio = new RatioEx(surfacePeptide.getRatio(), null, null, "ratio",
							AggregationLevel.PSM);
					ProteinReportWriter.printReportForPsm(fw, ratio, String.valueOf(surfacePeptide.hashCode()),
							peptideSequence, surfaceAccesibilityProteinReport, aa, false);
				}
			}

		}

	}

	private List<Protein> filterDNABindingProteins(Set<Protein> proteinAccs, UniprotProteinLocalRetriever uplr) {

		final List<Protein> ret = new ArrayList<Protein>();
		final File dnaBindingProteins = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\DNABindingProteins.csv");
		try {
			FileWriter fw = null;

			final Set<String> accs = new THashSet<String>();
			if (dnaBindingProteins.exists()) {
				Stream<String> stream = null;
				try {
					stream = Files.lines(Paths.get(dnaBindingProteins.toURI()));
					accs.addAll(stream.collect(Collectors.toSet()));
				} finally {
					if (stream != null) {
						stream.close();
					}
				}

			}

			fw = new FileWriter(dnaBindingProteins, true);

			final ProgressCounter counter = new ProgressCounter(proteinAccs.size(),
					ProgressPrintingType.PERCENTAGE_STEPS, 1);
			for (final Protein surfaceProtein : proteinAccs) {
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					System.out.println(printIfNecessary);
				}
				boolean accessibleRatio = false;
				for (final Peptide peptide : surfaceProtein.getPeptides()) {
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
					final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, surfaceProtein.getAcc());
					if (annotatedProtein.containsKey(surfaceProtein.getAcc())) {
						final Entry entry = annotatedProtein.get(surfaceProtein.getAcc());
						final List<String> query = JAXBXPathQuery.query(entry, "dbReference$id=GO:0003677", "$id");
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
		} catch (final IOException e) {

		}
		return ret;
	}

	@Test
	public void surfaceAccesibilityFromRemoteFilesMaxPerSite() {

		final File outputFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\C_max_per_site_surface.csv");
		final File uniprotReleasesFolder = new File("z:\\share\\Salva\\data\\uniprotKB");
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		UniprotProteinRetrievalSettings.getInstance(uniprotReleasesFolder, true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");
		final File fastaIndexFolder = new File("C:\\Users\\Salva\\Desktop\\tmp\\Pint\\fasta");
		final File projectConfigFile = new File(
				"z:\\share\\Salva\\data\\PINT projects\\molecular painting\\HS_Hek293.xml");
		final ImportCfgFileReader cfgReader = new ImportCfgFileReader();

		final Project project = cfgReader.getProjectFromCfgFile(projectConfigFile, fastaIndexFolder);
		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;
		int psmsNoAA = 0;
		int psmsWithAANoRatios = 0;
		int discardedContainingMoreThanOneK = 0;
		int validPSMs = 0;
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);

			final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
					removeOtherMolecules, true, pdbFolder);
			final Set<Protein> proteinAccs = getProteinAccsFromProject(project);
			final List<String> ratioScoreNames = getRatioScoreNamesFromProject(project);
			printHeader(fw, ratioScoreNames);
			final SurfaceProteinReportManager manager = new SurfaceProteinReportManager(calc);
			// everything is already calculated, so not go for it if not present
			manager.setCalculateIfNotPresent(false);
			final Map<String, SurfaceProteinReport> surfaceAccesibilityFromProteins = manager
					.getReportsFromProteins(proteinAccs);
			final Set<String> psmIds = new THashSet<String>();
			final Set<String> peptideSequences = new THashSet<String>();
			final Set<String> peptideSequencesValid = new THashSet<String>();
			final Set<String> uniquePositionsValid = new THashSet<String>();
			final Set<String> proteinsWithNoPDBInformation = new THashSet<String>();
			final Set<String> totalProteins = new THashSet<String>();
			final Set<Condition> conditions = project.getConditions();
			for (final Condition condition : conditions) {
				if (condition.getName().startsWith("HS")) {
					continue;
				}
				final Set<PSM> psMs = condition.getPSMs();
				for (final PSM psm : psMs) {
					final String peptideSequence = psm.getSequence();

					peptideSequences.add(peptideSequence);
					if (!psmIds.contains(psm.getIdentifier())) {
						psmIds.add(psm.getIdentifier());
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
						for (final Ratio ratio2 : psm.getRatios()) {
							if (ratio2.getDescription().equals("AREA_RATIO")) {
								ratio = ratio2;
							}
						}
						if (ratio == null) {
							psmsWithAANoRatios++;
							continue;
						}
						validPSMs++;

						final Set<String> accs = getProteinAccessions(psm.getProteins());
						totalProteins.addAll(accs);
						for (final String acc : accs) {
							if (surfaceAccesibilityFromProteins.containsKey(acc)) {
								peptideSequencesValid.add(peptideSequence);
								final SurfaceProteinReport surfaceAccesibilityProteinReport = surfaceAccesibilityFromProteins
										.get(acc);
								ProteinReportWriter.printReportForPsm(fw, ratio, psm.getIdentifier(), peptideSequence,
										surfaceAccesibilityProteinReport, aa, true);
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
		} catch (final IOException e) {
			e.printStackTrace();
			fail();
		}
	}

	private void printHeader(FileWriter fw, List<String> ratioScoreNames) throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append("ACC").append("\t").append("PSMID").append("\t").append("SEQ").append("\t")
				.append("Position in peptide").append("\t").append(JMolAtomReport.getStaticHeaders()).append("\t")
				.append("AREA_RATIO").append("\t");
		if (ratioScoreNames != null && !ratioScoreNames.isEmpty()) {
			for (final String ratioScoreName : ratioScoreNames) {
				sb.append(ratioScoreName).append("\t");
			}
		}

		sb.append("RATIO").append("\t");
		if (ratioScoreNames != null) {
			for (final String ratioScoreName : ratioScoreNames) {
				sb.append(ratioScoreName).append("\t");
			}
		}
		sb.append("\n");
		fw.write(sb.toString());
	}

	private Set<String> getProteinAccessions(Set<edu.scripps.yates.utilities.proteomicsmodel.Protein> proteins) {
		final Set<String> proteinAccessionsFromAccessions = new THashSet<String>();

		for (final edu.scripps.yates.utilities.proteomicsmodel.Protein protein : proteins) {
			proteinAccessionsFromAccessions.add(protein.getAccession());
		}
		return proteinAccessionsFromAccessions;
	}

	private Set<String> getSurfaceProteinAccessions(Set<Protein> proteins) {
		final Set<String> proteinAccessionsFromAccessions = new THashSet<String>();

		for (final Protein protein : proteins) {
			proteinAccessionsFromAccessions.add(protein.getAcc());
		}
		return proteinAccessionsFromAccessions;
	}

	private Set<Protein> getProteinAccsFromProject(Project project) {
		final Set<Protein> proteinAccessionsFromProject = new THashSet<Protein>();
		final Set<Condition> conditions = project.getConditions();
		final Map<String, Peptide> surfacePeptideMap = new THashMap<String, Peptide>();
		for (final Condition condition : conditions) {
			final Set<edu.scripps.yates.utilities.proteomicsmodel.Protein> proteins = condition.getProteins();
			for (final edu.scripps.yates.utilities.proteomicsmodel.Protein protein : proteins) {
				final Protein surfaceProtein = new Protein(FastaParser.getUniProtACC(protein.getAccession()));
				for (final PSM psm : protein.getPSMs()) {
					if (psm.getRatios() != null && !psm.getRatios().isEmpty()) {
						Peptide surfacePeptide = null;
						if (surfacePeptideMap.containsKey(psm.getIdentifier())) {
							surfacePeptide = surfacePeptideMap.get(psm.getIdentifier());
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
		final Set<String> ret = new THashSet<String>();
		final Set<Condition> conditions = project.getConditions();
		for (final Condition condition : conditions) {
			final Set<PSM> psms = condition.getPSMs();
			for (final PSM psm : psms) {
				final Set<Ratio> ratios = psm.getRatios();
				if (ratios != null) {
					for (final Ratio ratio : ratios) {
						if (ratio.getAssociatedConfidenceScore() != null) {
							ret.add(ratio.getAssociatedConfidenceScore().getScoreName());
						}
					}
				}
			}
		}
		final List<String> ratioScoreNames = new ArrayList<String>();
		ratioScoreNames.addAll(ret);
		Collections.sort(ratioScoreNames);
		return ratioScoreNames;
	}

	@Test
	public void getSurfaceAccebilityFromPDBSites() {
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
				new File("z:\\share\\Salva\\data\\uniprotKB"), true);
		final File pdbFolder = new File("z:\\share\\Salva\\data\\pdb");

		final char aa = 'K';
		final AtomType atomType = AtomType.NZ;
		final String protein = "5A5T";
		final String chain = null;
		final SurfaceCalculator calc = new SurfaceCalculator(uplr, aa, atomType, removeOtherChains,
				removeOtherMolecules, true, pdbFolder);
		final SurfaceProteinReport surfaceAccesibilityProteinReport = calc.getReportFromPDBModel(protein, chain);
		final List<SurfaceReport> reports = surfaceAccesibilityProteinReport.getReports();
		boolean first = true;
		for (final SurfaceReport report : reports) {
			if (first) {
				System.out.println(report.getToStringHeaders());
				first = false;
			}
			System.out.println(report);
		}

	}
}
