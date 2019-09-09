package edu.scripps.yates.pdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.springframework.core.io.ClassPathResource;

import com.compomics.dbtoolkit.io.EnzymeLoader;
import com.compomics.util.protein.Enzyme;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.util.FastaDigestionConfiguration;
import edu.scripps.yates.pdb.distance.DistanceCalculator;
import edu.scripps.yates.pdb.distance.DistanceProteinReport;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.read.InputFileReader;
import edu.scripps.yates.pdb.surface.SurfaceCalculator;
import edu.scripps.yates.pdb.surface.SurfaceProteinReport;
import edu.scripps.yates.pdb.util.PropertiesReader;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.properties.PropertiesUtil;

public class Main {
	private final static Logger log = Logger.getLogger(Main.class);
	private static AppVersion version;

	public static void main(String[] args) {
		final AppVersion version = Main.getVersion();
		System.out.println("Running PAC-calculator version " + version.toString());
		BufferedWriter writer = null;
		Path outputFile = null;
		if (args.length == 0) {
			log.error("No parameters provided");
			System.exit(-1);
		}
		final File propertiesFile = new File(args[0]);
		try {
			PropertiesReader.getProperties(propertiesFile);

			CalculationType calculationType = null;
			ProgressCounter counter = null;
			Map<String, Protein> proteins = null;
			try {
				calculationType = CalculationType
						.fromValue(PropertiesReader.getPropertyValue(PropertiesReader.CALCULATION_TYPE));
			} catch (final IllegalArgumentException e) {
				log.error("Invalid " + PropertiesReader.CALCULATION_TYPE + " value. Valid values are: "
						+ CalculationType.getValues());
				System.exit(-1);
			}
			final String aaAndAtomTypeString = PropertiesReader.getPropertyValue(PropertiesReader.AMINOACIDS);
			final Map<Character, List<AtomType>> atomTypeMap = new HashMap<Character, List<AtomType>>();
			final List<String> aaAndAtomTypeList = new ArrayList<String>();
			if (aaAndAtomTypeString.contains(",")) {
				for (final String aaanaAtomTypeTmp : aaAndAtomTypeString.split(",")) {
					aaAndAtomTypeList.add(aaanaAtomTypeTmp.trim());
				}
			} else {
				aaAndAtomTypeList.add(aaAndAtomTypeString);
			}
			for (final String aaAndAtomType : aaAndAtomTypeList) {
				if (aaAndAtomType.contains(" ")) {
					final String[] split = aaAndAtomType.split(" ");
					if (split[0].length() != 1) {
						log.error(
								"AMINOACIDS parameter is incorrect. A pair of aminacid character and atom type (separated by ' ') is incorrect. Aminoacid must be one letter.");
						System.exit(-1);
					}
					final char aa = split[0].charAt(0);
					final AtomType atomType = AtomType.getByName(split[1]);
					if (atomTypeMap.containsKey(aa)) {
						atomTypeMap.get(aa).add(atomType);
					} else {
						final List<AtomType> list = new ArrayList<AtomType>();
						list.add(atomType);
						atomTypeMap.put(aa, list);
					}
				} else {
					log.error("AMINOACIDS parameter is incorrect. It must contain a list (separated by '|') "
							+ " of pairs (separated by ' ') of aminacid character and atom type");
					System.exit(-1);
				}
			}

			final String pdbFolderString = PropertiesReader.getPropertyValue(PropertiesReader.PDB_FOLDER).trim();
			final File parentPDBFolder = new File(pdbFolderString);

			if (calculationType == CalculationType.PDB_SURFACE) {
				final List<String> pdbIDList = new ArrayList<String>();
				final String pdbIDsString = PropertiesReader.getPropertyValue(PropertiesReader.PDB_IDS);
				if (pdbIDsString != null && !"".equals(pdbIDsString)) {

					if (pdbIDsString.contains(",")) {
						for (final String pdbID : pdbIDsString.split(",")) {
							pdbIDList.add(pdbID.trim());
						}
					} else {
						pdbIDList.add(pdbIDsString);
					}
				} else {
					log.error(PropertiesReader.PDB_IDS + " parameter is required when "
							+ PropertiesReader.CALCULATION_TYPE + "=" + CalculationType.PDB_SURFACE);
					System.exit(-1);
				}
				outputFile = Paths.get(propertiesFile.getCanonicalPath() + File.separator + "PDB_SURFACE_REPORT.txt");
				writer = Files.newBufferedWriter(outputFile, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
						StandardOpenOption.WRITE);
				writer.write(JMolAtomReport.getStaticHeaders() + "\n");
				final SurfaceCalculator surfaceCalculator = new SurfaceCalculator(atomTypeMap, true, true,
						parentPDBFolder);
				for (final String pdbID : pdbIDList) {
					final SurfaceProteinReport surfaceAccesibilityReport = surfaceCalculator
							.getReportFromPDBModel(pdbID);
					if (surfaceAccesibilityReport == null) {
						log.info("No info for PDB model " + pdbID);
						continue;
					}

					ProteinReportWriter.printReportForPDB(writer, surfaceAccesibilityReport);

				}

				return;
			}

			final String uniprotFolder = PropertiesReader.getPropertyValue(PropertiesReader.UNIPROT_FOLDER);
			final File uniprotReleasesFolder = new File(uniprotFolder);
			String uniprotVersion = PropertiesReader.getPropertyValue(PropertiesReader.UNIPROT_VERSION);
			if ("".equals(uniprotVersion)) {
				uniprotVersion = null;
			}
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);

			final String fastaFileName = PropertiesReader.getPropertyValue(PropertiesReader.FASTA_FILE);

			final String enzymeName = PropertiesReader.getPropertyValue(PropertiesReader.ENZYME_NAME);
			char[] enzymeArray = null;
			int missedCleavages = 0;
			try {
				missedCleavages = Integer.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.MISSEDCLEAVAGES));
				if (enzymeName != null && !"".equals(enzymeName)) {
					final Enzyme enzyme = EnzymeLoader.loadEnzyme(enzymeName, String.valueOf(missedCleavages));

					if (enzyme != null) {
						enzymeArray = enzyme.getCleavage();
					}
				} else {
					String enzymeArrayString = PropertiesReader.getPropertyValue(PropertiesReader.ENZYME_ARRAY);
					if (enzymeArrayString.contains(",")) {
						enzymeArrayString = enzymeArrayString.replace(",", "").trim();
					}
					enzymeArray = enzymeArrayString.trim().toCharArray();
				}
			} catch (final NumberFormatException e) {
				if (fastaFileName != null) {
					throw new IllegalArgumentException(
							PropertiesReader.MISSEDCLEAVAGES + " parameter is missing when using FASTA file");
				}
			}
			Boolean semicleavage = null;
			Boolean ignorePeptidesNotFoundInDB = null;
			try {
				semicleavage = Boolean.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.SEMICLEAVAGE));
				ignorePeptidesNotFoundInDB = Boolean
						.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.IGNORE_PEPTIDE_NOT_FOUND_IN_DB));
			} catch (final Exception e) {
				throw new IllegalArgumentException("Error in parameters " + PropertiesReader.SEMICLEAVAGE + " or "
						+ PropertiesReader.IGNORE_PEPTIDE_NOT_FOUND_IN_DB + ". They are needed when using FASTA file");
			}
			final String inputFileName = PropertiesReader.getPropertyValue(PropertiesReader.INPUT_FILE);
			final File inputFile = new File(inputFileName);
			if (!inputFile.exists()) {
				log.error("Input file " + inputFileName + " doesnt exist");
				System.exit(-1);
			}

			String separatorString = PropertiesReader.getPropertyValue(PropertiesReader.INPUT_FILE_SEPARATOR);
			if ("TAB".equalsIgnoreCase(separatorString)) {
				separatorString = "\t";
			} else if ("COMMA".equalsIgnoreCase(separatorString)) {
				separatorString = ",";
			} else {
				log.error(PropertiesReader.INPUT_FILE_SEPARATOR + " property can only be COMMA or TAB");
				System.exit(-1);
			}
			boolean skipHeader = false;
			final String skipHeaderString = PropertiesReader.getPropertyValue(PropertiesReader.SKIP_FIRST_LINE);
			if ("TRUE".equalsIgnoreCase(skipHeaderString)) {
				skipHeader = true;
			} else if ("FALSE".equalsIgnoreCase(skipHeaderString)) {
				skipHeader = false;
			} else {
				log.error(PropertiesReader.SKIP_FIRST_LINE + " property can only be TRUE or FALSE");
				System.exit(-1);
			}
			boolean oneModelPerPRotein = true;
			final String oneModelPerProteinString = PropertiesReader
					.getPropertyValue(PropertiesReader.ONE_MODEL_PER_PROTEIN);
			if (!"".equals(oneModelPerProteinString)) {
				oneModelPerPRotein = Boolean.valueOf(oneModelPerProteinString);
			}

			int peptideSequenceColumnIndex = -1;
			int peptideRatioColumnIndex = -1;
			int proteinAccessionColumnIndex = -1;
			try {
				peptideSequenceColumnIndex = Integer
						.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.PEPTIDE_SEQUENCE_COLUMN_INDEX));
			} catch (final NumberFormatException e) {
				log.error(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX + " property has to be a number");
				System.exit(-1);
			}
			try {
				peptideRatioColumnIndex = Integer
						.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX));
			} catch (final NumberFormatException e) {
				log.error(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX + " property has to be a number");
				System.exit(-1);
			}

			final String peptideFilterRegexp = PropertiesReader
					.getPropertyValue(PropertiesReader.PEPTIDE_FILTER_REGEXP);
			// if fasta file is provided, then ignore the protein accession
			// column index
			FastaDigestionConfiguration fastaDigestion = null;

			File fastaFile = null;
			if (fastaFileName != null) {
				fastaFile = new File(fastaFileName);
				if (!fastaFile.exists()) {
					log.error("Fasta file " + fastaFileName + " doesnt exist");
					System.exit(-1);
				} else {
					fastaDigestion = new FastaDigestionConfiguration(fastaFile, enzymeArray, missedCleavages,
							semicleavage, peptideFilterRegexp, ignorePeptidesNotFoundInDB);
				}
			} else {
				try {
					proteinAccessionColumnIndex = Integer
							.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.PROTEIN_ACC_COLUMN_INDEX));
				} catch (final NumberFormatException e) {
					log.error(PropertiesReader.PROTEIN_ACC_COLUMN_INDEX + " property has to be a number");
					System.exit(-1);
				}
			}

			switch (calculationType) {
			case DISTANCE:
				final DistanceCalculator distanceCalculator = new DistanceCalculator(uplr, atomTypeMap, true, true,
						oneModelPerPRotein, parentPDBFolder, 2.0);
				distanceCalculator.setUniprotVersion(uniprotVersion);
				distanceCalculator.setDigestionConfiguration(fastaDigestion);
				// read input file
				proteins = InputFileReader.readInputFile(inputFile, separatorString, skipHeader, fastaDigestion,
						peptideSequenceColumnIndex, peptideRatioColumnIndex, proteinAccessionColumnIndex);
				log.info(proteins.size() + " proteins read");
				outputFile = Paths.get(inputFile.getParent() + File.separator + FilenameUtils.getBaseName(inputFileName)
						+ "_DISTANCE_REPORT.txt");
				writer = Files.newBufferedWriter(outputFile, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
						StandardOpenOption.WRITE);
				writer.write(ProteinReportWriter.getPeptideReportHeader(JMolAtomReport.getStaticHeaders()) + "\n");
				counter = new ProgressCounter(proteins.size(), ProgressPrintingType.PERCENTAGE_STEPS, 0);
				distanceCalculator.getUplr().getAnnotatedProteins(uniprotVersion, proteins.keySet());

				for (final Protein protein : proteins.values()) {
					counter.increment();
					for (final char aa : atomTypeMap.keySet()) {
						final DistanceProteinReport proteinReport = distanceCalculator.getReportFromProtein(protein,
								null, null);

						final String percentage = counter.printIfNecessary();
						if (percentage != null && !"".equals(percentage)) {
							log.info(percentage);
						}
						if (proteinReport == null) {
							log.info("No info for protein " + protein.getAcc());
							continue;
						}
						final Set<Peptide> peptides = protein.getPeptides();
						for (final Peptide peptide : peptides) {

							ProteinReportWriter.printReportForPeptide(writer, peptide, proteinReport, aa, false);
						}
					}
				}
				break;
			case SURFACE:
				final SurfaceCalculator surfaceCalculator = new SurfaceCalculator(uplr, atomTypeMap, true, true,
						oneModelPerPRotein, parentPDBFolder);
				surfaceCalculator.setUniprotVersion(uniprotVersion);
				surfaceCalculator.setDigestionConfiguration(fastaDigestion);
				// read input file
				proteins = InputFileReader.readInputFile(inputFile, separatorString, skipHeader, fastaDigestion,
						peptideSequenceColumnIndex, peptideRatioColumnIndex, proteinAccessionColumnIndex);
				log.info(proteins.size() + " proteins read");
				outputFile = Paths.get(inputFile.getParent() + File.separator + FilenameUtils.getBaseName(inputFileName)
						+ "_SURFACE_REPORT.txt");
				writer = Files.newBufferedWriter(outputFile, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
						StandardOpenOption.WRITE);
				writer.write(ProteinReportWriter.getPeptideReportHeader(JMolAtomReport.getStaticHeaders()) + "\n");
				counter = new ProgressCounter(proteins.size(), ProgressPrintingType.PERCENTAGE_STEPS, 0);
				surfaceCalculator.getUplr().getAnnotatedProteins(null, proteins.keySet());

				for (final Protein protein : proteins.values()) {
					counter.increment();
					for (final char aa : atomTypeMap.keySet()) {
						final SurfaceProteinReport surfaceAccesibilityFromProtein = surfaceCalculator
								.getReportFromProtein(protein, null, null);
						final String percentage = counter.printIfNecessary();
						if (percentage != null && !"".equals(percentage)) {
							log.info(percentage);
						}
						if (surfaceAccesibilityFromProtein == null) {
							log.info("No info for protein " + protein.getAcc());
							continue;
						}
						final Set<Peptide> peptides = protein.getPeptides();
						for (final Peptide peptide : peptides) {
							ProteinReportWriter.printReportForPeptide(writer, peptide, surfaceAccesibilityFromProtein,
									aa, false);
						}
					}
				}
				break;
			default:
				break;
			}

		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
			System.exit(-1);
		} finally {

			if (writer != null) {
				try {
					writer.close();
				} catch (final IOException e) {
					e.printStackTrace();
					log.error(e.getMessage());
					System.exit(-1);
				}
			}

		}

	}

	private static File getJarFolder() throws URISyntaxException {
		return new File(Main.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getParentFile();
	}

	public static AppVersion getVersion() {
		if (version == null) {
			try {
				final String tmp = PropertiesUtil
						.getProperties(new ClassPathResource(AppVersion.APP_PROPERTIES).getInputStream())
						.getProperty("assembly.dir");
				if (tmp.contains("v")) {
					version = new AppVersion(tmp.split("v")[1]);
				} else {
					version = new AppVersion(tmp);
				}
			} catch (final Exception e) {
				e.printStackTrace();
			}
		}
		return version;

	}
}
