package edu.scripps.yates.pdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import com.compomics.dbtoolkit.io.EnzymeLoader;
import com.compomics.util.protein.Enzyme;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.dbindex.util.PeptideFilter;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.SurfacePeptide;
import edu.scripps.yates.pdb.model.SurfaceProtein;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityCalculator;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityProteinReport;
import edu.scripps.yates.pdb.util.FastaDigestionConfiguration;
import edu.scripps.yates.pdb.util.PropertiesReader;
import edu.scripps.yates.pdb.util.SurfaceAccesibilityReportWriter;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.files.FileUtils;

public class PTM2PDB {
	private final static Logger log = Logger.getLogger(PTM2PDB.class);

	public static void main(String[] args) {

		BufferedWriter writer = null;
		if (args.length == 0) {
			log.error("No parameters provided");
			System.exit(-1);
		}
		File propertiesFile = new File(args[0]);
		try {
			PropertiesReader.getProperties(propertiesFile);
			String pdbFolder = PropertiesReader.getPropertyValue(PropertiesReader.PDB_FOLDER);
			File parentPDBFolderContainer = new File(pdbFolder);

			String uniprotFolder = PropertiesReader.getPropertyValue(PropertiesReader.UNIPROT_FOLDER);
			File uniprotReleasesFolder = new File(uniprotFolder);
			UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
			String aas = PropertiesReader.getPropertyValue(PropertiesReader.AMINOACIDS);
			AtomType atomType = AtomType.NZ;
			String enzymeName = PropertiesReader.getPropertyValue(PropertiesReader.ENZYME_NAME);

			int missedCleavages = Integer.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.MISSEDCLEAVAGES));
			Enzyme enzyme = EnzymeLoader.loadEnzyme(enzymeName, String.valueOf(missedCleavages));
			char[] enzymeArray = enzyme.getCleavage();
			boolean semicleavage = Boolean.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.SEMICLEAVAGE));
			boolean ignorePeptidesNotFoundInDB = Boolean
					.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.IGNORE_PEPTIDE_NOT_FOUND_IN_DB));
			FastaDigestionConfiguration fastaDigestion = null;
			String inputFileName = PropertiesReader.getPropertyValue(PropertiesReader.INPUT_FILE);
			File inputFile = new File(inputFileName);
			if (!inputFile.exists()) {
				log.error("Input file " + inputFileName + " doesnt exist");
				System.exit(-1);
			}
			Path outputFile = Paths.get(inputFile.getParent() + File.separator
					+ FilenameUtils.getBaseName(inputFileName) + "_SURFACE_REPORT.txt");
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
			String skipHeaderString = PropertiesReader.getPropertyValue(PropertiesReader.SKIP_FIRST_LINE);
			if ("TRUE".equalsIgnoreCase(skipHeaderString)) {
				skipHeader = true;
			} else if ("FALSE".equalsIgnoreCase(skipHeaderString)) {
				skipHeader = false;
			} else {
				log.error(PropertiesReader.SKIP_FIRST_LINE + " property can only be TRUE or FALSE");
				System.exit(-1);
			}
			int peptideSequenceColumnIndex = -1;
			int peptideRatioColumnIndex = -1;
			int proteinAccessionColumnIndex = -1;
			try {
				peptideSequenceColumnIndex = Integer
						.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.PEPTIDE_SEQUENCE_COLUMN_INDEX));
			} catch (NumberFormatException e) {
				log.error(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX + " property has to be a number");
				System.exit(-1);
			}
			try {
				peptideRatioColumnIndex = Integer
						.valueOf(PropertiesReader.getPropertyValue(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX));
			} catch (NumberFormatException e) {
				log.error(PropertiesReader.PEPTIDE_RATIO_COLUMN_INDEX + " property has to be a number");
				System.exit(-1);
			}

			String peptideFilterRegexp = PropertiesReader.getPropertyValue(PropertiesReader.PEPTIDE_FILTER_REGEXP);
			// if fasta file is provided, then ignore the protein accession
			// column index
			String fastaFileName = PropertiesReader.getPropertyValue(PropertiesReader.FASTA_FILE);
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
				} catch (NumberFormatException e) {
					log.error(PropertiesReader.PROTEIN_ACC_COLUMN_INDEX + " property has to be a number");
					System.exit(-1);
				}
			}
			SurfaceAccessibilityCalculator sac = new SurfaceAccessibilityCalculator(uplr, aas, atomType, true, true,
					parentPDBFolderContainer);
			// read input file
			Map<String, SurfaceProtein> surfaceProteins = readInputFile(inputFile, separatorString, skipHeader,
					fastaDigestion, peptideSequenceColumnIndex, peptideRatioColumnIndex, proteinAccessionColumnIndex);
			log.info(surfaceProteins.size() + " proteins read");
			writer = Files.newBufferedWriter(outputFile, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
					StandardOpenOption.WRITE);
			writer.write(SurfaceAccesibilityReportWriter.getReportHeader() + "\n");
			int percentage = 0;
			int count = 0;
			for (SurfaceProtein protein : surfaceProteins.values()) {
				final SurfaceAccessibilityProteinReport surfaceAccesibilityFromProtein = sac
						.getSurfaceAccesibilityFromProtein(protein);
				final int currentPercentage = count++ / surfaceProteins.size();
				if (currentPercentage != percentage) {
					percentage = currentPercentage;
					log.info(percentage + "% proteins (" + count + "/" + surfaceProteins.size() + ")");
				}
				if (surfaceAccesibilityFromProtein == null) {
					log.info("No info for protein " + protein.getAcc());
					continue;
				}
				final Set<SurfacePeptide> peptides = protein.getPeptides();
				for (SurfacePeptide peptide : peptides) {

					SurfaceAccesibilityReportWriter.printReportForPeptide(writer, peptide,
							surfaceAccesibilityFromProtein, aas, false);

					writer.flush();
				}

			}

		} catch (Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
			System.exit(-1);
		} finally {

			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
					log.error(e.getMessage());
					System.exit(-1);
				}
			}

		}

	}

	private static Map<String, SurfaceProtein> readInputFile(File inputFile, String separator, boolean skipHeader,
			FastaDigestionConfiguration fastaDigestion, int peptideSequenceColumnIndex, int peptideRatioColumnIndex,
			int proteinAccessionColumnIndex) {
		int numProteins = 0;
		int numPeptides = 0;
		String ratioString = null;
		int row = 0;
		try {
			Map<String, SurfaceProtein> ret = new HashMap<String, SurfaceProtein>();

			List<String> peptideSequences = FileUtils.readColumnFromTextFile(inputFile, separator,
					peptideSequenceColumnIndex, skipHeader);
			List<String> peptideRatios = null;
			if (peptideRatioColumnIndex > -1) {
				peptideRatios = FileUtils.readColumnFromTextFile(inputFile, separator, peptideRatioColumnIndex,
						skipHeader);
			}

			if (proteinAccessionColumnIndex > -1 && fastaDigestion == null) {
				List<String> proteinAccs = FileUtils.readColumnFromTextFile(inputFile, separator,
						proteinAccessionColumnIndex, skipHeader);

				for (row = 0; row < peptideSequences.size(); row++) {
					String acc = proteinAccs.get(row);
					String sequence = peptideSequences.get(row);
					Double ratio = null;
					if (peptideRatios != null) {
						ratioString = peptideRatios.get(row);
						if ("'-Infinity".equals(ratioString)) {
							ratio = Double.NEGATIVE_INFINITY;
						} else if ("'Infinity".equals(ratioString)) {
							ratio = Double.POSITIVE_INFINITY;
						} else {
							ratio = Double.valueOf(ratioString);
						}
					}
					SurfaceProtein protein = null;
					if (ret.containsKey(acc)) {
						protein = ret.get(acc);
					} else {
						protein = new SurfaceProtein(acc);
						numProteins++;
						ret.put(acc, protein);
					}
					SurfacePeptide peptide = new SurfacePeptide(sequence, ratio);
					numPeptides++;
					protein.getPeptides().add(peptide);
				}
			} else {
				DBIndexInterface dbindex = getFastaDBIndex(fastaDigestion.getFasta(), fastaDigestion.getEnzymeArray(),
						fastaDigestion.getNumMisscleavages(), fastaDigestion.isSemiCleavage(),
						fastaDigestion.getPeptideFilter());
				for (row = 0; row < peptideSequences.size(); row++) {

					String sequence = peptideSequences.get(row);
					Double ratio = null;
					if (peptideRatios != null) {
						ratioString = peptideRatios.get(row);
						if ("'-Infinity".equals(ratioString)) {
							ratio = Double.NEGATIVE_INFINITY;
						} else if ("'Infinity".equals(ratioString)) {
							ratio = Double.POSITIVE_INFINITY;
						} else {
							ratio = Double.valueOf(ratioString);
						}
					}
					// remove Ptms or pre and post peptide characteres
					String cleanSequence = FastaParser.cleanSequence(sequence);
					final Set<IndexedProtein> proteins = dbindex.getProteins(cleanSequence);
					if (proteins == null || proteins.isEmpty()) {
						log.warn("Peptide '" + sequence + "' (row " + row
								+ ") is not found in FASTA with the provided digestion parameters: " + fastaDigestion);
						if (!fastaDigestion.isIgnorePeptidesNotFoundInDB()) {
							throw new PeptideNotFoundInDBIndexException(sequence
									+ " not found in Fasta DB. If you want to ignore that peptide, set ignorePeptidesNotFoundInDB=true in the input parameters");
						}
						continue;
					}
					for (IndexedProtein indexedProtein : proteins) {
						String acc = indexedProtein.getAccession();
						SurfaceProtein protein = null;
						if (ret.containsKey(acc)) {
							protein = ret.get(acc);
						} else {
							protein = new SurfaceProtein(acc);
							numProteins++;
							ret.put(acc, protein);
						}
						SurfacePeptide peptide = new SurfacePeptide(sequence, ratio);
						numPeptides++;
						protein.getPeptides().add(peptide);
					}

				}
			}
			return ret;
		} catch (NumberFormatException e) {
			e.printStackTrace();
			log.error("Error parsing some ratios. A not valid entry has been detected: '" + ratioString + "' at row "
					+ row);
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			log.error("Error reading some file: " + e.getMessage());
			System.exit(-1);
		} finally {
			log.info(numProteins + " proteins and " + numPeptides + " peptides readed");
		}
		return null;
	}

	private static DBIndexInterface getFastaDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages,
			boolean semicleavage, PeptideFilter peptideFilter) {
		if (fastaFile != null) {

			DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			String fastaIndexKey = IndexUtil.createFullIndexFileName(defaultDBIndexParams);

			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);

			if (peptideFilter != null) {
				((DBIndexSearchParamsImpl) defaultDBIndexParams).setPeptideFilter(peptideFilter);
			}
			DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);

			return dbIndex;
		}
		return null;
	}
}
