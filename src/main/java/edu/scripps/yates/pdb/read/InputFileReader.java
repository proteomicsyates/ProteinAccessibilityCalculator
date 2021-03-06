package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import com.compomics.dbtoolkit.io.implementations.FASTADBLoader;
import com.compomics.dbtoolkit.io.interfaces.Filter;

import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.dbindex.util.FastaDigestionConfiguration;
import edu.scripps.yates.dbindex.util.PeptideNotFoundInDBIndexException;
import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.map.hash.THashMap;

public class InputFileReader {
	private final static Logger log = Logger.getLogger(InputFileReader.class);

	public static Map<String, Protein> readInputFile(File inputFile, String separator, boolean skipHeader,
			FastaDigestionConfiguration fastaDigestion, int peptideSequenceColumnIndex, int peptideRatioColumnIndex,
			int proteinAccessionColumnIndex) {
		int numProteins = 0;
		int numPeptides = 0;
		String ratioString = null;
		int row = 0;
		try {
			final Map<String, Protein> ret = new THashMap<String, Protein>();

			final List<String> peptideSequences = FileUtils.readColumnFromTextFile(inputFile, separator,
					peptideSequenceColumnIndex, skipHeader);
			List<String> peptideRatios = null;
			if (peptideRatioColumnIndex > -1) {
				peptideRatios = FileUtils.readColumnFromTextFile(inputFile, separator, peptideRatioColumnIndex,
						skipHeader);
			}

			if (proteinAccessionColumnIndex > -1 && fastaDigestion == null) {
				final List<String> proteinAccs = FileUtils.readColumnFromTextFile(inputFile, separator,
						proteinAccessionColumnIndex, skipHeader);

				for (row = 0; row < peptideSequences.size(); row++) {
					final String acc = proteinAccs.get(row);
					// get rawsequence. it can be a peptide node string such as
					// PEPTIDE1_PEPTIDE2_PEPTIDE3
					final String rawSequence = peptideSequences.get(row);
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
					Protein protein = null;
					if (ret.containsKey(acc)) {
						protein = ret.get(acc);
					} else {
						protein = new Protein(acc);
						numProteins++;
						ret.put(acc, protein);
					}
					final List<String> sequences = new ArrayList<String>();
					if (rawSequence.contains("_")) {
						for (final String seq : rawSequence.split("_")) {
							sequences.add(seq);
						}
					} else {
						sequences.add(rawSequence);
					}
					for (final String seq : sequences) {
						final Peptide peptide = new Peptide(seq, ratio);
						numPeptides++;
						protein.getPeptides().add(peptide);
					}

				}
			} else {
				if (fastaDigestion.getEnzymeArray() != null) {
					final DBIndexImpl dbindex = FastaDigestionConfiguration.getFastaDBIndex(fastaDigestion);
					for (row = 0; row < peptideSequences.size(); row++) {

						final String rawSequence = peptideSequences.get(row);
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
						final List<String> sequences = new ArrayList<String>();
						if (rawSequence.contains("_")) {
							for (final String seq : rawSequence.split("_")) {
								sequences.add(seq);
							}
						} else {
							sequences.add(rawSequence);
						}
						for (final String sequence : sequences) {
//							
//						
//						// remove Ptms or pre and post peptide characteres
//						// if the peptide is a peptide node it could have several peptide sequences
//						// referring to the same quantitation site
//						// so we need to know if contains several sequences and in that case, just take
//						// one
//						if (sequence.contains("_")) {
//							sequence = sequence.split("_")[0];
//						}
							final String cleanSequence = FastaParser.cleanSequence(sequence);
							final Set<IndexedProtein> proteins = dbindex.getProteins(cleanSequence);
							if (proteins == null || proteins.isEmpty()) {
								log.warn("Peptide '" + sequence + "' (row " + row
										+ ") is not found in FASTA with the provided digestion parameters: "
										+ fastaDigestion);
								if (!fastaDigestion.isIgnorePeptidesNotFoundInDB()) {
									throw new PeptideNotFoundInDBIndexException(sequence
											+ " not found in Fasta DB. If you want to ignore that peptide, set ignorePeptidesNotFoundInDB=true in the input parameters");
								}
								continue;
							}
							for (final IndexedProtein indexedProtein : proteins) {
								final String acc = indexedProtein.getAccession();
								Protein protein = null;
								if (ret.containsKey(acc)) {
									protein = ret.get(acc);
								} else {
									protein = new Protein(acc);
									numProteins++;
									ret.put(acc, protein);
								}
								final Peptide peptide = new Peptide(sequence, ratio);
								numPeptides++;
								protein.getPeptides().add(peptide);
							}
						}
					}
				} else {
					log.info(
							"No digestion enzyme provided. Reading the entire Fasta to find peptide to proteins relationships");
					final File fastaFile = fastaDigestion.getFasta();
					final Map<String, String> proteinsByProteinSequences = loadFastaInMemory(fastaFile);
					final List<String> proteinSequences = new ArrayList<String>();
					proteinSequences.addAll(proteinsByProteinSequences.keySet());
					final ProgressCounter counter = new ProgressCounter(peptideSequences.size(),
							ProgressPrintingType.PERCENTAGE_STEPS, 0, true);
					for (row = 0; row < peptideSequences.size(); row++) {
						counter.increment();
						final String printIfNecessary = counter.printIfNecessary();
						if (!"".contentEquals(printIfNecessary)) {
							log.info(printIfNecessary + " peptides processed");
						}
						final String rawSequence = peptideSequences.get(row);
						final List<String> sequences = new ArrayList<String>();
						if (rawSequence.contains("_")) {
							for (final String seq : rawSequence.split("_")) {
								sequences.add(seq);
							}
						} else {
							sequences.add(rawSequence);
						}
						for (final String sequence : sequences) {
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
							final Set<String> proteinAccs = new HashSet<String>();
							final String cleanSequence = FastaParser.cleanSequence(sequence);
							for (final String proteinSequence : proteinSequences) {
								if (proteinSequence.contains(cleanSequence)) {
									proteinAccs.add(proteinsByProteinSequences.get(proteinSequence));
								}
							}
							if (proteinAccs == null || proteinAccs.isEmpty()) {
								log.warn("Peptide '" + sequence + "' (row " + row
										+ ") is not found in FASTA (brute force search of peptide)");
								if (!fastaDigestion.isIgnorePeptidesNotFoundInDB()) {
									throw new PeptideNotFoundInDBIndexException(sequence
											+ " not found in Fasta DB. If you want to ignore that peptide, set ignorePeptidesNotFoundInDB=true in the input parameters");
								}
								continue;
							}
							for (final String acc : proteinAccs) {
								Protein protein = null;
								if (ret.containsKey(acc)) {
									protein = ret.get(acc);
								} else {
									protein = new Protein(acc);
									numProteins++;
									ret.put(acc, protein);
								}
								final Peptide peptide = new Peptide(sequence, ratio);
								numPeptides++;
								protein.getPeptides().add(peptide);
							}
						}
					}
				}
			}
			return ret;
		} catch (final NumberFormatException e) {
			e.printStackTrace();
			log.error("Error parsing some ratios. A not valid entry has been detected: '" + ratioString + "' at row "
					+ row);
			System.exit(-1);
		} catch (final IOException e) {
			e.printStackTrace();
			log.error("Error reading some file: " + e.getMessage());
			System.exit(-1);
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			log.error("Error reading some file: " + e.getMessage());
			System.exit(-1);
		} finally {
			log.info(numProteins + " proteins and " + numPeptides + " peptides readed");
		}
		return null;
	}

	/**
	 * Reads a fasta file, and returns a map in which the key is the protein
	 * sequence and the value is the accession.<br>
	 * We assume that all proteins have different sequences.
	 * 
	 * @param fastaFile
	 * @return
	 * @throws IOException
	 */
	private static Map<String, String> loadFastaInMemory(File fastaFile) throws IOException {
		final long t1 = System.currentTimeMillis();
		log.info("Loading entire Fasta file in memory...");
		final FASTADBLoader fastadbLoader = new FASTADBLoader();
		if (!fastadbLoader.canReadFile(fastaFile)) {
			throw new IllegalArgumentException(fastaFile.getAbsolutePath() + " cannot be read");
		}
		fastadbLoader.load(fastaFile.getAbsolutePath());
		final Map<String, String> map = new HashMap<String, String>();
		final Filter decoyFilter = new Filter() {

			@Override
			public boolean passesFilter(HashMap aEntry) {
				// TODO Auto-generated method stub
				return false;
			}

			@Override
			public boolean passesFilter(String aEntry) {
				if (!FastaParser.isReverse(aEntry)) {
					return true;
				}
				return false;
			}
		};
		String protein = null;
		final ProgressCounter counter = new ProgressCounter(
				Long.valueOf(fastadbLoader.countNumberOfEntries()).intValue(), ProgressPrintingType.PERCENTAGE_STEPS, 0,
				true);
		while ((protein = fastadbLoader.nextFilteredRawEntry(decoyFilter)) != null) {
			counter.increment();
			final String printIfNecessary = counter.printIfNecessary();
			if (!"".equals(printIfNecessary)) {
				log.info(printIfNecessary + " entries loaded");
			}
			final String acc = FastaParser.getACC(protein.split("\n")[0]).getAccession();
			final String sequence = protein.substring(protein.indexOf("\n") + 1);
			map.put(sequence, acc);
		}
		log.info("Fasta file loaded in memory in "
				+ DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		return map;
	}
}
