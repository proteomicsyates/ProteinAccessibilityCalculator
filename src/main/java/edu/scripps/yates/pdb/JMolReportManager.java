package edu.scripps.yates.pdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
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

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.surface.SurfaceCalculator;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public abstract class JMolReportManager<R extends ProteinReport<T>, T extends JMolAtomReport> {
	protected final static Logger log = Logger.getLogger(JMolReportManager.class);
	protected final Calculator<R, T> calculator;
	private final File file;
	private boolean loaded;
	private final Map<String, R> reportsByKey = new THashMap<String, R>();

	public JMolReportManager(Calculator<R, T> calculator) {
		this.calculator = calculator;
		this.calculator.setManager(this);
		file = new File(calculator.getPdbParserManager().getPdbFileManager().getParentPath().getAbsolutePath()
				+ File.separator + getFileName());
		log.info("Reports file located at: " + file.getAbsolutePath());
	}

	/**
	 * Gets the name of the file
	 * 
	 * @return
	 */
	public abstract String getFileName();

	public Map<String, R> getProteinReportsByProteins(Collection<Protein> proteins, Entry entry) {
		final Map<String, R> ret = new THashMap<String, R>();
		for (final Protein protein : proteins) {
			final R proteinReportByProtein = getProteinReportByProtein(protein.getAcc(), null, entry);
			if (proteinReportByProtein != null) {
				ret.put(protein.getAcc(), proteinReportByProtein);
			}
		}
		return ret;
	}

	public R getProteinReportByProtein(Protein protein, String pdbID, Entry entry) {
		return getProteinReportByProtein(protein.getAcc(), pdbID, entry);
	}

	public R getProteinReportByProtein(Protein protein, String pdbID) {
		return getProteinReportByProtein(protein.getAcc(), pdbID, null);
	}

	public R getProteinReportByProtein(String proteinAcc, String pdbID) {
		return getProteinReportByProtein(proteinAcc, pdbID, null);
	}

	public R getProteinReportByProtein(String proteinAcc, String pdbID, Entry entry) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map
		log.debug("Getting reports for protein: " + proteinAcc);

		final R proteinReport = calculator.getReportFromProtein(proteinAcc, pdbID, entry);
		if (proteinReport != null) {

			final TIntObjectHashMap<Set<T>> reportsByPositionInUniprotSeq = proteinReport
					.getReportsByPositionInUniprotSeq();
			for (final Set<T> reports : reportsByPositionInUniprotSeq.valueCollection()) {
				for (final T report : reports) {
					addReport(report);
				}
			}
			if (!proteinReport.getReportsByPositionInUniprotSeq().isEmpty()) {
				appendReportToFile(proteinReport);
			}
			return proteinReport;
		}
		return null;

	}

	public R getReportByKey(String reportKey) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map

		return reportsByKey.get(reportKey);

	}

	private void appendReportToFile(R proteinReport) {
		PrintWriter out = null;
		try {
			boolean writeHeader = false;
			if (!file.exists() || file.length() == 0l) {
				// write the header = true
				writeHeader = true;
			}
			final FileWriter fw = new FileWriter(file, true);
			final BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);

			boolean firstOne = true;
			final TIntObjectHashMap<Set<T>> positions = proteinReport.getReportsByPositionInUniprotSeq();
			final List<Integer> sortedPositions = new ArrayList<Integer>();
			for (final int position : positions.keys()) {
				sortedPositions.add(position);
			}
			Collections.sort(sortedPositions);
			for (final Integer position : sortedPositions) {
				final Set<T> reports = positions.get(position);
				for (final T report : reports) {
					if (!report.isStored()) {
						if (firstOne) {
							log.info("Appending report of protein " + proteinReport.getUniprotACC() + " to file");
							firstOne = false;
							if (writeHeader) {
								out.println(report.getToStringHeaders());
							}
						}
						out.println(report.toString());
					}
				}
			}
			if (!firstOne) {
				log.info("Report appended to file");
			}
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			out.close();
		}

	}

	/**
	 * Write all the information in memory to the file, overriding the content of
	 * the file
	 */
	public void dumpToFile() {
		PrintWriter out = null;
		try {
			log.info("Writting into file reports for " + reportsByKey.size() + " proteins");
			final FileWriter fw = new FileWriter(file, true);
			final BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);

			final List<String> keys = new ArrayList<String>();
			keys.addAll(reportsByKey.keySet());
			Collections.sort(keys);
			boolean firstOne = true;
			for (final String key : keys) {
				final R proteinReport = reportsByKey.get(key);
				final TIntObjectHashMap<Set<T>> positions = proteinReport.getReportsByPositionInUniprotSeq();
				final List<Integer> sortedPositions = new ArrayList<Integer>();
				for (final int position : positions.keys()) {
					sortedPositions.add(position);
				}

				Collections.sort(sortedPositions);
				for (final Integer position : sortedPositions) {
					final Set<T> reports = positions.get(position);
					for (final JMolAtomReport surfaceAccessibilityReport : reports) {
						if (firstOne) {
							log.info("Overriding report of protein " + proteinReport.getUniprotACC() + " to file");
							// write header for the first one
							out.println(surfaceAccessibilityReport.getToStringHeaders());
							firstOne = false;
						}
						out.println(surfaceAccessibilityReport.toString());
					}
				}
			}
			log.info("Reports writed at: " + file.getAbsolutePath());
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			out.close();
		}
	}

	private void loadReportsFromFile() {
		if (!loaded) {

			BufferedReader br = null;
			try {
				if (file.exists()) {
					// first read the file to get only the uniprot accessions
					// then call to annotate those proteins in uniprot
					// this will avoid the annotation of proteins one by one
					// when a new version of uniprot is released
					annotateAllProteinsInReports();
					final FileInputStream fstream = new FileInputStream(file);
					final DataInputStream in = new DataInputStream(fstream);
					br = new BufferedReader(new InputStreamReader(in));
					String strLine;
					// Read File Line By Line
					boolean firstLine = true;
					while ((strLine = br.readLine()) != null) {
						if (firstLine) {
							firstLine = false;
							continue;
						}
						final T report = parseFromLine(strLine);
						if (report != null) {
							report.setStored(true);
							addReport(report);
						}
					}
					loaded = true;
					log.info(" accesibilities numbers loaded from local file for " + reportsByKey.size() + " sites");
				}
			} catch (final IOException e) {
				e.printStackTrace();
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (final IOException e) {
						e.printStackTrace();
					}
				}
			}

		}

	}

	private void annotateAllProteinsInReports() {
		try {
			final UniprotProteinLocalRetriever uplr = calculator.getUplr();
			if (uplr != null) {
				log.info(
						"Getting Uniprot annotations for all the proteins in the stored reports, in order to get the protein sequences");
				if (file.exists()) {
					Stream<String> stream = null;
					try {
						stream = Files.lines(Paths.get(file.getAbsolutePath()));
						final Set<String> accs = stream.map(l -> l.split("\t")[4]).collect(Collectors.toSet());
						uplr.getAnnotatedProteins(null, accs);
					} finally {
						if (stream != null) {
							stream.close();
						}
					}

				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private void addReport(T report) {
		final String key = report.getReportKey();
		if (reportsByKey.containsKey(key)) {
			reportsByKey.get(key).addReport(report);
		} else {
			final Map<String, Entry> annotatedProtein = calculator.getUplr().getAnnotatedProtein(null,
					report.getUniprotACC());
			if (annotatedProtein.containsKey(report.getUniprotACC())) {
				final Entry entry = annotatedProtein.get(report.getUniprotACC());
				final String proteinSequence = SurfaceCalculator.getUniprotProteinSequence(entry);
				if (proteinSequence != null) {
					final R proteinReport = createProteinReport(report.getUniprotACC(), proteinSequence);
					proteinReport.addReport(report);
					reportsByKey.put(key, proteinReport);
				}
			}
		}
	}

	public abstract R createProteinReport(String proteinacc, String proteinSequence);

	public abstract T parseFromLine(String strLine);

	public Map<String, R> getReportsFromProteins(Collection<Protein> proteins) {
		final List<Protein> list = new ArrayList<Protein>();
		list.addAll(proteins);
		return getReportsFromProteins(list, null);
	}

	public Map<String, R> getReportsFromProteins(List<Protein> proteins, List<Entry> entries) {
		final Map<String, R> ret = new THashMap<String, R>();

		final ProgressCounter counter = new ProgressCounter(proteins.size(), ProgressPrintingType.PERCENTAGE_STEPS, 1);
		int index = 0;
		for (final Protein protein : proteins) {
			counter.increment();
			final String printIfNecessary = counter.printIfNecessary();
			if (!"".equals(printIfNecessary)) {
				log.info(printIfNecessary);
			}

			final R report = getProteinReportByProtein(protein, null, entries.get(index));
			if (report != null) {
				ret.put(protein.getAcc(), report);
				// dumpToFile();
			}
			index++;
		}
		return ret;
	}

	/**
	 * @param calculateIfNotPresent the calculateIfNotPresent to set
	 */
	public void setCalculateIfNotPresent(boolean calculateIfNotPresent) {
	}
}
