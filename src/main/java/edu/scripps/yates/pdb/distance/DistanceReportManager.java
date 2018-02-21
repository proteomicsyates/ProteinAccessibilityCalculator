package edu.scripps.yates.pdb.distance;

import edu.scripps.yates.pdb.JMolReportManager;

/**
 * Class that manages the distances reports in a persistence layer, in this
 * case, a text file. Any query over the distance will be looked up in the
 * persistence layer. If that is not present, it will call to the
 * {@link DistanceCalculator}
 *
 * @author Salva
 *
 */
public class DistanceReportManager extends JMolReportManager<DistanceProteinReport, DistanceReport> {
	private static final String FILE_NAME = "distances.csv";

	public DistanceReportManager(DistanceCalculator calc) {
		super(calc);

	}

	@Override
	public DistanceReport parseFromLine(String strLine) {
		return DistanceReport.getFromString(strLine);
	}

	@Override
	public String getFileName() {
		return FILE_NAME;
	}

	@Override
	public DistanceProteinReport createProteinReport(String proteinAcc, String proteinSequence) {
		return new DistanceProteinReport(proteinAcc, proteinSequence);
	}
}
