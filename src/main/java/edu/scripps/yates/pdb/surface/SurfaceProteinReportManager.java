package edu.scripps.yates.pdb.surface;

import edu.scripps.yates.pdb.JMolReportManager;

/**
 * Class that manages the surface accessibility reports in a persistence layer,
 * in this case, a text file. Any query over the surface accessibility will be
 * looked up in the persistence layer. If that is not present, it will call to
 * the {@link SurfaceAccesibilityCalculator}
 *
 * @author Salva
 *
 */
public class SurfaceProteinReportManager extends JMolReportManager<SurfaceProteinReport, SurfaceReport> {
	private static final String FILE_NAME = "surfaces_accessibilities.csv";

	public SurfaceProteinReportManager(SurfaceCalculator calc) {
		super(calc);

	}

	@Override
	public SurfaceReport parseFromLine(String strLine) {
		return SurfaceReport.getFromString(strLine);
	}

	@Override
	public String getFileName() {
		return FILE_NAME;
	}

	@Override
	public SurfaceProteinReport createProteinReport(String proteinAcc, String proteinSequence) {
		return new SurfaceProteinReport(proteinAcc, proteinSequence);
	}

}
