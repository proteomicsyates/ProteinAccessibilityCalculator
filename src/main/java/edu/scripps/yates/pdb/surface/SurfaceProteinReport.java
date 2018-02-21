package edu.scripps.yates.pdb.surface;

import edu.scripps.yates.pdb.ProteinReport;

/**
 * Stores all SurfaceAccesibilityReport from a protein, stored by position in
 * its sequence
 */
public class SurfaceProteinReport extends ProteinReport<SurfaceReport> {

	public SurfaceProteinReport(String uniprotACC, String uniprotProteinSeq) {
		super(uniprotACC, uniprotProteinSeq);
	}

}
