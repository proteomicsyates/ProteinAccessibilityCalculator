package edu.scripps.yates.pdb.surface;

import java.io.File;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pdb.CalculationType;
import edu.scripps.yates.pdb.Calculator;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.read.PDBParser;
import edu.scripps.yates.pdb.util.InputParameters;

public class SurfaceCalculator extends Calculator<SurfaceProteinReport, SurfaceReport> {
	protected final static Logger log = Logger.getLogger(SurfaceCalculator.class);

	public SurfaceCalculator(Map<Character, List<AtomType>> atomTypeMap, boolean removeOtherChains,
			boolean removeOtherMolecules, File parentPDBFolderContainer) {
		super(atomTypeMap, removeOtherChains, removeOtherMolecules, parentPDBFolderContainer);
		setManager(new SurfaceProteinReportManager(this));
	}

	public SurfaceCalculator(UniprotProteinLocalRetriever uplr, Map<Character, List<AtomType>> atomTypeMap,
			boolean removeOtherChains, boolean removeOtherMolecules, boolean oneModelPerProtein,
			File parentPDBFolderContainer) {
		super(uplr, atomTypeMap, removeOtherChains, removeOtherMolecules, oneModelPerProtein, parentPDBFolderContainer);
		setManager(new SurfaceProteinReportManager(this));
	}

	public SurfaceCalculator(UniprotProteinLocalRetriever uplr, Character aa, AtomType atomType,
			boolean removeOtherChains, boolean removeOtherMolecules, boolean oneModelPerProtein,
			File parentPDBFolderContainer) {
		super(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, oneModelPerProtein,
				parentPDBFolderContainer);
		setManager(new SurfaceProteinReportManager(this));
	}

	@Override
	public SurfaceReport calculateReport(PDBParser parser, InputParameters inputParameters, Atom3D atom,
			int positionInPDB, Float resolution) {

		Double accesibility = parser.getSurfaceAccessibilityOfAtom(atom, inputParameters.isRemoveOtherChains(),
				inputParameters.isRemoveOtherMolecules());
		if (accesibility != null) {
			SurfaceReport report = new SurfaceReport(accesibility, inputParameters.getPdbID(), atom,
					inputParameters.getUniprotACC(), inputParameters.getPositionInUniprotProtein(), resolution,
					inputParameters.isRemoveOtherChains(), inputParameters.isRemoveOtherMolecules(),
					parser.getMutation(), parser.getExperimentalMethod());
			return report;
		}
		return null;
	}

	@Override
	public SurfaceProteinReport createProteinReportObject(String acc, String proteinSequence) {
		return new SurfaceProteinReport(acc, proteinSequence);
	}

	@Override
	public CalculationType getCalculationType() {
		return CalculationType.SURFACE;
	}

	@Override
	public boolean isParseCoordinates() {
		return false;
	}
}
