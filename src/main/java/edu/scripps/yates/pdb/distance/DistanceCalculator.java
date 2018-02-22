package edu.scripps.yates.pdb.distance;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pdb.CalculationType;
import edu.scripps.yates.pdb.Calculator;
import edu.scripps.yates.pdb.distance.model.Distance;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.read.PDBParser;
import edu.scripps.yates.pdb.util.InputParameters;

public class DistanceCalculator extends Calculator<DistanceProteinReport, DistanceReport> {
	protected final static Logger log = Logger.getLogger(DistanceCalculator.class);
	private final double distanceThreshold;

	public DistanceCalculator(UniprotProteinLocalRetriever uplr, Map<Character, List<AtomType>> atomTypeMap,
			boolean removeOtherChains, boolean removeOtherMolecules, boolean oneModelPerProtein,
			File parentPDBFolderContainer, double distanceThreshold) {
		super(uplr, atomTypeMap, removeOtherChains, removeOtherMolecules, oneModelPerProtein, parentPDBFolderContainer);
		this.distanceThreshold = distanceThreshold;
		setManager(new DistanceReportManager(this));

	}

	public DistanceCalculator(UniprotProteinLocalRetriever uplr, char aa, AtomType atomType, boolean removeOtherChains,
			boolean removeOtherMolecules, boolean oneModelPerProtein, File parentPDBFolderContainer,
			double distanceThreshold) {
		super(uplr, aa, atomType, removeOtherChains, removeOtherMolecules, oneModelPerProtein,
				parentPDBFolderContainer);
		this.distanceThreshold = distanceThreshold;
		setManager(new DistanceReportManager(this));

	}

	@Override
	public DistanceReport calculateReport(PDBParser parser, InputParameters inputParameters, Atom3D atom,
			int positionInPDB, Float resolution) {
		Set<Distance> distances = new HashSet<Distance>();
		Collection<Distance> distances2 = parser.getDistancesOfAtom(distanceThreshold, atom, "D", AtomType.OD1);
		distances.addAll(distances2);
		Collection<Distance> distances3 = parser.getDistancesOfAtom(distanceThreshold, atom, "D", AtomType.OD2);
		distances.addAll(distances3);
		Collection<Distance> distances4 = parser.getDistancesOfAtom(distanceThreshold, atom, "E", AtomType.OE1);
		distances.addAll(distances4);
		Collection<Distance> distances5 = parser.getDistancesOfAtom(distanceThreshold, atom, "E", AtomType.OE2);
		distances.addAll(distances5);
		if (distances != null && !distances.isEmpty()) {
			DistanceReport report = new DistanceReport(distances, inputParameters.getPdbID(), atom,
					inputParameters.getUniprotACC(), inputParameters.getPositionInUniprotProtein(), resolution,
					inputParameters.isRemoveOtherChains(), inputParameters.isRemoveOtherMolecules(),
					parser.getMutation(), parser.getExperimentalMethod());
			return report;
		}
		return null;
	}

	@Override
	public DistanceProteinReport createProteinReportObject(String acc, String proteinSequence) {
		return new DistanceProteinReport(acc, proteinSequence);
	}

	@Override
	public CalculationType getCalculationType() {
		return CalculationType.DISTANCE;
	}

	@Override
	public boolean isParseCoordinates() {
		return true;
	}
}
