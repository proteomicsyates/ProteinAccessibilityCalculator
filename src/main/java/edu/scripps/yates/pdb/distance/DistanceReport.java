package edu.scripps.yates.pdb.distance;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.distance.model.Distance;
import edu.scripps.yates.pdb.model.Atom3D;

public class DistanceReport extends JMolAtomReport {
	private final List<Distance> distances = new ArrayList<Distance>();

	public DistanceReport(Collection<Distance> distances, String pdbID, Atom3D atom, String uniprotACC,
			int positionInUniprot, Float resolution, boolean otherChainsRemoved, boolean otherMoleculesRemoved,
			boolean mutation, String method) {
		super(pdbID, atom, uniprotACC, positionInUniprot, resolution, otherChainsRemoved, otherMoleculesRemoved,
				mutation, method);
		if (distances != null) {
			this.distances.addAll(distances);
		}
	}

	public DistanceReport(Collection<Distance> distances, JMolAtomReport report) {
		this(distances, report.getPdbID(), report.getAtom(), report.getUniprotACC(), report.getPositionInUniprot(),
				report.getResolution(), report.isOtherChainsRemoved(), report.isOtherMoleculesRemoved(),
				report.isMutation(), report.getMethod());
	}

	/**
	 * @return the distances
	 */
	public List<Distance> getDistances() {
		return distances;
	}

	/**
	 * Parse a string as opposite function of toString()
	 * sb.append("pdb_ID").append(sep).append("uniprot_ACC").append(sep).append("position_in_uniprot").append(sep)
	 * .append("resolution").append(sep).append(Atom3D.getToStringHeaders()).append(sep)
	 * .append("OtherChainsRemoved").append(sep).append("OtherMoleculesRemoved").append(sep).append("Mutated")
	 * .append(sep).append("Method");
	 *
	 * @param string
	 * @return
	 */
	public static DistanceReport getFromString(String string) {
		try {
			List<String> lines = new ArrayList<String>();
			if (string.contains("\n")) {
				for (String split : string.split("\n")) {
					lines.add(split);
				}
			} else {
				lines.add(string);
			}
			DistanceReport report = null;
			boolean first = true;
			List<Distance> distances = new ArrayList<Distance>();
			for (String line : lines) {
				if (first) {

					final String[] split = string.split(sep);
					int index = 0;
					String pdbID = split[index++];
					String uniprotACC = split[index++];
					int positionInUniprot = Integer.valueOf(split[index++]);
					Float resolution = null;
					try {
						resolution = Float.valueOf(split[index++]);
					} catch (Exception e) {
					}

					// // get the rest of the splitted items to construct the
					// Atom object
					StringBuilder sb = new StringBuilder();
					for (int i = 0; i < Atom3D.numElementsInPrint(); i++) {
						sb.append(split[index++]).append(sep);
					}

					Atom3D atom = Atom3D.getFromString(sb.toString());
					Boolean removeOtherChains = Boolean.valueOf(split[index++]);
					Boolean removeOtherMolecules = Boolean.valueOf(split[index++]);
					Boolean containsMutation = Boolean.valueOf(split[index++]);
					String method = split[index++];

					report = new DistanceReport(null, pdbID, atom, uniprotACC, positionInUniprot, resolution,
							removeOtherChains, removeOtherMolecules, containsMutation, method);

					String distanceString = line.substring(report.toString().length());
					Distance distance = Distance.fromString(distanceString);
					distances.add(distance);
				}
			}
			DistanceReport ret = new DistanceReport(distances, report);
			return ret;
		} catch (

		Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public String getToStringHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append(JMolAtomReport.getStaticHeaders()).append(sep).append(Distance.getStaticHeaders());
		return sb.toString();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Distance distance : distances) {
			sb.append(super.toString()).append(sep).append(distance).append("\n");
		}
		return sb.toString();
	}

}
