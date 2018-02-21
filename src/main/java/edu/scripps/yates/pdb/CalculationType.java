package edu.scripps.yates.pdb;

public enum CalculationType {
	SURFACE, // surface accesibility from input Uniprot proteins
	PDB_SURFACE, // surface accessibility from PDB models
	DISTANCE;

	public static CalculationType fromValue(String propertyValue) {
		for (CalculationType calculationType : values()) {
			if (calculationType.name().equalsIgnoreCase(propertyValue)) {
				return calculationType;
			}
		}
		throw new IllegalArgumentException(
				propertyValue + " is not recognized as a valid value. Valid values: " + getValues());
	}

	public static String getValues() {
		StringBuilder sb = new StringBuilder();
		for (CalculationType calculationType : values()) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append(calculationType.name());
		}
		return sb.toString();
	}
}
