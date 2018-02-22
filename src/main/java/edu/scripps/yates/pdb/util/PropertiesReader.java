package edu.scripps.yates.pdb.util;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by Salva on 3/9/2016.
 */
public class PropertiesReader {

	public static final String PDB_FOLDER = "pdb_files_folder";
	public static final String UNIPROT_FOLDER = "uniprot_files_folder";
	public static final String UNIPROT_VERSION = "uniprot_version";
	public static final String AMINOACIDS = "AAs";
	public static final String ENZYME_NAME = "enzyme_name";
	public static final String MISSEDCLEAVAGES = "missedCleavages";
	public static final String SEMICLEAVAGE = "semiCleavage";
	public static final String INPUT_FILE = "input_file";
	public static final String INPUT_FILE_SEPARATOR = "input_file_separator";
	public static final String PEPTIDE_SEQUENCE_COLUMN_INDEX = "peptide_sequence_column_index";
	public static final String PEPTIDE_RATIO_COLUMN_INDEX = "peptide_ratio_column_index";
	public static final String PROTEIN_ACC_COLUMN_INDEX = "protein_acc_column_index";
	public static final String FASTA_FILE = "fasta_file";
	public static final String SKIP_FIRST_LINE = "skip_first_line";
	public static final String PEPTIDE_FILTER_REGEXP = "peptide_filter_regexp";
	public static final String IGNORE_PEPTIDE_NOT_FOUND_IN_DB = "ignore_peptide_not_found_in_db";
	public static final String CALCULATION_TYPE = "calculation_type";
	public static final String ONE_MODEL_PER_PROTEIN = "one_model_pre_protein";
	public static final String PDB_IDS = "pdb_ids";

	private static File file;

	private PropertiesReader() {

	}

	public static java.util.Properties getProperties(File file) throws IOException {
		PropertiesReader.file = file;
		return getProperties();
	}

	public static java.util.Properties getProperties() throws IOException {

		FileReader is = new FileReader(file);

		java.util.Properties prop = new java.util.Properties();

		prop.load(is);

		return prop;
	}

	public static String getPropertyValue(String propertyName) {
		try {
			return getProperties().getProperty(propertyName);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
