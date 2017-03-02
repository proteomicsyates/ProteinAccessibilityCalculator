package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLConnection;

import org.apache.log4j.Logger;

/**
 * Retrieves PDB files by pdb ID from the repository
 *
 * @author Salva
 *
 */
public class PDBFileRetriever {
	private final static Logger log = Logger.getLogger(PDBFileRetriever.class);
	private final static String PDB_URL_GZIP = "http://www.rcsb.org/pdb/files/";
	public final static String PDB_FILE_GIZP_EXTENSION = ".pdb.gz";
	private final static String PDB_URL_NO_COMPRESSED = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=";

	/**
	 * Gets a file with the PDB file in a temp file
	 *
	 * @param pdbID
	 * @return
	 */
	public static File getPDBFile(String pdbID) {
		return sendRequest(getPDBNoCompressedURL(pdbID));
	}

	/**
	 * Gets a file with the PDB file as GZIP in a temp file
	 *
	 * @param pdbID
	 * @return
	 */
	public static File getPDBGZipFile(String pdbID) {
		return sendRequest(getPDBGZipURL(pdbID));
	}

	private static String getPDBGZipURL(String pdbID) {
		return PDB_URL_GZIP + pdbID + PDB_FILE_GIZP_EXTENSION;
	}

	private static String getPDBNoCompressedURL(String pdbID) {
		return PDB_URL_NO_COMPRESSED + pdbID;
	}

	private static File sendRequest(String urlString) {

		try {
			// log.info("Submitting String= " + urlString + "...");
			URL url = new URL(urlString).toURI().toURL();
			log.info("Submitting URL= " + url + "...");
			long t1 = System.currentTimeMillis();
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
			HttpURLConnection.setFollowRedirects(true);
			conn.setDoInput(true);
			conn.connect();

			int status = conn.getResponseCode();
			while (true) {
				int wait = 0;
				String header = conn.getHeaderField("Retry-After");
				if (header != null)
					wait = Integer.valueOf(header);
				if (wait == 0)
					break;
				log.info("Waiting (" + wait + ")...");
				conn.disconnect();
				Thread.sleep(wait * 1000);
				conn = (HttpURLConnection) new URL(urlString).openConnection();
				conn.setDoInput(true);
				conn.connect();
				status = conn.getResponseCode();
			}
			if (status == HttpURLConnection.HTTP_OK) {
				long t2 = System.currentTimeMillis();
				log.info("Got a OK reply in " + (t2 - t1) / 1000 + "sg");
				InputStream is = conn.getInputStream();
				URLConnection.guessContentTypeFromStream(is);
				final File response = parseResponse(is);
				return response;
			} else
				log.error("Failed, got " + conn.getResponseMessage() + " for " + urlString);
			conn.disconnect();

		} catch (MalformedURLException e) {
			e.printStackTrace();
			log.warn(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			log.warn(e.getMessage());
		} catch (InterruptedException e) {
			e.printStackTrace();
			log.warn(e.getMessage());
		} catch (URISyntaxException e) {
			e.printStackTrace();
			log.warn(e.getMessage());
		}
		return null;
	}

	private static File parseResponse(InputStream is) {

		log.debug("Processing response");
		OutputStream outputStream = null;
		try {

			final File createTempFile = File.createTempFile("pdbTMP", "pdb");
			createTempFile.deleteOnExit();
			// read this file into InputStream

			// write the inputStream to a FileOutputStream
			outputStream = new FileOutputStream(createTempFile);

			int read = 0;
			byte[] bytes = new byte[1024];

			while ((read = is.read(bytes)) != -1) {
				outputStream.write(bytes, 0, read);
			}
			log.debug("Response parsed succesfully. File saved '" + createTempFile.getAbsolutePath() + "'");

			return createTempFile;
		} catch (IOException e2) {
			e2.printStackTrace();
		} finally {
			if (is != null) {
				try {
					is.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (outputStream != null) {
				try {
					// outputStream.flush();
					outputStream.close();
				} catch (IOException e) {
					e.printStackTrace();
				}

			}
		}
		return null;
	}
}
