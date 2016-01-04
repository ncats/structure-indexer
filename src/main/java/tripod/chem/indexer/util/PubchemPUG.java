// $Id: PubchemPUG.java 2278 2008-05-29 22:27:45Z nguyenda $
// 06.19.07 pubchem's pug interface

package tripod.chem.indexer.util;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.OutputStream;
import java.io.ByteArrayInputStream;
import java.net.URL;
import java.net.URLConnection;
import java.net.HttpURLConnection;
import java.util.Vector;
import java.util.Collection;
import java.util.zip.GZIPInputStream;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;
import chemaxon.util.MolHandler;


public class PubchemPUG {
    public static final String PUG_URL = 
	"http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi";

    /*
     * all relevant paths
     */
    protected static final String PUG_PATH_OUTPUT = 
	"PCT-Data/PCT-Data_output/PCT-OutputData";
    protected static final String PUG_PATH_STATUS = PUG_PATH_OUTPUT 
	+ "/PCT-OutputData_status/PCT-Status-Message/PCT-Status-Message_status/PCT-Status";
    protected static final String PUG_PATH_REQID = PUG_PATH_OUTPUT
	+ "/PCT-OutputData_output/PCT-OutputData_output_waiting/PCT-Waiting/PCT-Waiting_reqid";
    protected static final String PUG_PATH_DOWNLOAD_URL = PUG_PATH_OUTPUT
	+ "/PCT-OutputData_output/PCT-OutputData_output_download-url/PCT-Download-URL/PCT-Download-URL_url";


    private String m_reqid = null;
    private String m_status = null;
    private String m_result = null;

    public PubchemPUG () {
    }

    protected static String downloadRequest (String database, 
					     Collection<Integer> uids) {
	return downloadRequest (database, uids, "sdf", "gzip");
    }

    // valid database: "pccompound", "pcsubstance", or "pcassay"
    // valid format: text-asn, binary-asn, xml, or sdf
    // valid compression: none, gzip, or bzip2
    protected static String downloadRequest 
	(String database, Collection<Integer> uids, 
	 String format, String compression) {
	StringBuffer request = new StringBuffer(
"<?xml version=\"1.0\"?>\n"+
"<!DOCTYPE PCT-Data PUBLIC \"-//NCBI//NCBI PCTools/EN\" \"http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd\">\n"+
"<PCT-Data>\n"+
"  <PCT-Data_input>\n"+
"    <PCT-InputData>\n"+
"      <PCT-InputData_download>\n"+
"        <PCT-Download>\n"+
"          <PCT-Download_uids>\n"+
"            <PCT-QueryUids>\n"+
"              <PCT-QueryUids_ids>\n"+
"                <PCT-ID-List>\n"+
"                  <PCT-ID-List_db>" + database + "</PCT-ID-List_db>\n"+
"                  <PCT-ID-List_uids>\n");
	for (Integer uid : uids) {
	    request.append 
		("                    <PCT-ID-List_uids_E>" 
		 + uid + "</PCT-ID-List_uids_E>\n");
	}
	request.append( 
"                  </PCT-ID-List_uids>\n"+
"                </PCT-ID-List>\n"+
"              </PCT-QueryUids_ids>\n"+
"            </PCT-QueryUids>\n"+
"          </PCT-Download_uids>\n"+
"          <PCT-Download_format value=\"" + format + "\"/>\n"+
"          <PCT-Download_compression value=\"" + compression + "\"/>\n"+
"        </PCT-Download>\n"+
"      </PCT-InputData_download>\n"+
"    </PCT-InputData>\n"+
"  </PCT-Data_input>\n"+
"</PCT-Data>\n");

	return request.toString();
    }

    protected static String pollRequest (String reqid) {
	String request = 
"<?xml version=\"1.0\"?>\n"+
"<!DOCTYPE PCT-Data PUBLIC \"-//NCBI//NCBI PCTools/EN\" \"http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd\">\n"+
"<PCT-Data>\n"+
"  <PCT-Data_input>\n"+
"    <PCT-InputData>\n"+
"      <PCT-InputData_request>\n"+
"        <PCT-Request>\n"+
"          <PCT-Request_reqid>";
	request += reqid + "</PCT-Request_reqid>\n";
	request +=
"          <PCT-Request_type value=\"status\"/>\n"+
"        </PCT-Request>\n"+
"      </PCT-InputData_request>\n"+
"    </PCT-InputData>\n"+
"  </PCT-Data_input>\n"+
"</PCT-Data>\n";
	return request;
    }

    protected void reset () {
	m_status = null;
	m_reqid = null;
	m_result = null;
    }

    protected void parseResponse (String response) {
	XmlTwig twig = new XmlTwig
	    (new ByteArrayInputStream (response.getBytes()));

	/* status is one of the following:

	   success	-  task completed successfully
	   server-error	-  request completion general server failure
	   hit-limit	-  success, but hit limit reached
	   time-limit	-  success, but time limit reached
	   input-error	-  input errors problem with input options
	   data-error	-  problem with input data
	   stopped	-  request management status request was stopped
	   running	-  request is running
	   queued	-  request is queued */
	m_status = twig.getElementAttrValue(PUG_PATH_STATUS, "value");
	m_reqid = twig.getElementValue(PUG_PATH_REQID);
	m_result = twig.getElementValue(PUG_PATH_DOWNLOAD_URL);
	if (System.getProperty("debug") != null) {
	    System.err.println(response);
	}
	System.err.println
	    ("status=" + m_status + " reqid=" + m_reqid + " url=" + m_result);
    }


    protected static Molecule[] downloadCompounds (String urlftp) 
	throws Exception {

	URL url = new URL (urlftp);
	URLConnection ftp = url.openConnection();

	MolImporter mi = new MolImporter 
	    (new GZIPInputStream (ftp.getInputStream()));
	Vector<Molecule> compounds = new Vector<Molecule>();
	for (Molecule mol; (mol = mi.read()) != null; ) {
	    //mol.hydrogenize(false);
	    //mol.aromatize();
	    //mol.calcHybridization();
	    compounds.add(mol);
	}
	mi.close();
	
	return compounds.toArray(new Molecule[0]);
    }

    protected static void downloadResults (OutputStream os, String urlftp) 
	throws Exception {
	URL url = new URL (urlftp);
	URLConnection ftp = url.openConnection();

	InputStream is = new GZIPInputStream (ftp.getInputStream());
	byte[] buf = new byte[2046]; 
	for (int nb; (nb = is.read(buf)) > 0; ) {
	    os.write(buf, 0, nb);
	}
	is.close();
    }

    protected String sendRequest (String request) throws Exception {
	String response = "";

	URL url = new URL (PUG_URL);
	HttpURLConnection http = (HttpURLConnection)url.openConnection();

	http.setRequestProperty("Content-Length", 
				String.valueOf(request.length()));
	http.setRequestProperty("Content-Type", "text/xml; charset=utf-8");
	http.setRequestMethod("POST");
	http.setDoOutput(true);
	http.setDoInput(true);

	if (System.getProperty("debug") != null) {
	    System.err.println(request);
	}
	OutputStream os = http.getOutputStream();
	os.write(request.getBytes());
	os.close();

	BufferedReader br = new BufferedReader 
	    (new InputStreamReader (http.getInputStream()));
	for (String line; (line = br.readLine()) != null; ) {
	    response += line + "\n";
	}
	br.close();

	http.disconnect();

	return response;
    }


    // get pubchem compound given list of IDs
    public void getResults (OutputStream os, String database, 
			    Collection<Integer> uids)  {
	try {
	    reset ();

	    String response = sendRequest 
		(downloadRequest (database, uids));
	    do {
		parseResponse (response);
		if (m_reqid != null) { // the job is queued..
		    System.err.println("waiting for reqid " + m_reqid);
		    Thread.sleep(500);
		    response = sendRequest (pollRequest (m_reqid));
		}
		else if (m_result != null) {
		    System.err.println("downloading from " + m_result);
		    downloadResults (os, m_result);

		    break;
		}
		else if (!"success".equals(m_status)) {
		    break;
		}
	    } while (true);
	}
	catch (Exception ex) {
	    ex.printStackTrace();
	}
    }

    public Molecule[] getCompounds (Collection<Integer> uids)  {
	Molecule[] mols = {};
	try {
	    reset ();

	    String response = sendRequest 
		(downloadRequest ("pccompound", uids));
	    do {
		parseResponse (response);
		if (m_reqid != null) { // the job is queued..
		    System.err.println("waiting for reqid " + m_reqid);
		    Thread.sleep(500);
		    response = sendRequest (pollRequest (m_reqid));
		}
		else if (m_result != null) {
		    System.err.println("downloading from " + m_result);
		    mols = downloadCompounds (m_result);

		    break;
		}
		else if (!"success".equals(m_status)) {
		    break;
		}
	    } while (true);
	}
	catch (Exception ex) {
	    ex.printStackTrace();
	}
	
	return mols;
    }

    public static void main (String argv[]) throws Exception {
	if (argv.length < 2) {
	    System.err.println("usage: PubchemPUG DB[pccompound,pcsubstance,pcassay] [-f FILE] [UIDS...]");
	    System.exit(1);
	}

	PubchemPUG pug = new PubchemPUG ();

	String db = argv[0];

	Vector<Integer> uids = new Vector<Integer>();
	int arg = 1;
	if (argv[arg].equals("-f")) {
	    if (argv.length < 3) {
		System.err.println("Not enough argument specified!");
		System.exit(1);
	    }
	    BufferedReader br = new BufferedReader
		(new FileReader (argv[++arg]));
	    for (String line; (line = br.readLine()) != null; ) {
		uids.add(Integer.parseInt(line));
	    }
	    br.close();
	    ++arg;
	}

	for (int i = arg; i < argv.length; ++i) {
	    uids.add(Integer.parseInt(argv[i]));
	}

	pug.getResults(System.out, db, uids);
	/*
	Molecule mol[] = pug.getCompounds(uids);
	for (int i = 0; i < mol.length; ++i) {
	    System.out.println(i + " " + mol[i].toFormat("smiles:u"));
	}
	*/
	
    }
}
