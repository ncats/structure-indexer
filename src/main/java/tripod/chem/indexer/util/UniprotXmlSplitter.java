
package tripod.chem.indexer.util;

import java.io.*;
import java.util.zip.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import org.xml.sax.helpers.DefaultHandler;
import javax.xml.parsers.*;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;


/**
 * Split the massive Uniprot XML file into N entries per file.
 */
public class UniprotXmlSplitter extends DefaultHandler {
    static private final Logger logger = Logger.getLogger
	(UniprotXmlSplitter.class.getName());

    static final String XMLHEAD = "<?xml version='1.0' encoding='UTF-8'?>\n"
	+"<uniprot xmlns=\"http://uniprot.org/uniprot\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd\">\n";

    static final String XMLTAIL = "</uniprot>\n";

    int entriesPerFile;
    String prefix;
    File dir = new File (".");

    // parsing parameters
    long count; // total number of entries processed 
    int files;
    boolean entry;
    PrintStream out; // current output stream
    StringBuilder buffer;

    public UniprotXmlSplitter () {
	this (1000, "");
    }

    public UniprotXmlSplitter (String prefix) {
	this (1000, prefix);
    }

    public UniprotXmlSplitter (int entriesPerFile, String prefix) {
	setEntriesPerFile (entriesPerFile);
	setPrefix (prefix);
    }

    public void setOutputDir (File dir) { this.dir = dir; }
    public File getOutputDir () { return dir; }

    public void setEntriesPerFile (int entriesPerFile) {
	this.entriesPerFile = entriesPerFile;
    }
    public int getEntriesPerFile () { return entriesPerFile; }
    public void setPrefix (String prefix) { this.prefix = prefix; }
    public String getPrefix () { return prefix; }

    public void parse (InputStream is) throws Exception {
	SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
	parser.parse(is, this);
    }

    protected String getNextFile () {
	StringBuilder sb = new StringBuilder ();
	if (prefix != null && prefix.length() > 0) {
	    sb.append( prefix);
	}
	sb.append(String.format("%1$05d", ++files));
	sb.append(".xml.gz");
	return sb.toString();
    }

    protected void nextChunk () {
	if (out != null) {
	    out.print(XMLTAIL);
	    out.close();
	}

	File file = new File (dir, getNextFile ());
	try {
	    out = new PrintStream (new GZIPOutputStream 
				   (new FileOutputStream (file)), true);
	    out.print(XMLHEAD);
	    logger.info("processing new chunk " + count +" "+file.getName());
	}
	catch (IOException ex) {
	    logger.log(Level.SEVERE, 
		       "Can't open file "+file+" for writing!", ex);
	}
    }

    @Override
    public void startDocument () {
	count = 0;
	entry = false;
	files = 0;
	buffer = new StringBuilder ();
	nextChunk ();
    }

    @Override
    public void endDocument () {
	if (out != null) {
	    out.print(XMLTAIL);
	    out.close();
	}

	logger.info(count+" entries splitted across "
		    +((count+entriesPerFile-1)/entriesPerFile)+" file(s)!");
	out = null;
	buffer = null;
	entry = false;
    }

    @Override
    public void characters (char[] ch, int start, int len) {
	if (entry) {
	    buffer.append(escape (new String (ch, start, len)));
	}
    }

    @Override
    public void startElement (String uri, String localName, String qName, 
			      org.xml.sax.Attributes attrs) {
	if (qName.equals("entry")) {
	    if (count > 0 && count % entriesPerFile == 0) {
		nextChunk ();
	    }
	    ++count;
	    entry = true;
	}

	if (entry) {
	    buffer.append("<"+qName);
	    for (int i = 0; i < attrs.getLength(); ++i) {
		String s = escape (attrs.getValue(i));
		buffer.append(" " + attrs.getQName(i)+"=\""+s+"\"");
	    }
	    buffer.append(">");
	}
    }
    
    static String escape (String s) {
	return s.replaceAll("&", "&amp;")
	    .replaceAll("<", "&lt;")
	    .replaceAll(">", "&gt;")
	    .replaceAll("\"", "&quot;");
    }
    
    @Override
    public void endElement (String uri, String localName, String qName) {
	if (entry) {
	    buffer.append("</"+qName+">");
	}
	
	if (qName.equals("entry")) {
	    processChunk (buffer.toString());
	    buffer = new StringBuilder ();
	    entry = false;
	}
    }

    // override by subclass
    protected void processChunk (String chunk) {
	//System.out.println(chunk);
	out.println(chunk);
    }
    
    public long getCount () { return count; }

    public static void main (String[] argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("Usage: UniprotXmlSplitter FILE [OPTIONS]");
	    System.out.println
		("where OPTIONS can be one or more of the following:");
	    System.out.println("  split=N - split N entries per file (default: N = 1000)");
	    System.out.println("  prefix=STR - use STR as prefix to generated files (default: none)");
	    System.out.println("  dir=DIR - use DIR as output directory (default: .)");
	    System.exit(1);
	}
	
	UniprotXmlSplitter splitter = new UniprotXmlSplitter ();

	String file = argv[0];
	for (int i = 1; i < argv.length; ++i) {
	    String[] args = argv[i].split("=");
	    if (args.length == 2) {
		if (args[0].equalsIgnoreCase("split")) {
		    splitter.setEntriesPerFile(Integer.parseInt(args[1]));
		}
		else if (args[0].equalsIgnoreCase("prefix")) {
		    splitter.setPrefix(args[1]);
		}
		else if (args[0].equalsIgnoreCase("dir")) {
		    splitter.setOutputDir(new File (args[1]));
		}
	    }
	}

	logger.info("Input="+file+"; split="+splitter.getEntriesPerFile()
		    +"; prefix="+splitter.getPrefix()+"; dir="
		    +splitter.getOutputDir());
	splitter.parse(new GZIPInputStream (new FileInputStream (file)));
    }
}
