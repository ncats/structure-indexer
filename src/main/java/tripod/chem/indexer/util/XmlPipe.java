// $Id: XmlPipe.java 3212 2009-09-04 20:24:24Z nguyenda $

package tripod.chem.indexer.util;

import java.io.*;
import java.net.*;
import javax.xml.parsers.*;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.*;
import java.util.zip.*;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;

public class XmlPipe extends PipedInputStream {
    static final Logger logger = Logger.getLogger(XmlPipe.class.getName());

    PrintWriter out;
    String stopElement = null;

    public XmlPipe () throws Exception {
	out = new PrintWriter (new PipedOutputStream (this), false);
    }

    public void parse (InputStream is) throws Exception {
	SAXParser parser = SAXParserFactory.newInstance().newSAXParser();
	parser.parse(is, new XmlHandler ());
	out.close();
    }

    public void setStopElement (String stopElement) { 
	this.stopElement = stopElement; 
    }
    public String getStopElement () { return stopElement; }

    class XmlHandler extends DefaultHandler {
	Stack<String> stack = new Stack<String>();
	boolean done = false;

	@Override
	public void characters (char[] ch, int start, int length) {
	    if (!done) {
		out.write(ch, start, length);
	    }
	}
	
	public void endElement (String uri, String local, String qname) {
	    //logger.info("uri="+uri+" local="+local+" qname="+qname);
	    if (!done) {
		out.print("</"+(uri.length() >0?(uri+":"):"")+qname+">");
		String pname = stack.pop(); // pname = qname
	    }
	}

	@Override
	public void notationDecl (String name, String publicId, String sysId) {
	    logger.info("name="+name+" public="+publicId+" system="+sysId);
	}

	@Override
	public void processingInstruction (String target, String data) {
	    logger.info("target="+target+" data="+data);
	}

	@Override
	public void startDocument () {
	    out.println("<?xml version=\"1.0\"?>");
	    done = false;
	}

	@Override
	public void endDocument () {
	    out.flush();
	}

	@Override
	public void startElement (String uri, String local, 
				  String qname, Attributes attrs) {
	    if (done) return;

	    if (stopElement != null) {
		if (qname.equals(stopElement)) {
		    while (!stack.isEmpty()) {
			qname = stack.pop();
			out.println("</"+qname+">");
		    }
		    done = true;
		    return;
		}
	    }

	    StringBuilder sb = new StringBuilder ();
	    sb.append("<"+(uri.length()>0?(uri+":"):"")+qname);
	    for (int i = 0; i < attrs.getLength(); ++i) {
		sb.append(" "+attrs.getQName(i)+"=\""+attrs.getValue(i)+"\"");
	    }
	    sb.append(">");
	    stack.push(qname);
	    out.print(sb.toString());
	    //logger.info("uri="+uri+" local="+local+" qname="+qname);
	}

	@Override
	public void unparsedEntityDecl (String name, String publicId, 
					String sysId, String notation) {
	    logger.info("name="+name+" publicId="+publicId+" systemId="+sysId
			+" notation="+notation);
	}

	@Override
	public void error (SAXParseException ex) {
	    logger.log(Level.SEVERE, "SAX Exception", ex);
	}
    }

    public static void main (String[] argv) throws Exception {

	final XmlPipe xml = new XmlPipe ();
	xml.setStopElement("PC-AssaySubmit_data");

	final URL url = new URL
	    ("ftp://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/XML/1708.xml.gz");
	final ExecutorService service = Executors.newSingleThreadExecutor();
	service.execute(new Runnable () {
		public void run () {
		    try {
			xml.parse(new GZIPInputStream (url.openStream()));
		    }
		    catch (Exception ex) {
			ex.printStackTrace();
		    }
		}
	    });

	Reader reader = new InputStreamReader (xml);
	try {
	    char[] buf = new char[1024];
	    for (int nc; (nc = reader.read
			  (buf, 0, buf.length)) > 0; ) {
		System.out.print(new String (buf, 0, nc));
	    }
	}
	catch (IOException ex) {
	    //ex.printStackTrace();
	}

	service.shutdown();
    }
}
