// $Id: PubChemParseAssayContainer.java 2278 2008-05-29 22:27:45Z nguyenda $
package tripod.chem.indexer.util;

import java.io.InputStream;
import java.io.ByteArrayInputStream;
import java.util.Vector;
import java.util.zip.*;
import java.io.IOException;
import java.io.FileInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.Document;

// parse pubchem assay description xml
public class PubChemParseAssayContainer {
    private XmlTwig m_twig;

    public PubChemParseAssayContainer (InputStream is) {
	m_twig = new XmlTwig (is);
    }

    public int getAID () {
	String aid = m_twig.getElementValue("PC-AssayContainer/PC-AssaySubmit/PC-AssaySubmit_assay/PC-AssaySubmit_assay_descr/PC-AssayDescription/PC-AssayDescription_aid/PC-ID/PC-ID_id");
	if (aid != null) {
	    return Integer.parseInt(aid);
	}
	return -1;
    }

    public String getSource () {
	return m_twig.getElementValue("PC-AssayContainer/PC-AssaySubmit/PC-AssaySubmit_assay/PC-AssaySubmit_assay_descr/PC-AssayDescription/PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_name");
    }

    public String getSourceID () {
	return m_twig.getElementValue("PC-AssayContainer/PC-AssaySubmit/PC-AssaySubmit_assay/PC-AssaySubmit_assay_descr/PC-AssayDescription/PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_source-id/Object-id/Object-id_str");
    }

    public String getName () {
	return m_twig.getElementValue("PC-AssayContainer/PC-AssaySubmit/PC-AssaySubmit_assay/PC-AssaySubmit_assay_descr/PC-AssayDescription/PC-AssayDescription_name");
    }

    public String getDescription () {
	return getText ("PC-AssayDescription_description_E");
    }

    public String getProtocol () {
	return getText ("PC-AssayDescription_protocol_E");
    }

    public String getComment () {
	return getText ("PC-AssayDescription_comment_E");
    }

    protected String getText (String tag) {
	StringBuffer buf = new StringBuffer ();
	NodeList elms = m_twig.getDocument().getElementsByTagName(tag);
	for (int i = 0; i < elms.getLength(); ++i) {
	    Node node = (Element)elms.item(i);
	    NodeList childs = node.getChildNodes();
	    for (int j = 0; j < childs.getLength(); ++j) {
		Node text = childs.item(j);
		if (Node.TEXT_NODE == text.getNodeType()) {
		    buf.append(text.getNodeValue());
		}
	    }
	    buf.append("\n");
	}
	return buf.toString();
    }

    public int[] getXrefPMIDList () {
	String text[] = getText ("PC-XRefData_pmid").split("\n");
	int pmids[] = new int[text.length];
	for (int i = 0; i < text.length; ++i) {
	    pmids[i] = Integer.parseInt(text[i]);
	}
	return pmids;
    }
    public String getXrefPMID () { return getText ("PC-XRefData_pmid"); }

    public int[] getXrefAIDList () {
	String text[] = getText ("PC-XRefData_aid").split("\n");
	int aids[] = new int[text.length];
	for (int i = 0; i < text.length; ++i) {
	    aids[i] = Integer.parseInt(text[i]);
	}
	return aids;
    }
    public String getXrefAID () { return getText ("PC-XRefData_aid"); }


    public String getXrefUrl () {
	return getText ("PC-XRefData_dburl");
    }

    public String getTargetName () {
	return m_twig.getElementValue("PC-AssayContainer/PC-AssaySubmit/PC-AssaySubmit_assay/PC-AssaySubmit_assay_descr/PC-AssayDescription/PC-AssayDescription_target/PC-AssayTargetInfo/PC-AssayTargetInfo_name");
    }

    public static void main (String argv[]) throws Exception {
	if (argv.length == 0) {
	    System.out.println("usage: PubChemParseAssayContainer ASSAY_DESC...");
	    System.exit(1);
	}

	for (int i = 0; i < argv.length; ++i) {
	    PubChemParseAssayContainer assay;
	    try {
		assay = new PubChemParseAssayContainer 
		    (new GZIPInputStream (new FileInputStream (argv[i])));
	    }
	    catch (Exception ex) {
		try {
		    assay = new PubChemParseAssayContainer 
			(new ZipInputStream (new FileInputStream (argv[i])));
		}
		catch (Exception exx) {
		    assay = new PubChemParseAssayContainer 
			(new FileInputStream (argv[i]));
		}
	    }


	    System.out.println("*** " + argv[i] + " ****");
	    System.out.println("AID: " + assay.getAID());
	    System.out.println("Name: " + assay.getName());
	    System.out.println("Source: " + assay.getSource());
	    System.out.println("SourceID: " + assay.getSourceID());
	    System.out.println("Description: \n" + assay.getDescription());
	    System.out.println("Protocol:\n" + assay.getProtocol());
	    System.out.println("URL: " + assay.getXrefUrl());
	    System.out.println("Xref PMID: \n" + assay.getXrefPMID());
	    System.out.println("Xref AID: \n" + assay.getXrefAID());
	    System.out.println("Target: " + assay.getTargetName());
	}
    }
}
