// $Id: InChIKey.java 2580 2009-05-20 15:40:19Z nguyenda $

package tripod.chem.indexer.util;

import java.io.*;
import java.util.*;

import chemaxon.struc.Molecule;
import chemaxon.formats.MolImporter;

public class InChIKey {

    static PrintStream errStream = new PrintStream (System.err);
    static PrintStream keyStream = new PrintStream (System.out);
    static PrintStream grpStream = new PrintStream (System.out);


    static boolean hasIsotopeLayer (String inchi) {
	String[] layers = inchi.split("/");
	for (String layer : layers) {
	    if (layer.charAt(0) == 'i') {
		return true;
	    }
	}
	return false;
    }

    static void doInput (InputStream is) throws Exception {
	BufferedReader br = new BufferedReader (new InputStreamReader (is));

	Map<String, Vector<String>> eqclass = new 
	    HashMap<String, Vector<String>>();
	
	for (String line; (line = br.readLine()) != null; ) {
	    String[] toks = line.split("\\t");

	    String smiles = "";
	    int index = toks[0].indexOf("SMILES=");
	    if (index >= 0) {
		smiles = toks[0].substring(index+7);
	    }
	    String inchi = toks[1].substring(6);
	    String ikey = toks[2].substring(9);

	    keyStream.println(smiles + "\t" + ikey);

	    String key = null;
	    /*
	    if (hasIsotopeLayer (inchi)) {
		// compare the whole key...
		key = ikey;
	    }
	    else*/ { // otherwise use the first 14 characters...
		key = ikey.split("-")[0];
	    }

	    Vector<String> mb = eqclass.get(key);
	    if (mb == null) {
		eqclass.put(key, mb = new Vector<String>());
	    }

	    toks = smiles.split("\\s");
	    if (toks.length > 1) {
		mb.add(toks[1]);
	    }
	}

	grpStream.println(eqclass.size() + " equivalence classes!");
	for (Map.Entry<String, Vector<String>> e : eqclass.entrySet()) {
	    grpStream.print(e.getKey());
	    for (String mb : e.getValue()) {
		grpStream.print(" " + mb);
	    }
	    grpStream.println();
	}
    }

    static void doError (InputStream is) throws Exception {
	BufferedReader br = new BufferedReader (new InputStreamReader (is));
	for (String line; (line = br.readLine()) != null; ) {
	    errStream.println(line);
	}
    }

    static public void main (String[] argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("Usage: InChIKey FILES...");
	    System.exit(1);
	}

	String basename = "InchIKey";
	/*
	keyStream = new PrintStream 
	    (new FileOutputStream (new File (basename + "_key.smi")));
	grpStream = new PrintStream
	    (new FileOutputStream (new File (basename + "_grp.txt")));
	errStream = new PrintStream
	    (new FileOutputStream (new File (basename + "_err.txt")));
	*/
	ProcessBuilder pb = new ProcessBuilder
	    ("stdinchi-1", "-AuxNone", "-Key", "-Tabbed", 
	     "-STDIO", "-SDF:SMILES");
	final Process InChI = pb.start();

	new Timer().schedule(new TimerTask () {
		public void run () {
		    try {
			doInput (InChI.getInputStream());
		    }
		    catch (Exception ex) { 
			errStream.println(ex.getMessage());
			ex.printStackTrace(); 
		    }
		}
	    }, 0);

	new Timer().schedule(new TimerTask () {
		public void run () {
		    try {
			doError (InChI.getErrorStream());
		    }
		    catch (Exception ex) { 
			errStream.println(ex.getMessage());
			ex.printStackTrace(); 
		    }
		}
	    }, 0);

	PrintStream ps = new PrintStream (InChI.getOutputStream());
	for (int i = 0; i < argv.length; ++i) {
	    MolImporter mi = new MolImporter (new FileInputStream (argv[i]));
	    for (Molecule mol = new Molecule (); mi.read(mol); ) {
		String name = mol.getName();
		MolStandardizer.removeSaltOrSolvent(mol);
		if (mol.getAtomCount() > 0) {
		    mol.setProperty
			("SMILES", mol.toFormat("smiles:q") + "\t" + name);
		    ps.print(mol.toFormat("sdf:-a"));
		    ps.flush();
		}
	    }
	    mi.close();
	}
	ps.close();
	/*
	keyStream.close();
	grpStream.close();
	errStream.close();
	*/
    }
}
