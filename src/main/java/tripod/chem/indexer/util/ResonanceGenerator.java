// $Id: ResonanceGenerator.java 2716 2009-06-23 21:34:00Z nguyenda $

package tripod.chem.indexer.util;

import java.util.*;
import java.io.PrintStream;
import java.util.logging.Logger;

import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.formats.MolImporter;
import chemaxon.util.MolHandler;

/**
 * simple resonance generator
 */
public class ResonanceGenerator {
    private static boolean debug = false;
    static {
	try {
	    debug = Boolean.getBoolean("resonance.debug");
	}
	catch (Exception ex) {
	}
    }

    private static final Logger logger = Logger.getLogger
	(ResonanceGenerator.class.getName());

    private long timeout = 30000; // 30s timeout
    private long start = 0;

    private Molecule mol; 
    private Molecule origmol; // original input
    private MolAtom[] atoms;
    private MolBond[] bonds;

    private Vector<Molecule> resonances = new Vector<Molecule>();

    public ResonanceGenerator () {
    }

    public ResonanceGenerator (long timeout) {
	setTimeout (timeout);
    }

    public void setMolecule (Molecule m) {
	this.origmol = m;
	this.mol = m.cloneMolecule();
	atoms = mol.getAtomArray();
	bonds = mol.getBondArray();
	resonances.clear();
    }

    public void run () {
	if (mol == null) {
	    throw new IllegalStateException ("No input molecule specified!");
	}

	boolean neutral = true;
	for (int i = 0; i < bonds.length && neutral; ++i) {
	    MolBond b = bonds[i];
	    if (b.getAtom1().getCharge() != 0 
		|| b.getAtom2().getCharge() != 0) {
		neutral = false;
	    }
	}

	if (!neutral) {
	    start = System.currentTimeMillis();
	    generateResonances ();
	}
    }

    public void setTimeout (long timeout) { this.timeout = timeout; }
    public long getTimeout () { return timeout; }

    static protected boolean hasPiBond (MolAtom atom) {
	for (int i = 0; i < atom.getBondCount(); ++i) {
	    MolBond b = atom.getBond(i);
	    if (b.getType() == 2 || b.getType() == 3) {
		return true;
	    }
	}
	return false;
    }

    protected void generateResonances () {
	BitSet aset = new BitSet ();
	BitSet pi = new BitSet ();
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond b = bonds[i];
	    // we should also consider triple bond too?
	    if (b.getType() == 2) { 
		int a1 = mol.indexOf(b.getAtom1());
		int a2 = mol.indexOf(b.getAtom2());
		if (!aset.get(a1)) {
		    aset.set(a1);
		    atoms[a1].setCharge(atoms[a1].getCharge()-1);
		}
		if (!aset.get(a2)) {
		    aset.set(a2);
		    atoms[a2].setCharge(atoms[a2].getCharge()-1);
		}
		//b.setFlags(1, MolBond.TYPE_MASK);
		b.setType(1);
		pi.set(i);
	    }
	}

	if (debug) {
	    System.err.println("## number of resonance pi bonds: "  
			       + pi.cardinality());
	}

	/*
	// now generate a list of all possible candidate pi bonds that are 
	//  connected from these pi bonds
	BitSet bset = new BitSet ();
	{
	    BitSet visited = new BitSet ();
	    for (int i = pi.nextSetBit(0); i >= 0; i = pi.nextSetBit(i+1)) {
		getConnectedBonds (i, visited, bset, 0, 5);
	    }
	}

	if (debug) {
	    for (int i = bset.nextSetBit(0); i >= 0; 
		 i = bset.nextSetBit(i+1)) {
		MolAtom a1 = bonds[i].getAtom1();
		MolAtom a2 = bonds[i].getAtom2();
		System.err.println("## bond " + (mol.indexOf(a1)+1) + "-"
				   + (mol.indexOf(a2)+1) 
				   + " is candidate");
		_atom (System.err, mol.indexOf(a1)+1,a1);
		_atom (System.err, mol.indexOf(a2)+1,a2);
	    }
	    System.err.println("** " + bset.cardinality() 
			       + " candidate pi bonds!");
	}
	*/

	BitSet bset = new BitSet ();
	for (int i = 0; i < bonds.length; ++i) {
	    MolAtom a1 = bonds[i].getAtom1();
	    MolAtom a2 = bonds[i].getAtom2();

	    if (ElementData.checkValence(a1.getAtno(), a1.getCharge()+1,
					 a1.getValence()+1, 
					 a1.getImplicitHcount())
		&& ElementData.checkValence(a2.getAtno(), a2.getCharge()+1,
					    a2.getValence()+1,
					    a2.getImplicitHcount())) {
		if (debug) {
		    System.err.println("## bond " + (mol.indexOf(a1)+1) + "-"
				       + (mol.indexOf(a2)+1) 
				       + " is candidate");
		    _atom (System.err, mol.indexOf(a1)+1,a1);
		    _atom (System.err, mol.indexOf(a2)+1,a2);
		}
		bset.set(i);
	    }
	}

	if (debug) {
	    System.err.println("** " + bset.cardinality() 
			       + " candidate pi bonds!");
	}

	final int[] cands = new int[bset.cardinality()];
	for (int i = bset.nextSetBit(0), j = 0; 
	     i >= 0; i = bset.nextSetBit(i+1)) {
	    cands[j++] = i;
	}

	final Set<BitSet> configs = new HashSet<BitSet>();
	final int hcount = pi.cardinality();

	GrayCode gray = GrayCode.createBinaryGrayCode(cands.length);
	gray.addObserver(new Observer () {
		public void update (Observable o, Object arg) {
		    int[] code = (int[])arg;
		    BitSet atoms = new BitSet ();
		    BitSet c = new BitSet ();
		    for (int i = 0; i < code.length; ++i) {
			if (code[i] != 0) {
			    MolBond b = bonds[cands[i]];
			    int a1 = mol.indexOf(b.getAtom1());
			    int a2 = mol.indexOf(b.getAtom2());
			    if (atoms.get(a1) || atoms.get(a2)) {
				// bad configuration
			    }
			    else {
				atoms.set(a1);
				atoms.set(a2);
				c.set(cands[i]);
			    }
			}
		    }

		    // for now only consider resonance that have the same
		    //   hydrogen count as the original
		    if (c.cardinality() == hcount) {
			configs.add(c);
		    }

		    if (timeout > 0 
			&& (System.currentTimeMillis() - start > timeout)) {
			if (debug) {
			    System.err.println("## timeout (" +
					       timeout+") reached; "
					       +"search truncated!");
			}
			// if we've reached the timeout, then don't bother
			//  to process intermediate results
			configs.clear();
			o.deleteObserver(this);
		    }
		}
	    });
	gray.generate();

	// valid configuration..
	final Map<Molecule, Integer> scores = 
	    new IdentityHashMap<Molecule, Integer>();

	for (BitSet c : configs) {
	    Molecule m = mol.cloneMolecule();
	    MolBond[] bnds = m.getBondArray();
	    for (int i = c.nextSetBit(0); i >= 0; i = c.nextSetBit(i+1)) {
		MolAtom a1 = bnds[i].getAtom1();
		MolAtom a2 = bnds[i].getAtom2();
		a1.setCharge(a1.getCharge()+1);
		a2.setCharge(a2.getCharge()+1);
		bnds[i].setType(2);
	    }

	    int score = ChemUtil.calcMolScore(m);
	    scores.put(m, score);
	    resonances.add(m);

	    if (debug) {
		System.err.print("## valid configuration:");
		for (int i = c.nextSetBit(0); i >= 0; i = c.nextSetBit(i+1)) {
		    MolBond b = bonds[i];
		    System.err.print
			(" " + (mol.indexOf(b.getAtom1())+1)
			 +"-"+(mol.indexOf(b.getAtom2())+1));
		}
		System.err.println();
	    }
	}
	Collections.sort(resonances, new Comparator<Molecule>() {
			     public int compare (Molecule m1, Molecule m2) {
				 return scores.get(m2) - scores.get(m1);
			     }
			 });
	if (debug) {
	    for (Molecule m : resonances) {
		System.err.println(m.toFormat("smiles:q") 
				   + " " + scores.get(m));
	    }
	}
    }

    protected void getConnectedBonds 
	(int b, BitSet visited, BitSet connected, int depth, int maxdepth) {
	if (depth >= maxdepth) {
	    return;
	}

	visited.set(b);
	connected.set(b);
	System.err.println("depth: " + depth);

	MolAtom a1 = bonds[b].getAtom1();
	for (int i = 0; i < a1.getBondCount(); ++i) {
	    MolBond bnd = a1.getBond(i);
	    int xb = mol.indexOf(bnd);
	    if (!visited.get(xb) && xb != b) {
		MolAtom xa = bnd.getOtherAtom(a1);
		if (ElementData.checkValence(xa.getAtno(), xa.getCharge()+1,
					     xa.getValence()+1, 
					     xa.getImplicitHcount())) {
		    getConnectedBonds (xb, visited, connected, 
				       depth+1, maxdepth);
		}
		visited.set(xb);
	    }
	}

	// now do the other end
	MolAtom a2 = bonds[b].getAtom2();
	for (int i = 0; i < a2.getBondCount(); ++i) {
	    MolBond bnd = a2.getBond(i);
	    int xb = mol.indexOf(bnd);
	    if (!visited.get(xb) && xb != b) {
		MolAtom xa = bnd.getOtherAtom(a2);
		if (ElementData.checkValence(xa.getAtno(), xa.getCharge()+1,
					     xa.getValence()+1, 
					     xa.getImplicitHcount())) {
		    getConnectedBonds (xb, visited, connected, 
				       depth+1, maxdepth);
		}
		visited.set(xb);
	    }
	}
    }

    protected static void _atom (PrintStream ps, int index, MolAtom a) {
	ps.println("  " + String.format("%1$3d", index)
		   +"[a="+a.getAtno()
		   +",m="+a.getAtomMap()
		   +",c="+a.getCharge()
		   +",h="+a.getImplicitHcount() 
		   +",H="+a.getExplicitHcount()
		   +",r="+a.getRadical()
		   +",q="+a.getQuerystr()
		   +",l="+a.getQueryLabel()
		   +",x="+a.getExtraLabel()
		   +",v="+a.getValence()
		   +",s="+a.getSymbol()
		   +"]");
    }

    public Enumeration<Molecule> getResonances () {
	return resonances.elements();
    }

    // return canonical resonance
    public Molecule getCanonicalResonance () {
	return resonances.isEmpty() ? origmol : resonances.firstElement();
    }

    public int getResonanceCount () { return resonances.size(); }

    protected static void outputResonances (ResonanceGenerator resgen, 
					    java.io.PrintStream os, 
					    Molecule mol) 
	throws Exception {
	String name = mol.getName();
	if (name == null || name.equals("")) {
	    name = mol.getProperty("field_0");
	    mol.setName(name);
	}

	mol.dearomatize();
	resgen.setMolecule(mol);
	resgen.run();

	System.out.println(name + " has " 
			   + resgen.getResonanceCount() 
			   + " resonance(s)");
	int cnt = 1;
	for (Enumeration<Molecule> res = resgen.getResonances();
	     res.hasMoreElements(); ++cnt) {
	    Molecule t = res.nextElement();
	    int score = ChemUtil.calcMolScore(t);
	    System.out.println(t.toFormat("smiles:q") + " "+cnt
			       + " " + score);
	}
	System.out.println(resgen.getCanonicalResonance().toFormat("smiles:q")
			   + " Canonical" + " " + ChemUtil.calcMolScore
			   (resgen.getCanonicalResonance()));
    }

    public static void main (String argv[]) throws Exception {

	ResonanceGenerator resgen = new ResonanceGenerator (3000);
	if (argv.length == 0) {
	    System.err.println("** reading from STDIN...");

	    MolImporter mi = new MolImporter (System.in);
	    for (Molecule mol = new Molecule (); mi.read(mol); ) {
		outputResonances (resgen, System.out, mol);
	    }
	}
	else {
	    for (int i = 0; i < argv.length; ++i) {
		MolImporter mi = new MolImporter (argv[i]);
		for (Molecule mol = new Molecule (); mi.read(mol); ) {
		    outputResonances (resgen, System.out, mol);
		}
		mi.close();
	    }
	}
    }
}
