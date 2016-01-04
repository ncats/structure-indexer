// $Id$

package tripod.chem.indexer.util;

import java.io.InputStream;
import java.io.FileInputStream;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;


public class MolPath2 implements Comparator<String> {
    private static final Logger logger = 
	Logger.getLogger(MolPath.class.getName());

    public interface PathVisitor {
	void visit (Molecule mol, List<MolBond> path);
    }

    public interface Annotator {
	String getLabel (MolAtom a);
	String getLabel (MolBond b);
    }

    public class Path implements Comparable<Path> {
	MolBond[] bPath;
	MolAtom[] aPath;
	String path;
	int[] aseq;

	public int hashCode () { return path.hashCode(); }
	public int compareTo (Path p) { return path.compareTo(p.path); }
	public boolean equals (Object obj) {
	    if (obj == this) {
		return true;
	    }
	    if (obj instanceof Path) {
		return ((Path)obj).compareTo(this) == 0;
	    }
	    return false;
	}

	void setPath (MolAtom[] aPath) {
	    aseq = new int[aPath.length];
	    for (int i = 0; i < aPath.length; ++i) {
		aseq[i] = mol.indexOf(aPath[i]);
	    }
	    this.aPath = aPath;
	}
	void setPath (MolBond[] bPath) {
	    this.bPath = bPath;
	}
	void setPath (String path) {
	    this.path = path;
	}

	public int[] getAtoms () { return aseq; }

	public MolPath2 getMolPath () { return MolPath2.this; }
	public MolBond[] getBondPath () { return bPath; }
	public MolAtom[] getAtomPath () { return aPath; }
	public int getBondCount () { return bPath.length; }
	public int getAtomCount () { return aPath.length; }
	public String toString () { return path; }
    }

    public static class AllPathVisitor implements PathVisitor {
	ArrayList<MolBond[]> paths = new ArrayList<MolBond[]>();

	public void visit (Molecule m, List<MolBond> path) {
	    paths.add(path.toArray(new MolBond[0]));
	}

	public MolBond[][] getPaths () {
	    return paths.toArray(new MolBond[0][]);
	}
    }

    private Molecule mol;
    private MolAtom[] atoms;
    private MolBond[] bonds;
    private Path[][][] paths;
    // lookup of paths by its canonical string
    private Map<String, Path[]> pathLUT;
    private Map<String, Integer> sv; // sparse vector of path counts

    public MolPath2 () {
    }

    public MolPath2 (Molecule mol) {
	setMolecule (mol);
    }

    public int compare (String s1, String s2) {
	int d = s2.length()  - s1.length();
	if (d == 0) {
	    d = s1.compareTo(s2);
	}
	return d;
    }

    public void setMolecule (Molecule mol) {
	this.mol = mol;
	atoms = this.mol.getAtomArray();
	bonds = this.mol.getBondArray();
	paths = new Path[atoms.length][atoms.length][];

	Map<String, List<Path>> lut = new TreeMap<String, List<Path>> (this);
	
	int size = atoms.length;
	for (int i = 0; i < size; ++i) {
	    for (int j = i+1; j < size; ++j) {
		Path[] paths = getPaths (i, j);
		for (int k = 0; k < paths.length; ++k) {
		    String s = paths[k].toString();
		    List<Path> lv = lut.get(s);
		    if (lv == null) {
			lut.put(s, lv = new ArrayList<Path>());
		    }
		    lv.add(paths[k]);
		}
	    }
	}
	
	pathLUT = new TreeMap<String, Path[]>(this);
	sv = new TreeMap<String, Integer>(this);
	
	for (Map.Entry<String, List<Path>> e : lut.entrySet()) {
	    List<Path> lv = e.getValue();
	    pathLUT.put(e.getKey(), lv.toArray(new Path[0]));
	    sv.put(e.getKey(), lv.size());
	}
    }

    public Molecule getMolecule () { return mol; }
    public MolAtom[] getAtoms () { return atoms; }
    public MolBond[] getBonds () { return bonds; }
    public int getAtomCount () { return atoms != null ? atoms.length : -1; }
    public int getBondCount () { return bonds != null ? bonds.length : -1; }
    public MolAtom getAtom (int i) { return atoms[i]; }
    public MolBond getBond (int j) { return bonds[j]; }

    /*
     * Given a starting and ending atoms, this method returns
     * all unique paths that connect the two atoms
     */
    public void walkPaths (PathVisitor visitor, 
			   MolAtom startAtom, MolAtom endAtom) {
	if (mol == null) {
	    throw new IllegalArgumentException ("No input molecule defined");
	}
	int start = mol.indexOf(startAtom);
	int end = mol.indexOf(endAtom);

	if (start < 0) {
	    throw new IllegalArgumentException
		("Starting atom is not part of molecule");
	}

	if (end < 0) {
	    throw new IllegalArgumentException
		("Ending atom is not part of molecule");
	}

	walkPaths (visitor, start, end);
    }


    public void walkPaths (PathVisitor visitor, int start, int end) {
	walkPaths (visitor, -1, start, end);
    }

    /*
     * length <= atom count
     */
    public void walkPaths (PathVisitor visitor, int length, 
			   int start, int end) {
	if (visitor == null) {
	    throw new IllegalArgumentException ("Path visitor is null");
	}

	if (start < 0 || start >= atoms.length) {
	    throw new IllegalArgumentException
		("Invalid starting atom specified");
	}

	if (end < 0 || end >= atoms.length) {
	    throw new IllegalArgumentException
		("Invalid ending atom specified");
	}

	BitSet visited = new BitSet (atoms.length);

	Stack<MolBond> path = new Stack<MolBond>();
	dfs (visitor, path, visited, length, start, end);
    }

    public Map<String, Integer> getSparseVector () { return sv; }

    /*
     * Generate a fingerprint where numInts is the length of the 
     *  fingerprint in 32-bit units.  That is, if an 1024-bit fingerprint
     *  is desired, then numInts is 32.
     */
    public int[] generateFingerprint (int numInts) {
	return generateFingerprint (new int[numInts], 0, numInts);
    }

    public int[] generateFingerprint (int[] fp, int offset, int length) {
	for (int i = 0; i < length; ++i) {
	    fp[offset+i] = 0;
	}

	int nbits = length * 32;
	Random rand = new Random ();
	for (String s : sv.keySet()) {
	    long hash = (long)s.hashCode() & 0xffffffffl;
	    int bit =  (int)(hash % nbits); 
	    //System.out.println(s + " => " + bit);
	    fp[offset+bit/32] |= 1 << (31 - bit%32);
	    if (s.indexOf('.') > 0) {
		// multiple path, then turn on additional bit
		rand.setSeed(hash);
		bit = rand.nextInt(nbits);
		fp[offset+bit/32] |= 1 << (31 - bit%32);
	    }
	}
	return fp;
    }


    /*
     * Default to return all paths
     */
    public Path[] getPaths (int a, int b) {
	Path[] p = paths[a][b];
	if (p == null) {
	    AllPathVisitor visitor = new AllPathVisitor ();
	    walkPaths (visitor, a, b);

	    MolBond[][] path = visitor.getPaths();
	    p = new Path[path.length];
	    //System.out.println("Paths: "+(a+1)+" "+(b+1));
	    for (int i = 0; i < path.length; ++i) {
		p[i] = createPath (path[i]);

		/*
		int[] q = p[i].getAtoms();
		System.out.print(" "+String.format("%1$2d", q[0]+1));
		for (int j = 1; j < q.length; ++j) {
		    System.out.print(" " + String.format("%1$2d", q[j]+1));
		}
		System.out.println(" "+p[i]);
		*/
	    }
	    paths[a][b] = paths[b][a] = p;
	}
	return p;
    }

    public Map<String, Path[]> getPaths () { 
	// return all paths
	return pathLUT;
    }

    protected static String getDefaultLabel (MolAtom atom) {
	return atom.getSymbol() + (atom.hasAromaticBond() ? "ar" 
				   : "sp"+atom.getHybridizationState());
    }

    private Annotator pathAnnotator = createDefaultAnnotator ();
    public Annotator getAnnotator () { return pathAnnotator; }
    public void setAnnotator (Annotator annotator) { 
	this.pathAnnotator = annotator; 
    }

    public static Annotator createDefaultAnnotator () {
	return new Annotator () {
		public String getLabel (MolAtom a) {
		    return MolPath.getDefaultLabel(a);
		}
		public String getLabel (MolBond b) {
		    return "";
		}
	    };
    }

    Path createPath (MolBond[] path) {
	if (path == null || path.length == 0) {
	    return null;
	}

	Annotator anno = getAnnotator ();

	if (path.length == 1) {
	    String s1 = anno.getLabel(path[0].getAtom1());
	    String s2 = anno.getLabel(path[0].getAtom2());
	    Path p = new Path ();
	    p.setPath(path);
	    if (s1.compareTo(s2) <= 0) {
		p.setPath(s1+anno.getLabel(path[0])+s2);
		p.setPath(new MolAtom[]{
			      path[0].getAtom1(), path[0].getAtom2()
			  });
	    }
	    else {
		p.setPath(s2+anno.getLabel(path[0])+s1);
		p.setPath(new MolAtom[]{
			      path[0].getAtom2(), path[0].getAtom1()
			  });
	    }
	    return p;
	}

	MolAtom head = null; // figuring out direction
	if (path[0].getAtom1() == path[1].getAtom1()
	    || path[0].getAtom1() == path[1].getAtom2()) {
	    head = path[0].getAtom2();
	}
	else if (path[0].getAtom2() == path[1].getAtom1()
		 || path[0].getAtom2() == path[1].getAtom2()) {
	    head = path[0].getAtom1();
	}
	else {
	    throw new IllegalArgumentException ("Path is disconnected");
	}

	MolAtom tail = null;
	{ int i = path.length - 1;
	    if (path[i].getAtom1() == path[i-1].getAtom1()
		|| path[i].getAtom1() == path[i-1].getAtom2()) {
		tail = path[i].getAtom2();
	    }
	    else if (path[i].getAtom2() == path[i-1].getAtom1()
		     || path[i].getAtom2() == path[i-1].getAtom2()) {
		tail = path[i].getAtom1();
	    }
	    else {
		throw new IllegalArgumentException ("Path is disconnected");
	    }
	}

	StringBuilder forward = new StringBuilder ();
	StringBuilder backward = new StringBuilder ();

	forward.append(anno.getLabel(head));
	backward.append(anno.getLabel(tail));
	int dir = head.getAtno() - tail.getAtno();

	MolBond[] rpath = new MolBond[path.length];
	MolAtom[] fwd = new MolAtom[path.length+1];
	MolAtom[] rev = new MolAtom[path.length+1];

	fwd[0] = head;
	rev[0] = tail;
	for (int i = 0, j = path.length-1; i < path.length; ++i, --j) {
	    head = path[i].getOtherAtom(head);
	    if (head == null) {
		throw new IllegalArgumentException ("Path is disconnected");
	    }
	    forward.append(anno.getLabel(path[i]));
	    forward.append(anno.getLabel(head));

	    tail = path[j].getOtherAtom(tail);
	    if (tail == null) {
		throw new IllegalArgumentException ("Path is disconnected");
	    }
	    backward.append(anno.getLabel(path[j]));
	    backward.append(anno.getLabel(tail));

	    if (dir == 0) {
		dir = head.getAtno() - tail.getAtno();
	    }
	    rpath[i] = path[j];
	    fwd[i+1] = head;
	    rev[i+1] = tail;
	}

	Path p = new Path ();
	if (dir <= 0) {
	    p.setPath(forward.toString());
	    p.setPath(path);
	    p.setPath(fwd);
	}
	else {
	    p.setPath(backward.toString());
	    p.setPath(rpath);
	    p.setPath(rev);
	}
	return p;
    }

    protected void dfs (PathVisitor visitor, Stack<MolBond> path, 
			BitSet visited, int length, int a, int end) {
	/*
	if (visited.get(a)) {
	    return;
	}
	*/

	if (a == end) {
	    visitor.visit(getMolecule (), path);
	    return;
	}

	if (length < 0 || path.size() <= length) {
	    visited.set(a);
	    MolAtom atom = atoms[a];
	    for (int b = 0; b < atom.getBondCount(); ++b) {
		MolBond bond = atom.getBond(b);
		int xa = mol.indexOf(bond.getOtherAtom(atom));
		if (!visited.get(xa)) {
		    path.push(bond);
		    dfs (visitor, path, visited, length, xa, end);
		    path.pop();
		}
	    }
	    visited.clear(a);
	}
    }
    
    static void genpath (Molecule mol) {
	mol.aromatize();
	mol.calcHybridization();
	MolPath2 mp2 = new MolPath2 (mol);
	
	Map<String, Path[]> paths = mp2.getPaths();
	System.out.println("** " + mol.getName() + " has "+paths.size() 
			   + " unique paths!");
	for (Map.Entry<String, Path[]> e : paths.entrySet()) {
	    System.out.println(e.getKey());
	    for (Path p : e.getValue()) {
		int[] a = p.getAtoms();
		System.out.print("  "+String.format("%1$2d", a[0]+1));
		for (int i = 1; i < a.length; ++i) {
		    System.out.print(" " + String.format("%1$2d", a[i]+1));
		}
		System.out.println();
	    }
	}
    }

    public static void main (String[] argv) throws Exception {
	for (String f : argv) {
	    MolImporter mi = new MolImporter (f);
	    for (Molecule mol = new Molecule (); mi.read(mol); ) {
		genpath (mol);
	    }
	    mi.close();
	}
    }
}
