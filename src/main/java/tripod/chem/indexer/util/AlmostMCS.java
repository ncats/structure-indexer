// $Id: AlmostMCS.java 3504 2009-10-29 16:02:45Z nguyenda $

package tripod.chem.indexer.util;

import java.util.*;
import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.*;

public class AlmostMCS {
    static final int EMPTY = -1;
    static final int NOT_MATCHED = -2;

    public interface AtomComparator {
	public boolean match (MolAtom a, MolAtom b);
    }

    public interface BondComparator {
	public boolean match (MolBond a, MolBond b);
    }

    static class DefaultAtomComparator implements AtomComparator {
	public boolean match (MolAtom a, MolAtom b) {
	    int atnoA = a.getAtno(), atnoB = b.getAtno();
	    
	    return (atnoA == MolAtom.ANY 
		    || atnoB == MolAtom.ANY 
		    || ((atnoA == atnoB)
			&& (a.getHybridizationState() 
			    == b.getHybridizationState())));
	}
    }

    static class DefaultBondComparator implements BondComparator {
	public boolean match (MolBond a, MolBond b) {
	    int atype = a.getType(), btype = b.getType();
	    return atype == btype
		|| atype == MolBond.ANY
		|| btype == MolBond.ANY;
	}
    }

    private int [][]M = {}; // vertex matching 
    private Molecule query, target;
    private int[][] qCtab, tCtab;
    private int[][] qBtab, tBtab;
    private Random rand = new Random ();

    private int[] resultF, resultR;
    private Molecule resultCore;
    private int size;

    private AtomComparator atomComparator = new DefaultAtomComparator ();
    private BondComparator bondComparator = new DefaultBondComparator ();

    public AlmostMCS () {
    }

    public void setAtomComparator (AtomComparator comparator) {
	atomComparator = comparator;
    }
    public void setBondComparator (BondComparator comparator) {
	bondComparator = comparator;
    }

    public void setQuery (Molecule query) {
	this.query = query;
	qCtab = query.getCtab();
	qBtab = query.getBtab();
    }
    public Molecule getQuery () { return query; }
    public void setTarget (Molecule target) {
	this.target = target; 
	tCtab = target.getCtab();
	tBtab = target.getBtab();
    }
    public Molecule getTarget () { return target; }
    public void setMolecules (Molecule query, Molecule target) {
	setQuery (query);
	setTarget (target);
    }

    public boolean search () {
	if (query == null) {
	    throw new IllegalStateException ("No query molecule specified");
	}
	if (target == null) {
	    throw new IllegalStateException ("No target molecule specified");
	}

	MolAtom []qAtoms = query.getAtomArray();
	MolAtom []tAtoms = target.getAtomArray();

	M = new int[qAtoms.length][];
	for (int i = 0; i < qAtoms.length; ++i) {
	    Vector<Integer> m = new Vector<Integer>();
	    for (int j = 0; j < tAtoms.length; ++j) {
		if (atomComparator.match(qAtoms[i], tAtoms[j])) {
		    m.add(j);
		}
	    }
	    M[i] = new int[m.size()];
	    for (int j = 0; j < M[i].length; ++j) {
		M[i][j] = m.get(j);
	    }
	}

	int[] fwd = new int[qAtoms.length];
	for (int i = 0; i < fwd.length; ++i) {
	    fwd[i] = EMPTY;
	}
	int[] rev = new int[tAtoms.length];
	for (int i = 0; i < rev.length; ++i) {
	    rev[i] = EMPTY;
	}

	size = 0;
	resultF = new int[0];
	resultR = new int[0];
	resultCore = null;

	int score = 0;
	for (int v = 0; v < M.length; ++v) {
	    for (int j = 0; j < M[v].length; ++j) {
		int w = M[v][j];
		
		// mapping of query atoms to target's
		int[] F = new int[qAtoms.length];
		System.arraycopy(fwd, 0, F, 0, F.length);
		    
		// mapping of target atoms to query's
		int[] R = new int[tAtoms.length];
		System.arraycopy(rev, 0, R, 0, R.length);
		    
		F[v] = w;
		R[w] = v;
		    
		extendMatch (F, R);

		int s = 0;
		for (int k = 0; k < F.length; ++k) {
		    if (F[k] >= 0) {
			//System.out.print(" " + (k+1) + "<->" + (F[k]+1));
			++s;
		    }
		}
		//System.out.println();
		if (s >= size) {
		    // score this 
		    int myScore = 
			scoreMCS (createMolecule (query, F, false))
			+ scoreMCS (createMolecule (target, R, false));
		    if (myScore > score) {
			size = s;
			resultF = F;
			resultR = R;
			score = myScore;
		    }
		}

		if (s == qAtoms.length || s == tAtoms.length) {
		    // we're done...bail out early
		    return true;
		}

		F = null;
		R = null;
	    }
	}

	fwd = null;
	rev = null;

	return size > 0;
    }

    // simply sum up the number of rings in the core
    private int scoreMCS (Molecule core) {
	int[][][] rings = core.getAromaticAndAliphaticRings
	    (Molecule.AROM_BASIC, false, false, 0, 0);
	return rings[0].length + rings[1].length;
    }

    // return true if there exists a node that has been matched but one
    //  of its neighbor hasn't.  return false if there is no such candidate
    //  left.
    private boolean getNextCandidate (int edge[], int[] fwd) {
	Vector<Integer> best = new Vector<Integer>();

	for (int i = 0; i < fwd.length; ++i) {
	    // get all neighbors of this atom
	    if (fwd[i] >= 0) {
		Vector<Integer> cand = new Vector<Integer>();
		int nb[] = qCtab[i];
		for (int j = 0; j < nb.length; ++j) {
		    int h = nb[j];
		    // if this neighbor hasn't been matched and that there
		    //   exist at least one vertex matching
		    if (fwd[h] == EMPTY && M[h].length > 0) {
			cand.add(h);
		    }
		}

		if (best.size() < cand.size()) {
		    best = cand;
		    edge[0] = i;
		}
	    }
	}

	boolean hasNext = !best.isEmpty();
	if (hasNext) {
	    //Collections.shuffle(best);
	    edge[1] = best.firstElement();

	    for (Integer nb : best) {
		if (M[nb].length < M[edge[1]].length) {
		    edge[1] = nb;
		}
	    }
	}
	return hasNext;
    }

    private static MolBond getBond (Molecule mol, int atom1, int atom2) {
	MolAtom a = mol.getAtom(atom1);
	MolAtom b = mol.getAtom(atom2);

	int nb = a.getBondCount();
	for (int i = 0; i < nb; ++i) {
	    MolBond bond = a.getBond(i);
	    if (bond.getOtherAtom(a) == b) {
		return bond;
	    }
	}
	return null;
    }

    private int extendMatch (int[] fwd, int[] rev) {
	int size = 0; // matching bond count

	int edge[] = new int[2];
	while (getNextCandidate (edge, fwd)) {
	    int sv = edge[0], v = edge[1];

	    //MolBond qbond = getBond (query, sv, v);
	    MolBond qbond = query.getBond(qBtab[sv][v]);

	    int np = fwd[sv];
	    int[] nb = qCtab[v];

	    Vector<Integer> cand = new Vector<Integer>();
	    for (int i = 0; i < M[v].length; ++i) {
		int vp = M[v][i];

		if (rev[vp] == EMPTY) {
		    //MolBond tbond = getBond (target, vp, np);
		    int index = tBtab[vp][np];
		    
		    if (index >= 0 && bondComparator.match
			(qbond, target.getBond(index))) {

			// now check all neighbors of tnv to make sure 
			//  that the inclusion of this match is 
			//  consistent

			int p = 0, q = 0;
			for (int k = 0; k < nb.length; ++k) {
			    int h = fwd[nb[k]];
			    if (h >= 0) {
				//if (getBond (target, vp, h) != null) {
				if (tBtab[vp][h] >= 0) {
				    ++p;
				}
				++q;
			    }
			}

			if (p == q) {
			    cand.add(vp);
			}
		    }
		}
	    }

	    if (cand.size() > 0) {
		++size;
		// simply picks an edge (randomly or otherwise)
		//if (cand.size() > 1) {
		    //Collections.shuffle(cand);
		//}

		int n = cand.firstElement();
		fwd[v] = n;
		rev[n] = v;
	    }
	    else { // no match
		fwd[v] = NOT_MATCHED;
	    }
	}
	return size;
    }

    public int[] getResult () { return resultF; }
    public int getResultSize () { return size; }

    public Molecule getResultAsMolecule () {
	return getResultAsMolecule (true);
    }

    public Molecule getResultAsMolecule (boolean atomMapping) {
	if (resultF == null || resultF.length == 0) {
	    return null;
	}
	if (resultCore == null) {
	    resultCore = createMolecule (query, resultF, atomMapping);
	}
	return resultCore;
    }

    protected static Molecule createMolecule 
	(Molecule ref, int[] matches, boolean atomMapping) {

	Molecule mol = new Molecule ();
	MolBond[] bonds = ref.getBondArray();
	MolAtom[] atoms = new MolAtom[matches.length];
	for (int i = 0; i < bonds.length; ++i) {
	    MolBond bond = bonds[i];
	    int idx1 = ref.indexOf(bond.getAtom1());
	    int idx2 = ref.indexOf(bond.getAtom2());
	    if (matches[idx1] >= 0 && matches[idx2] >= 0) {
		// this bond is in the result structure
		MolAtom a1 = atoms[idx1];
		if (a1 == null) {
		    atoms[idx1] = a1 = (MolAtom)bond.getAtom1().clone();
		    a1.setImplicitHcount(0);
		    a1.setRadical(0);
		    if (atomMapping) {
			//a1.setAtomMap(idx1+1);
			a1.setAtomMap(matches[idx1]+1);
		    }
		    mol.add(a1);
		}
		MolAtom a2 = atoms[idx2];
		if (a2 == null) {
		    atoms[idx2] = a2 = (MolAtom)bond.getAtom2().clone();
		    a2.setImplicitHcount(0);
		    a2.setRadical(0);
		    if (atomMapping) {
			//a2.setAtomMap(idx2+1);
			a2.setAtomMap(matches[idx2]+1);
		    }
		    mol.add(a2);
		}
		MolBond b = bond.cloneBond(a1, a2);

		b.setFlags(0);
		// this forces the valences of a1 & a2 to be checked...
		b.setType(bond.getType());

		mol.add(b);
	    }
	}
	atoms = null;
	bonds = null;

	mol.valenceCheck();
	if (mol.getAtomCount() > 0) {
	    try {
		String smiles = mol.toFormat("smiles:u0-H");
		mol.clear();
		MolImporter.importMol(smiles, mol);
	    }
	    catch (MolFormatException ex) {
		ex.printStackTrace();
	    }
	}

	return mol;
    }


    public String getResultAsSmiles (boolean unique, boolean atomMap) {
	Molecule mol = getResultAsMolecule (atomMap);
	if (mol != null) {
	    try {
		return mol.toFormat("smiles" + (unique?":a_basu":""));
	    }
	    catch (Exception ex) {
		System.err.println("** query="+query.getName() + " target="
				   + target.getName());
		ex.printStackTrace();
	    }
	}
	return null;
    }

    public String getResultAsSmiles () { 
	return getResultAsSmiles (true, false);  
    }

    public static void main (String argv[]) throws Exception {
	if (argv.length == 0) {
	    System.out.println("usage: AlmostMCS FILES...");
	    System.exit(1);
	}

	Vector<Molecule> mols = new Vector<Molecule>();
	for (int i = 0; i < argv.length; ++i) {
	    MolImporter mi = new MolImporter (argv[i]);
	    for (Molecule m; (m = mi.read()) != null
		     /*&& mols.size() < 2000*/; ) {
		m.aromatize(Molecule.AROM_BASIC);
		m.calcHybridization();
		m.hydrogenize(false);
		String name = m.getProperty("field_0");
		if (name != null) {
		    m.setName(name);
		}
		mols.add(m);
	    }
	    mi.close();
	}
	//System.out.println(mols.size() + " molecules read!");
	
	AlmostMCS amcs = new AlmostMCS ();

	MCS mcs = new MCS ();
	mcs.setFastSearch(false);
	mcs.setMinimumCommonSize(2);

	for (int i = 0; i < mols.size(); ++i) {
	    Molecule query = mols.get(i);
	    amcs.setQuery(query);
	    mcs.setQuery(query);

	    for (int j = i; j < mols.size(); ++j) {
		Molecule target = mols.get(j);
		amcs.setTarget(target);
		mcs.setTarget(target);
		
		long start = System.currentTimeMillis();
		amcs.search();
		int[] result = amcs.getResult(); 
		long end = System.currentTimeMillis();

		System.out.print
		    (amcs.getResultAsSmiles(true, false) + "\tAlmostMCS");
		System.out.print
		    ("\t" + query.getName() + " vs. " + target.getName());
		System.out.print("\t" + amcs.getResultSize());
		for (int k = 0; k < result.length; ++k) {
		    if (result[k] >= 0) {
			System.out.print
			    (" " + (k+1) + "-" + (result[k]+1));
		    }
		}
		System.out.println();

		start = System.currentTimeMillis();
		boolean found = mcs.findFirst();
		end = System.currentTimeMillis();
		long total = end - start;

		if (found) {
		    result = mcs.getResult();
		    System.out.print(mcs.getResultAsSmiles(true)+"\tMCS");
		    System.out.print
			("\t" + query.getName() + " vs. " + target.getName());
		    System.out.print("\t" + mcs.getResultSize());
		    if (result != null) {
			for (int k = 0; k < result.length; ++k) {
			    if (result[k] >= 0) {
				System.out.print
				    (" " + (k+1) + "-" + (result[k]+1));
			    }
			}
			System.out.println();
		    }

		    start = System.currentTimeMillis();
		    found = mcs.findNext();
		    end = System.currentTimeMillis();
		    total += end-start;
		}
		/*
		System.out.println("... in " + total + "ms");
		System.out.println
		    ("** delta: " + (amcs.getResultSize() - mcs.getResultSize()));
		*/
	    }
	}
    }
}

