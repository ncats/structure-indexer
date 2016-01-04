// $Id: MCSMatrix.java 2278 2008-05-29 22:27:45Z nguyenda $

package tripod.chem.indexer.util;

import java.util.Map;
import java.util.HashMap;
import java.util.Vector;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;

public class MCSMatrix {
    private Molecule matrix[][];

    protected MCSMatrix (Molecule[] mols) {
	matrix = new Molecule[mols.length][];

	AlmostMCS mcs = new AlmostMCS ();
	//mcs.setRingBreaking(false);

	Map<String, Molecule> uniq = new HashMap<String, Molecule>();
	for (int i = 0; i < mols.length; ++i) {
	    Molecule mi = mols[i];
	    if (mi != null) {
		mcs.setQuery(mi);
		matrix[i] = new Molecule[mols.length-i];
		for (int j = i + 1; j < mols.length; ++j) {
		    Molecule mj = mols[j];
		    if (mj != null) {
			mcs.setTarget(mj);
			if (mcs.search()) {
			    Molecule newcore = mcs.getResultAsMolecule(false);
			    String smi = newcore.toFormat("smiles:ua_bas");
			    Molecule core = uniq.get(smi);
			    if (core == null) {
				uniq.put(smi, newcore);
				core = newcore;
			    }
			    matrix[i][j-i] = core;
			    newcore = null;
			}
		    }
		}
		matrix[i][0] = mi;
	    }
	}
	uniq.clear();
    }

    public Molecule get (int i, int j) {
	if (i > j) {
	    int t = j;
	    j = i;
	    i = t;
	}
	return matrix[i][j-i];
    }

    public static MCSMatrix createMatrix (Molecule[] mols) {
	return new MCSMatrix (mols);
    }

    public static void main (String[] argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("usage: MCSMatrix FILES...");
	    System.exit(1);
	}

	Vector<Molecule> mols = new Vector<Molecule>();
	for (int i = 0; i < argv.length; ++i) {
	    MolImporter mi = new MolImporter (argv[i]);
	    for (Molecule m; (m = mi.read()) != null; ) {
		m.aromatize(Molecule.AROM_BASIC);
		m.calcHybridization();

		String name = m.getName();
		// use the first non-empty property 
		for (int j = 0; j < m.getPropertyCount() 
			 && (name == null || name.equals("")); ++j) {
		    name = m.getProperty(m.getPropertyKey(j));
		}
		m.setName(name);
		mols.add(m);
	    }
	}

	Molecule[] array = mols.toArray(new Molecule[0]);
	MCSMatrix matrix = createMatrix (array);
	for (int i = 0; i < array.length; ++i) {
	    for (int j = 0; j < array.length; ++j) {
		Molecule c= matrix.get(i, j);
		if (c != null) {
		    System.out.println(c.toFormat("smiles:ua_bas") 
				       + "\t("+i+","+j+") => (" 
				       + array[i].getName() + "," 
				       + array[j].getName() + ")");
		}
	    }
	}
	
    }
}
