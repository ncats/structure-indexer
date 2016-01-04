package tripod.chem.indexer.util;

import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;


/**
 * A simple ring perception implementation based on 
 * J.D. Horton, A polynomial-time algorithm to find the shortest cycle
 *   basis of a graph. SIAM J. Comput., 16(2):358-366, 1987.
 */
public class Ring {
    private static final Logger logger = 
	Logger.getLogger(Ring.class.getName());

    static class Path {
        BitSet aset = new BitSet ();
        BitSet bset = new BitSet ();
        int[] atoms, bonds;

        Path (Collection<MolBond> path) {
            bonds = new int[path.size()];
            int i = 0;
            for (MolBond b : path) {
                Molecule mol = (Molecule)b.getParent();
                aset.set(mol.indexOf(b.getAtom1()));
                aset.set(mol.indexOf(b.getAtom2()));
                int idx = mol.indexOf(b);
                bset.set(idx);
                bonds[i++] = idx;
            }

            atoms = new int[aset.cardinality()];
            i = 0;
            for (int j = aset.nextSetBit(0); j >= 0; j = aset.nextSetBit(j+1))
                atoms[i++] = j;
        }

        public BitSet getAtomSet () { return aset; }
        public BitSet getBondSet () { return bset; }
        public int[] getAtoms () { return atoms; }
        public int[] getBonds () { return bonds; }
    }

    private Molecule mol;
    private int[][] cost;
    private int[][] btab;
    private MolAtom[] atoms;

    private int[][] rings;
    private Set<BitSet> ringsets = new HashSet<BitSet>();
    
    public Ring (Molecule mol) {
        this (mol, 0);
    }

    public Ring (Molecule mol, int maxsize) {
        this.mol = mol;
        atoms = mol.getAtomArray();

        int acount = atoms.length;
        cost = new int[acount][acount];
        btab = mol.getBtab();

        // init the path cost...
        for (int i = 0; i < acount; ++i) {
            cost[i][i] = 0;
            for (int j = i+1; j < acount; ++j) {
                cost[i][j] = cost[j][i] = btab[i][j] < 0 ? acount : 1;
            }
        }
        
        // 1. now perform floyd-warshall's all pairs shortest path
        for (int k = 0; k < acount; ++k) 
            for (int i = 0; i < acount; ++i) 
                for (int j = 0; j < acount; ++j) 
                    cost[i][j] = Math.min
                        (cost[i][j], cost[i][k]+cost[k][j]);

        // 2. find all cycles of at most <= maxsize
        MolBond[] bonds = mol.getBondArray();
        for (int i = 0; i < atoms.length; ++i) {
            for (int j = 0; j < bonds.length; ++j) {
                int a1 = mol.indexOf(bonds[j].getAtom1());
                int a2 = mol.indexOf(bonds[j].getAtom2());
                BitSet[] path1 = getPaths (i, a1);
                BitSet[] path2 = getPaths (i, a2);
                for (BitSet p1 : path1) {
                    for (BitSet p2 : path2) {
                        if (p1.intersects(p2)) {
                            BitSet bs = (BitSet)p1.clone();
                            bs.and(p2);
                            if (bs.cardinality() == 1 
                                && bs.nextSetBit(0) == i) {
                                bs.clear();
                                bs.or(p1);
                                bs.or(p2);
                                if (maxsize <= 2 
                                    || bs.cardinality() <= maxsize) {
                                    ringsets.add(bs);
                                }
                            }
                        }
                    }
                }
            }
        }

        /*
        logger.info(" "+urings.size()+" unique rings!");
        for (BitSet r : urings) {
            System.err.println(" "+r.cardinality()+" "+r);
        }
        */

        rings = new int[ringsets.size()][];
        int i = 0;
        for (BitSet bs : ringsets) {
            int[] r = new int[bs.cardinality()];
            for (int k = bs.nextSetBit(0), j = 0; 
                 k >= 0; k = bs.nextSetBit(k+1))
                r[j++] = k;
            rings[i++] = r;
        }
    }

    /**
     * Return all rings of sizes <= maxsize
     */
    public int[][] getRings () { return rings; }
    public int getRingCount () { return rings.length; }
    public BitSet[] getRingSets () { 
        return ringsets.toArray(new BitSet[0]); 
    }
    // return the number of rings that contains this atom
    public int getRingCount (int atom) { 
        int c = 0;
        for (BitSet r : ringsets) {
            if (r.get(atom)) ++c;
        }
        return c;
    }

    /**
     * Return all shortest paths between atoms start and end
     */
    private BitSet[] getPaths (int start, int end) {
        List<BitSet> paths = new ArrayList<BitSet>();
        BitSet visited = new BitSet (atoms.length);
        Stack<MolBond> path = new Stack<MolBond>();
        dfs (paths, new Stack<MolBond>(), visited, start, start, end);
        return paths.toArray(new BitSet[0]);
    }
 
    private void dfs (List<BitSet> paths, Stack<MolBond> path, 
                      BitSet visited, int start, int a, int end) {
        if (a == end) {
            BitSet p = new BitSet (atoms.length);
            for (MolBond b : path) {
                p.set(mol.indexOf(b.getAtom1()));
                p.set(mol.indexOf(b.getAtom2()));
            }
            paths.add(p);

            return;
        }
        
        visited.set(a);
        MolAtom atom = atoms[a];
        for (int b = 0; b < atom.getBondCount(); ++b) {
            MolBond bond = atom.getBond(b);
            int xa = mol.indexOf(bond.getOtherAtom(atom));

            if ((cost[start][xa] + cost[xa][end] <= cost[start][end])
                && !visited.get(xa)) {
                path.push(bond);
                dfs (paths, path, visited, start, xa, end);
                path.pop();
            }
        }
        visited.clear(a);
    }

    public static void main (String[] argv) throws Exception {
        for (String a : argv) {
            MolHandler mh = new MolHandler (a);
            int[][] rings = new Ring (mh.getMolecule()).getRings();
            System.out.println("## "+a);
            System.out.println("..."+rings.length+" rings!");
            for (int[] r : rings) {
                System.out.print(String.format("%1$3d", r.length)+": {"+r[0]);
                for (int i = 1; i < r.length; ++i) {
                    System.out.print(","+r[i]);
                }
                System.out.println("}");
            }
        }
    }
}