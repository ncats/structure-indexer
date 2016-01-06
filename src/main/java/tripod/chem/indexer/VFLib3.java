package tripod.chem.indexer;

import java.util.Map;
import java.util.HashMap;
import java.util.BitSet;
import java.util.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import gov.nih.ncgc.v3.api.Atom;
import gov.nih.ncgc.v3.api.Bond;
import gov.nih.ncgc.v3.api.Bond.BondType;
import gov.nih.ncgc.v3.api.BondTable;
import gov.nih.ncgc.v3.api.Chemical2;
import gov.nih.ncgc.v3.api.GraphInvariant;


/*
 * port of the VFLib3 c++ library package
 */
public class VFLib3 {
   
    static final int EMPTY = -1;

    static class MatchPair {
        int qn, tn; // query node, target node

        MatchPair () {
            this (EMPTY, EMPTY);
        }

        MatchPair (int qn, int tn) {
            this.qn = qn;
            this.tn = tn;
        }

        MatchPair (MatchPair mp) {
            this.qn = mp.qn;
            this.tn = mp.tn;
        }

        public String toString () {
            return "" + (qn + 1) + ":" + (tn + 1);
        }
    }

    public interface State {
        boolean isAtomCompatible (int qn /*query*/, int tn /*target*/);
       
        boolean next (MatchPair pair);
        boolean isFeasible (MatchPair pair);
        void add (MatchPair pair);
        boolean isGoal ();
        boolean isDead ();
        int getResultSize ();
        int[] getResult ();
        void backtrack ();
        Object clone ();
    }

    public interface MatchVisitor {
        boolean visit (int size, int[] result);
    }

    static class DefaultMatchVisitor implements MatchVisitor {
        List<int[]> hits = new ArrayList<int[]> ();
        Map<BitSet, int[]> unique = new HashMap<BitSet, int[]> ();

        DefaultMatchVisitor () {
        }

        public boolean visit (int size, int[] result) {
            BitSet bset = new BitSet ();
            for (int i = 0; i < result.length; ++i) {
                bset.set (result[i]);
            }
            if (!unique.containsKey (bset)) {
                unique.put (bset, result);
            }
            hits.add (result);
            return false;
        }

        public int[][] getAllHits () {
            return hits.toArray (new int[0][]);
        }

        public int[][] getUniqueHits () {
            return unique.values ().toArray (new int[0][]);
        }
    }

    /**
     * ***********************************************
     * Isomorphism
     * ************************************************
     */
    public static class VF2State implements State {
      
    	Chemical2 query, target;
        int origcorelen, corelen, qlen, tlen,
                qmaplen, tmaplen;
        int[] qcore, tcore;
        int[] qmap, tmap;
        int[] order = null;
        int qlastNode = EMPTY;

        BondTable tbtab;


        public VF2State (Chemical2 query, Chemical2 target) {
            this.query = query;
            this.target = target;
           
            initVars ();
        }

        public VF2State (VF2State s) { // copy constructor
            query = s.query;
            target = s.target;
           
            corelen = origcorelen = s.corelen;
            qlen = s.qlen;
            tlen = s.tlen;
            qmaplen = s.qmaplen;
            tmaplen = s.tmaplen;
            qcore = s.qcore;
            tcore = s.tcore;
            qmap = s.qmap;
            tmap = s.tmap;
            order = s.order;
            tbtab = s.tbtab;
            qlastNode = EMPTY;
        }

        void initVars () {
        	
        	
            qlen = query.getAtomCount ();
            tlen = target.getAtomCount ();

            qcore = new int[qlen];
            tcore = new int[tlen];
            qmap = new int[qlen];
            tmap = new int[tlen];

            for (int i = 0; i < qlen; ++i) {
                qcore[i] = EMPTY;
                qmap[i] = 0;
            }
            for (int i = 0; i < tlen; ++i) {
                tcore[i] = EMPTY;
                tmap[i] = 0;
            }

            origcorelen = 0;
            corelen = 0;
            qlastNode = EMPTY;

            tbtab = target.getBondTable();
        }

        public boolean isBondCompatible(Bond qb, Bond tb){
        	return qb.getBondType().equals(tb.getBondType());
        }

        public boolean isAtomCompatible (int qn, int tn) {
            Atom queryAtom = query.getAtom (qn);
			Atom targetAtom = target.getAtom (tn);

			return
					queryAtom.getAtomicNumber() ==targetAtom.getAtomicNumber()
                		&& queryAtom.hasAromaticBond() == targetAtom.hasAromaticBond();
        }

       

        /*
         * State interface
         */
        public boolean isGoal () {
            return corelen == qlen && corelen == tlen;
        }

        public boolean isDead () {
            return qlen != tlen || qmaplen != tmaplen;
        }

        public int getResultSize () {
            return corelen;
        }

        public int[] getResult () {
            return qcore;
        }

        public boolean next (MatchPair mp) {
            int qn = mp.qn;
            int tn = mp.tn;

            if (qn == EMPTY)
                qn = 0;

            if (tn == EMPTY)
                tn = 0;
            else
                tn++;

            if (qmaplen > corelen && tmaplen > corelen) {
                while (qn < qlen &&
                               (qcore[qn] != EMPTY || qmap[qn] == 0)) {
                    qn++;
                    tn = 0;
                }
            } else if (qn == 0 && order != null) {
                int i = 0;
                while (i < qlen && qcore[qn = order[i]] != EMPTY)
                    i++;
                if (i == qlen)
                    qn = qlen;
            } else {
                while (qn < qlen && qcore[qn] != EMPTY) {
                    qn++;
                    tn = 0;
                }
            }

            if (qmaplen > corelen && tmaplen > corelen) {
                while (tn < tlen &&
                               (tcore[tn] != EMPTY || tmap[tn] == 0)) {
                    tn++;
                }
            } else {
                while (tn < tlen && tcore[tn] != EMPTY) {
                    tn++;
                }
            }

            
            //System.out.println ("## next: " + (mp.qn + 1) + ":" + (mp.tn + 1) + " -> "
            //                            + (qn + 1) + ":" + (tn + 1));

            if (qn < qlen && tn < tlen) {
                mp.qn = qn;
                mp.tn = tn;

                return true;
            }

            return false;
        } // next ()


        public boolean isFeasible (MatchPair mp) {
        	//System.out.println ((mp.qn + 1) + ":" + (mp.tn + 1) + " ");
            if (!isAtomCompatible (mp.qn, mp.tn)) {
                //System.out.print ("incompat atoms ");
                return false;
            }

            BitSet qbs = new BitSet ();
            BitSet tbs = new BitSet ();

            Atom qa = query.getAtom (mp.qn);
            Atom ta = target.getAtom (mp.tn);

            int qc = 0, tc = 0;
            
            for (Bond qb : qa.getBonds()) {               
               
                int xa = query.indexOf (qb.getOtherAtom (qa));
               
                if (qcore[xa] != EMPTY) {
                    int a = qcore[xa];
                    if (!tbtab.bondExists(a,mp.tn) || !isBondCompatible(qb, tbtab.getBond(a,mp.tn))){
                   
                        //System.out.print ("target " + (a + 1) + " " + (mp.tn + 1) + " " + tbtab[a][mp.tn]);
                        return false;
                    }
                } else {
                    qbs.set (xa + 1);
                    if (qmap[xa] != 0){
                        //System.out.print (" qmap " + (xa + 1) + " = " + qmap[xa]);
                        ++qc;
                    }
                }
            }

            for (Bond tb : ta.getBonds()) {               
                int xa = target.indexOf (tb.getOtherAtom (ta));
                if (tcore[xa] == EMPTY) {                   
                    tbs.set (xa + 1);
                    if (tmap[xa] != 0){                      
                        //System.out.print (" tmap " + (xa + 1) + " = " + tmap[xa]);
                        ++tc;
                    }
                }
            }
            //System.out.print ("-- qc=" + qc + " tc=" + tc + " q=" + qbs + " t=" + tbs);

            // == isomorphism, <= sub graph isomorphism
            return qc <= tc;// && qnew <= tnew;
        }

        public void add (MatchPair mp) {
            ++corelen;
            qlastNode = mp.qn;

            if (qmap[mp.qn] == 0) {
                qmap[mp.qn] = corelen;
                ++qmaplen;
            }

            if (tmap[mp.tn] == 0) {
                tmap[mp.tn] = corelen;
                ++tmaplen;
            }

            qcore[mp.qn] = mp.tn;
            tcore[mp.tn] = mp.qn;

            Atom qa = query.getAtom (mp.qn);
            for (Bond b : qa.getBonds()) {
                int xa = query.indexOf (b.getOtherAtom (qa));
                if (qmap[xa] == 0) {
                    qmap[xa] = corelen;
                    ++qmaplen;
                }
            }

            Atom ta = target.getAtom (mp.tn);
            for (Bond b : ta.getBonds()) {
                int xa = target.indexOf (b.getOtherAtom (ta));
                if (tmap[xa] == 0) {
                    tmap[xa] = corelen;
                    ++tmaplen;
                }
            }
            
          
        }

        public void backtrack () {
            if (origcorelen >= corelen) {
                return;
            }

            // origcorelen < corelen
            if (qmap[qlastNode] == corelen) {
                qmap[qlastNode] = 0;
            }

            Atom qa = query.getAtom (qlastNode);
            for (Bond b : qa.getBonds()) {
                int xa = query.indexOf (b.getOtherAtom (qa));
                if (qmap[xa] == corelen) {
                    qmap[xa] = 0;
                }
            }

            int tlastNode = qcore[qlastNode];
            if (tmap[tlastNode] == corelen){
                tmap[tlastNode] = 0;
            }
            Atom ta = target.getAtom (tlastNode);
            for (Bond b : ta.getBonds()) {
                int xa = target.indexOf (b.getOtherAtom (ta));
                if (tmap[xa] == corelen) {
                    tmap[xa] = 0;
                }
            }

            qcore[qlastNode] = EMPTY;
            tcore[tlastNode] = EMPTY;
            corelen = origcorelen;
            qlastNode = EMPTY;
        }

        @Override
        public Object clone () {
            return new VF2State (this);
        }

        public String toString () {
            StringBuilder sb = new StringBuilder ();
            sb.append ("STATE: corelen:" + corelen + " mapping:");
            for (int i = 0; i < qcore.length && qcore[i] >= 0; ++i) {
                sb.append (" ");
                sb.append ((i + 1) + ":" + (qcore[i] + 1));
            }
            return sb.toString ();
        }
    } // VF2State


    /**
     * ********************************************
     * Subgraph isomorphism
     * *********************************************
     */
    public static class VF2SubState extends VF2State {

        public VF2SubState (Chemical2 query, Chemical2 target) {
            super (query, target);
        }

       

        public VF2SubState (VF2SubState s) {
            super (s);
        }

        @Override
        public boolean isGoal () {
            return corelen == qlen;
        }

        @Override
        public boolean isDead () {
            return qlen > tlen || qmaplen > tmaplen;
        }

        @Override
        public Object clone () {
            return new VF2SubState (this);
        }
    }

    /**
     * ********************************************
     * Automorphism
     * *********************************************
     */
    public static class VF2AutoState extends VF2State {
        private final GraphInvariant grinv;

        public VF2AutoState (Chemical2 mol) {
            super (mol, mol);
            grinv = mol.getGraphInvariant();
        }

        public VF2AutoState (VF2AutoState s) {
            super (s);
            grinv = s.grinv;
        }
        

        @Override
		public boolean isBondCompatible(Bond qb, Bond tb) {
			return true;
		}

		@Override
        public boolean isAtomCompatible (int qn, int tn) {
            return grinv.isAtomCompatible(qn, tn);
        }

     

        @Override
        public Object clone () {
            return new VF2AutoState (this);
        }
    }

    /**
     * main class VFLib3
     */

    protected State s;

    // can't instantiate this class directly; must use one of the static
    // methods
    protected VFLib3 (State s) {
        this.s = s;
    }

    static int[] clone (int[] result) {
        int[] r = new int[result.length];
        System.arraycopy (result, 0, r, 0, result.length);
        return r;
    }

    boolean find (State s, MatchVisitor visitor) {
        //System.out.print (">>> " + s);
        if (s.isGoal ()) {
            //System.out.println (" => GOAL!");
            return visitor.visit (s.getResultSize (), clone (s.getResult ()));
        }

        if (s.isDead ()) {
            //System.out.println (" => DEAD END!");
            return false;
        }

        //System.out.println (" => ...");

        MatchPair mp = new MatchPair ();
        boolean found = false;
        while (!found && s.next (mp)) {
            //System.out.print (" ** ");
            boolean cont = s.isFeasible (mp);
            //System.out.println (" " + s + " + " + mp + " => " + cont);

            if (cont) {
                State s1 = (State) s.clone ();
                s1.add (new MatchPair (mp));
                found = find (s1, visitor);
                s1.backtrack ();
                //System.out.println ("<<< " + s1);
            }
        }

        return found;
    }

    public boolean find (MatchVisitor visitor) {
        return find (s, visitor);
    }

    public boolean match () {
        return find (new MatchVisitor () {
            public boolean visit (int size, int[] result) {
                return true;
            }
        });
    }

    public int[][] findAll () {
        return findAll (false);
    }

    public int[][] findAll (boolean unique) {
        DefaultMatchVisitor visitor = new DefaultMatchVisitor ();
        find (visitor);
        return unique ? visitor.getUniqueHits () : visitor.getAllHits ();
    }



    public static VFLib3 subgraphIsomorphism
        (Chemical2 query, Chemical2 target) {
        return new VFLib3 (query == target ? new VF2AutoState (query)
                           : new VF2SubState (query, target));
    }

   

    public static VFLib3 automorphism (Chemical2 mol) {
        return new VFLib3 (new VF2AutoState (mol));
    }

   

    /*
     * *********************************************************
     * TEST CASES
     * *********************************************************
     *
    static void automorphism (String mol) throws Exception {
        MolHandler mh = new MolHandler(mol);
        mh.aromatize ();

        Molecule m = mh.getMolecule ();
        for (int i = 0; i < m.getAtomCount (); ++i) {
            MolAtom a = m.getAtom (i);
            a.setAtomMap (i + 1);
        }

        System.out.println(">>> "+m.toFormat("smiles:q"));
        System.out.println ("---- automorphism (VFLib3) ----");
        VFLib3 vf = automorphism (m);
        long start = System.currentTimeMillis ();
        int[][] hits = vf.findAll ();
        long end = System.currentTimeMillis ();
        for (int j = 0; j < hits.length; ++j) {
            int[] hit = hits[j];
            System.out.print ("Matched " + j + ":");
            for (int i = 0; i < hit.length; ++i) {
                System.out.print (" " + (i + 1) + ":" + (hit[i] + 1));
            }
            System.out.println ();
        }
        System.out.println ("## " + hits.length + " hit(s) found in "
                                    + String.format ("%1$3dms", end - start));


       
    }

    static void testIsomorphism () throws Exception {
        automorphism ("CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5");
    }

    static void testSubIsomorphism (String Q, String T) throws Exception {
        System.out.println ("---- subgraph isomorphism (VFLib3) ----");
        MolHandler mh = new MolHandler ();
        mh.setMolecule (Q);
        mh.aromatize ();

        Molecule query = mh.getMolecule ();
        mh.setMolecule (T);
        mh.aromatize ();
        Molecule target = mh.getMolecule ();

        for (int i = 0; i < query.getAtomCount (); ++i) {
            MolAtom a = query.getAtom (i);
            a.setAtomMap (i + 1);
        }
        for (int i = 0; i < target.getAtomCount (); ++i) {
            MolAtom a = target.getAtom (i);
            a.setAtomMap (i + 1);
        }

        VFLib3 vf = subgraphIsomorphism (query, target);
        long start = System.currentTimeMillis ();
        int[][] hits = vf.findAll (true);
        long end = System.currentTimeMillis ();
        for (int j = 0; j < hits.length; ++j) {
            int[] hit = hits[j];
            System.out.print ("** Matched " + (j + 1) + ":");
            for (int i = 0; i < hit.length; ++i) {
                System.out.print (" " + (i + 1) + ":" + (hit[i] + 1));
            }
            System.out.println ();
        }

        System.out.println ("## " + hits.length + " hit(s) found in "
                                    + String.format ("%1$3dms", end - start));

       

        System.out.println ("++  QUERY: " + query.toFormat ("smiles:q"));
        System.out.println ("++ TARGET: " + target.toFormat ("smiles:q"));
    }

    static void testSubIsomorphisms () throws Exception {
        String[][] tests = new String[][]{
            {"C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", "CC(=C)O[C@@]12CO[C@@H]1C[C@H](O)[C@]3(C)C2[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)C(NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](OC(=O)C)C3=O)C5(C)C)C"},
            {"C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", "CC(=O)O[C@H]1C[C@H]2OC[C@@]2(OC(=O)C)[C@H]3[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](O)C(=O)[C@]13C)C5(C)C)C"},
            {"C(C1=CC=CC=C1)C2=CC=CC=C2", "OC(=O)c1ccc(cc1)C2c3ccccc3c4ccccc24"},
            {"C(Cc1ccccc1)NCc2ccccc2","Oc1ccc(cc1)C2CN(Cc3cc(O)ccc23)C(=O)c4ccccc4"},
            {"O=C(NCc1ccccc1)c2ccccc2","Fc1ccc(cc1)C(=O)N2CCN3C(=O)c4ccccc4C23c5ccccc5"},
            {"C(Cc1ccccc1)N2CCNCC2","CN1CCN2CC(c3ccc(cc3)[N+]([O-])=O)c4ccccc4C2C1"},
            {"C(Oc1ccccc1)c2ccccc2","CCNc1ccc-2c(c1)C(Oc3cccc(OC)c-23)c4ccccc4"},
            {"C(Oc1ccccc1)c2ccccc2","CCNc1ccc-2c(c1)C(Oc3cccc(OC)c-23)c4ccccc4"},
            {"c1cn(cn1)C(c2ccccc2)c3ccccc3","c1cn(cn1)C2(c3ccccc3-c4ccccc24)c5ccccc5"},
            {"C(CNCCc1ccccc1)Cc2ccccc2","CC(Cc1ccccc1)NC2CC3c4c2cccc4CCc5ccccc35"},
            {"C(NSNCc1ccccc1)OCc2ccccc2","O=C1N(Cc2ccccc2)S(=O)(=O)N(COCc3ccccc3)c4ccccc14"},
            {"O=C(CCCCC1CCCCC1)NCc2ccccc2","CC12CC(O)C3C(CCC4=CC(=O)C=CC34C)C1CC[C@]2(O)C(O)C(=O)NCc5ccccc5"},
            {"C(Cc1ccccc1)Nc2ccccc2","O=C(Nc1ccccc1)C2C(=O)N3c4c2cccc4Sc5ccccc35"},
            {"C(c1ccccc1)[n+]2ccn3CCC(Cc23)c4ccoc4", "[Cl-].C(c1ccccc1)n2cc[n+]3C4CC(C(c23)c5ccccc45)(c6ccoc6)c7ccoc7"},
            {"C(CC1CCCC1)OC2COCC3C(CCCC23)C4CCC=CC4", "C[C@H]1[C@@H](O)[C@@]2(O)OC[C@@]34[C@@H](C[C@H]5C(=CC(=O)[C@@H](O)[C@]5(C)[C@@H]23)C)OC(=O)[C@H](OC(=O)CC6(O)C(C)(C)CCC6(C)C)[C@@H]14"},
            {"C(CC1CCCC1)OC2COCC3C(CCCC23)C4CCC=CC4", "C[C@H]1[C@@H](O)[C@@]2(O)OC[C@@]34[C@@H](C[C@H]5C(=CC(=O)[C@@H](O)[C@]5(C)[C@@H]23)C)OC(=O)[C@H](OC(=O)CC6(O)CCCC6)[C@@H]14"}
        };

        for (int i = 0; i < tests.length; ++i) {
            System.out.println (">>>>>>>>>>>>>>> TEST " + (i + 1));
            testSubIsomorphism (tests[i][0], tests[i][1]);
        }
    }
    
    public static void main (String[] argv) throws Exception {
        if (argv.length == 0) {
            testSubIsomorphisms ();
        }
        else if (argv.length == 1) {
            automorphism (argv[0]);
        }
        else {
            for (int i = 0; i < argv.length; ++i) 
                for (int j = i+1; j < argv.length; ++j)
                    testSubIsomorphism (argv[i], argv[j]);
        }
    }
    
    */
}
