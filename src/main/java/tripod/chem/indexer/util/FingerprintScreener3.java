package tripod.chem.indexer.util;

import java.io.PrintStream;

import java.util.Collections;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.BitSet;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * A generalized version of FingerprintScreener2
 * It's also slower!
 */
public class FingerprintScreener3 {
    static final Logger logger = 
	Logger.getLogger(FingerprintScreener3.class.getName());

    private static final double LN2 = Math.log(2.);
    static final int DEFAULT_MAX_DEPTH = 10;

    public interface ScreenVisitor {
	void matched (Object key, int[] fp);
    }

    // fixed-length bit vector 
    static class BitVec {
	int[] bits;
	int len; // length of bit vector (multiple of 32)

	BitVec (int len) {
	    bits = new int[len/32];
	    this.len = 32*bits.length;
	}
	BitVec (int[] bv) {
	    this.bits = bv;
	    this.len = 32*bv.length;
	}

	void clear () {
	    for (int i = 0; i < bits.length; ++i) {
		bits[i] = 0;
	    }
	}
	int length () { return len; }
	int cardinality () {
	    int card = 0;
	    for (int i = 0; i < bits.length; ++i) {
		card += Integer.bitCount(bits[i]);
	    }
	    return card;
	}
	boolean isEmpty () {
	    for (int i = 0; i < bits.length; ++i) {
		if (bits[i] != 0) return false;
	    }
	    return true;
	}
	boolean get (int pos) {
	    return (bits[pos/32] & (1<<(pos%32))) != 0;
	}
	void set (int pos) {
	    bits[pos/32] |= (1<<(pos%32));
	}
	boolean intersects (int[] bits) {
	    for (int i = 0; i < bits.length; ++i) {
		if ((this.bits[i] & bits[i]) != 0) 
		    return true;
	    }
	    return false;
	}
	// check this bit vector contains the given vector
	boolean contains (int[] bits) {
	    int N = Math.min(this.bits.length, bits.length);
	    for (int i = 0; i < N; ++i) {
		if ((this.bits[i] & bits[i]) != bits[i]) 
		    return false;
	    }
	    return true;
	}
	// check this bit vector is subset of the given vector
	boolean subset (int[] bits) {
	    int N = Math.min(this.bits.length, bits.length);
	    for (int i = 0; i < N; ++i) {
		if ((this.bits[i] & bits[i]) != this.bits[i]) 
		    return false;
	    }
	    return true;
	}

	boolean contains (BitVec bv) {
	    return contains (bv.bits);
	}
	boolean subset (BitVec bv) {
	    return subset (bv.bits);
	}

	public String toString () {
	    StringBuilder sb = new StringBuilder ();
	    if (!isEmpty ()) {
		for (int i = 0; i < len; ++i) {
		    if (get (i)) {
			if (sb.length() > 0) sb.append(",");
			sb.append(i);
		    }
		}
	    }
	    return "{"+sb.toString()+"}";
	}
	public int[] toArray () {
	    int[] ary = new int[cardinality ()];
	    for (int i = 0, j = 0; i < len; ++i) {
		if (get (i)) {
		    ary[j++] = i;
		}
	    }
	    return ary;
	}
    }

    static class SNode {
	BitVec bv;
	SNode left, right;
	SNode parent;
	List<FPV> values = new ArrayList<FPV>();

	SNode () {}
	SNode (FPV v) { values.add(v); }
	SNode (List<FPV> values) {
	    this.values.addAll(values);
	}

	boolean isLeaf () { return left == null && right == null; }
    }

    static class FPV {
	Object key;
	int[] fp;
	FPV (Object key, int[] fp) {
	    this.key = key;
	    this.fp = fp;
	}

	boolean get (int bit) { // in bit coordinate
	    return (fp[bit/32] & (1<<(bit%32))) != 0;
	}

	public String toString () {
	    StringBuilder sb = new StringBuilder ();
	    sb.append(key.toString()+"[");
	    sb.append(fp[0]+"="+Integer.toBinaryString(fp[0]));
	    for (int i = 1; i < fp.length; ++i) {
		sb.append(","+fp[i]+"="+Integer.toBinaryString(fp[i]));
	    }
	    sb.append("]");
	    return sb.toString();
	}
    }

    static class DeltaSorter implements Comparator<int[]> {
	public int compare (int[] a, int[] b) {
	    return b[0] - a[0];
	}
    }

    public static class SearchStats {
	int hitCount;
	int screenCount;
	List<int[]> signatures = new ArrayList<int[]>();

	SearchStats () {}
	SearchStats (int hitCount, int screenCount) {
	    this.hitCount = hitCount;
	    this.screenCount = screenCount;
	}

	public int getHitCount () { return hitCount; }
	public int getScreenCount () { return screenCount; }
	public int[][] getSignatures () {
	    return signatures.toArray(new int[0][]);
	}

	public static int signaturePos (int sig) {
	    return sig >> 1;
	}
	public static int signatureBit (int sig) {
	    return sig & 1;
	}
    }

    public static class ScreenStats {
	int nodeCount; // total node count = \sum_{k=0}^n 2^k
	int leafCount;
	int maxDepth;
	int minLeafSize, maxLeafSize;
	double avgLeafSize;
	int[] minsig, maxsig;

	ScreenStats () {}
	public int getLeafCount () { return leafCount; }
	public int getMaxDepth () { return maxDepth; }
	public int getMinLeafSize () { return minLeafSize; }
	public int getMaxLeafSize () { return maxLeafSize; }
	public double getAvgLeafSize () { return avgLeafSize; }

	public String toString () {
	    StringBuilder sb = new StringBuilder
		("##    NumNodes: "+nodeCount+"\n"
		 +"##    NumLeafs: "+leafCount+"\n"
		 +"##    MaxDepth: "+maxDepth+"\n"
		 +"## AvgLeafSize: "+avgLeafSize+"\n"
		 +"## MinLeafSize: "+minLeafSize+"\n"
		 +"## MaxLeafSize: "+maxLeafSize+"\n");
	    sb.append("##MinSignature:");
	    for (int i = 0; i < minsig.length; ++i) {
		sb.append(" ("+(minsig[i]>>1)+","+(minsig[i]&1)+")");
	    }
	    sb.append("\n##MaxSignature:");
	    for (int i = 0; i < maxsig.length; ++i) {
		sb.append(" ("+(maxsig[i]>>1)+","+(maxsig[i]&1)+")");
	    }
	    sb.append("\n");
	    return sb.toString();
	}
    }

    int dim;
    SNode root = null;
    int maxdepth;
    List<FPV> values = new ArrayList<FPV>();

    int[] _prof;
    DeltaSorter _sorter = new DeltaSorter ();

    public FingerprintScreener3 (int dim) {
	this (dim, DEFAULT_MAX_DEPTH);
    }

    public FingerprintScreener3 (int dim, int maxdepth) {
	if (dim < 32) {
	    throw new IllegalArgumentException
		("Bad dimension (< 32) specified: "+dim);
	}
	else if ((dim % 32) != 0) {
	    throw new IllegalArgumentException
		("Dimension "+dim+" is not a multiple of 32!");
	}
	this.dim = dim;
	this.maxdepth = maxdepth;

	_prof = new int[dim];
    }

    public void add (Object key, int[] fp) {
	if (fp.length*32 != dim) {
	    throw new IllegalArgumentException
		("Entry "+key+" has invalid dimension: "+(fp.length*32));
	}
	
	FPV val = new FPV (key, fp);
	if (root == null) {
	    root = new SNode (val);
	}
	else {
	    insert (val);
	}
	values.add(val);
    }

    public int size () { return values.size(); }

    static String toString (SNode node) {
	StringBuilder sb = new StringBuilder ();
	depthFirst (sb, 0, node, "**");
	return sb.toString();
    }

    static void depthFirst (StringBuilder sb, int depth, 
			    SNode n, String prefix) {
	if (n == null) {
	    return;
	}
	
	for (int i = 0; i <= depth; ++i) {
	    sb.append("  ");
	}

	if (prefix != null) {
	    sb.append(prefix);
	}
	if (n.isLeaf()) {
	    sb.append(" d="+depth);
	    for (FPV v : n.values) {
		sb.append(" "+v);
	    }
	    sb.append("\n");
	}
	else {
	    sb.append(" d="+depth+" bv="+n.bv+"\n");
	}
	depthFirst (sb, depth+1, n.left, "<");
	depthFirst (sb, depth+1, n.right,">");
    }

    void insert (FPV x) {
	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	int depth = 0;
	BitVec newvec = new BitVec (x.fp);

	while (!stack.isEmpty()) {
	    SNode v = stack.pop();
	    if (v.isLeaf()) {
		BitVec bv = new BitVec (dim);

		if (maxdepth <= 0 || depth < maxdepth) {
		    for (FPV w : v.values) {
			for (int i = 0; i < dim; ++i) {
			    if (w.get(i) != x.get(i)) {
				++_prof[i];
			    }
			}
		    }

		    List<int[]> deltas = new ArrayList<int[]>();
		    for (int i = 0; i < _prof.length; ++i) {
			if (_prof[i] > 0) {
			    deltas.add(new int[]{_prof[i], i});
			}
			_prof[i] = 0;
		    }
		    
		    if (!deltas.isEmpty()) {
			Collections.sort(deltas, _sorter);
			int[] d = deltas.get(0);
			bv.set(d[1]);
			if (deltas.size() > 1) {
			    d = deltas.get(deltas.size()-1);
			    bv.set(d[1]);
			}
		    }
		}

		if (bv.isEmpty()) {
		    v.values.add(x);
		}
		else {
		    // promote this leaf into internal node
		    v.bv = bv;
		    if (newvec.contains(bv)) {
			v.right = new SNode (x);
			v.left = new SNode (v.values);
		    }
		    else {
			v.right = new SNode (v.values);
			v.left = new SNode (x);
		    }
		    v.right.parent = v;
		    v.left.parent = v;
		    v.values.clear();
		}
	    }
	    else {
		SNode child = newvec.contains(v.bv) ? v.right : v.left;
		stack.push(child);
		++depth;
	    }
	}
    }

    public SearchStats search (int[] fp) {
	return search (fp, null);
    }

    int[] getSignature (SNode n) {
	List<Integer> sig = new ArrayList<Integer>();
	for (SNode p = n.parent; p != null; p = p.parent) {
	    int[] bits = p.bv.toArray();
	    for (int i = 0; i < bits.length; ++i) {
		sig.add((bits[i] << 1) | (p.left == n ? 0 : 1));
	    }
	    /*
	    if (p.left != n && p.right != n) {
		System.err.println("FATAL ERROR!");
		System.exit(1);
	    }
	    */
	    n = p;
	}
	int[] ps = new int[sig.size()];
	for (int i = 0; i < ps.length; ++i) { ps[i] = sig.get(i); }
	return ps;
    }

    public SearchStats search (int[] fp, ScreenVisitor visitor) {
	if (fp.length*32 != dim) {
	    throw new IllegalArgumentException
		("Query has invalid dimension: "+(fp.length*32));
	}
	if (root == null) {
	    throw new IllegalArgumentException
		("Signature tree hasn't been constructed yet!");
	}

	BitVec bv = new BitVec (fp);
	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	SearchStats stats = new SearchStats ();
	while (!stack.isEmpty()) {
	    SNode n = stack.pop();
	    if (n.isLeaf()) {
		int[] sig = getSignature (n);
		stats.signatures.add(sig);

		for (FPV v : n.values) {
		    int i = 0;
		    for (; i < fp.length; ++i) {
			if ((v.fp[i] & fp[i]) != fp[i]) {
			    break;
			}
		    }

		    boolean matched = i == fp.length;
		    if (matched) {
			if (visitor != null) {
			    visitor.matched(v.key, v.fp);
			}
			++stats.hitCount;
		    }
		}
		stats.screenCount += n.values.size();
	    }
	    else {
		stack.push(n.right);
		/*
		if ((fp[n.bit/32] & (1<<(n.bit%32))) == 0) {
		    stack.push(n.left);
		}
		*/
		if (!bv.contains(n.bv)) {
		    stack.push(n.left);
		}
	    }
	}

	return stats;
    }

    /*
     * perform linear scan
     */
    public int linear (int[] fp) { 
	return linear (fp, null);
    }

    public int linear (int[] fp, ScreenVisitor visitor) {
	int i, count = 0;
	for (FPV v : values) {
	    i = 0;
	    for (; i < fp.length; ++i) {
		if ((v.fp[i] & fp[i]) != fp[i]) {
		    break;
		}
	    }

	    boolean matched = i == fp.length;
	    if (matched) {
		if (visitor != null) {
		    visitor.matched(v.key, v.fp);
		}
		++count;
	    }
	}
	return count;
    }

    public void dump (PrintStream ps) {
	ps.println(toString (root));
    }

    public ScreenStats getScreenStats () {
	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	ScreenStats stats = new ScreenStats ();
	stats.minLeafSize = Integer.MAX_VALUE;

	while (!stack.isEmpty()) {
	    SNode n = stack.pop();
	    if (n.isLeaf()) {
		int d = stack.size();
		if (d > stats.maxDepth) {
		    stats.maxDepth = d;
		}

		int size = n.values.size();
		if (size < stats.minLeafSize) {
		    stats.minLeafSize = size;
		    stats.minsig = getSignature (n);
		}
		if (size > stats.maxLeafSize) {
		    stats.maxLeafSize = size;
		    stats.maxsig = getSignature (n);
		}
		stats.avgLeafSize += size;

		++stats.leafCount;
	    }
	    else {
		stack.push(n.right);
		stack.push(n.left);
	    }
	    ++stats.nodeCount;
	}
	stats.avgLeafSize /= stats.leafCount;

	return stats;
    }


    static void testRandomScreen (int dim, int size, int test) {
	java.util.Random rand = new java.util.Random (1l);
	FingerprintScreener3 screener = new FingerprintScreener3 (dim, 12);

	logger.info("generating "+size+" random "
		    +dim+"-bit fingerprints...");

	for (int i = 0; i < size; ++i) {
	    int[] fp = randomFp (rand, dim, rand.nextGaussian());
	    screener.add("key"+i, fp);
	}
	FingerprintScreener3.ScreenStats stats = screener.getScreenStats();
	logger.info("** Signature tree stats\n"+stats);

	logger.info("testing screening performance...");
	for (int i = 0; i < test; ++i) {
	    int[] fp = randomFp (rand, dim);
	    System.out.print("Query: ["+fp[0]);
	    for (int j = 1; j< fp.length; ++j) {
		System.out.print(","+fp[j]);
	    }
	    System.out.println("]");
	    double den = 0.;
	    for (int j = 0; j < fp.length; ++j) {
		den += Integer.bitCount(fp[j]);
	    }
	    den /= dim;
	    System.out.println
		("Query density: "+String.format("%1$.3f", den));

	    long start = System.currentTimeMillis();
	    SearchStats ss = screener.search(fp);
	    double time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("Signature screen took "
			       +String.format("%1$.3fs",time)
			       +" spanning "+ss.getSignatures().length
			       +" bucket(s); "+ss.getHitCount()+"/"
			       +ss.getScreenCount()+" found!");
	    start = System.currentTimeMillis();
	    int cnt2 = screener.linear(fp);
	    time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("Linear screen took "
			       +String.format("%1$.3fs",time)
			       +"; "+cnt2+" found!");
	    if (ss.getHitCount() != cnt2) {
		System.out.println("** FATAL: bug found for this query; "
				   +"please report this error!");
		System.exit(1);
	    }
	    System.out.println();
	}	
    }

    static int[] randomFp (java.util.Random rand, int dim) {
	return randomFp (rand, dim, rand.nextDouble());
    }

    static int[] randomFp (java.util.Random rand, int dim, double density) {
	int[] fp = new int[dim/32];
	if (density < 0) density *= -1.;
	int nb = (int)(density*dim + .5);
	for (int i = 0; i < nb; ++i) {
	    int b = rand.nextInt(dim); // uniformly turn on the bits
	    fp[b/32] |= 1<<(b%32);
	}
	return fp;
    }

    static void test () {
	FingerprintScreener3 fps = new FingerprintScreener3 (32);
	Random rand = new Random ();
	for (int i = 0; i < 10; ++i) {
	    int[] fp = new int[1];
	    fp[0] = rand.nextInt(20);
	    fps.add("key"+i, fp);
	}
	System.out.println("** Signature Tree **\n");
	fps.dump(System.out);
    }

    public static void main (String[] argv) throws Exception {
	testRandomScreen (512, 2000000, 100);
    }
}
