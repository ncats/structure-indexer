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

import java.util.logging.Logger;
import java.util.logging.Level;

public class FingerprintScreener2 {
    static final Logger logger = 
	Logger.getLogger(FingerprintScreener2.class.getName());

    private static final double LN2 = Math.log(2.);
    static final int DEFAULT_MAX_DEPTH = 10;

    public interface ScreenVisitor {
	void matched (Object key, int[] fp);
    }
    
    static class SNode {
	int bit;
	SNode left, right;
	SNode parent;
	List<FPV> values = new ArrayList<FPV>();

	SNode () {}
	SNode (FPV v) { values.add(v); }
	SNode (List<FPV> values) {
	    this.values = values;
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
    int[] _prof;
    int maxdepth;
    List<FPV> values = new ArrayList<FPV>();

    public FingerprintScreener2 (int dim) {
	this (dim, DEFAULT_MAX_DEPTH);
    }

    public FingerprintScreener2 (int dim, int maxdepth) {
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
	    sb.append(" d="+depth+" bit="+n.bit+"\n");
	}
	depthFirst (sb, depth+1, n.left, "<");
	depthFirst (sb, depth+1, n.right,">");
    }

    void insert (FPV x) {
	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	int depth = 0;
	while (!stack.isEmpty()) {
	    SNode v = stack.pop();
	    if (v.isLeaf()) {
		int k = dim;

		if (maxdepth <= 0 || depth < maxdepth) {
		    for (FPV w : v.values) {
			for (int i = 0; i < dim; ++i) {
			    if (w.get(i) != x.get(i)) {
				++_prof[i];
			    }
			}
		    }
		    int max = 0;
		    for (int i = 0; i < _prof.length; ++i) {
			if (_prof[i] > max) {
			    max = _prof[i];
			    k = i;
			}
			_prof[i] = 0;
		    }
		}

		if (k == dim) {
		    v.values.add(x);
		}
		else {
		    v.bit = k; // promote this leaf into internal node
		    if (x.get(k)) {
			v.right = new SNode (x);
			v.left = new SNode (v.values);
			/*
			v.left = createNode (v.values, 0);
			if (v.left != null) {
			    v.left.parent = v;
			}
			*/
		    }
		    else {
			v.right = new SNode (v.values);
			/*
			v.right = createNode (v.values, 0);
			if (v.right != null) {
			    v.right.parent = v;
			}
			*/
			v.left = new SNode (x);
		    }
		    v.right.parent = v;
		    v.left.parent = v;
		    v.values = null;
		}
	    }
	    else {
		SNode child = x.get(v.bit) ? v.right : v.left;
		stack.push(child);
		++depth;
	    }
	}
    }

    SNode createNode (List<FPV> subset, int depth) {
	if (subset == null || subset.isEmpty()) {
	    return null;
	}

	int[] freq = new int[dim];
	for (FPV v : subset) {
	    for (int i = 0; i < dim; ++i) {
		if (v.get(i)) {
		    ++freq[i];
		}
	    }
	}

	// find an index that split the set into as evenly in two
	// halves as possible
	int min = Integer.MAX_VALUE, split = -1, half = subset.size()/2;
	for (int i = 0; i < dim; ++i) {
	    int d = Math.abs(freq[i]-half);
	    if (d < min) {
		min = d;
		split = i;
	    }
	}

	SNode node = new SNode (subset);
	node.bit = split;
	if (min == half || split < 0 || depth >= maxdepth) {
	    // leaf node
	}
	else {
	    List<FPV> lsub = new ArrayList<FPV>();
	    List<FPV> rsub = new ArrayList<FPV>();
	    for (FPV v : subset) {
		if (v.get(split)) {
		    rsub.add(v);
		}
		else {
		    lsub.add(v);
		}
	    }
	    
	    node.left = createNode (lsub, depth+1);
	    if (node.left != null) {
		node.left.parent = node;
	    }

	    node.right = createNode (rsub, depth+1);
	    if (node.right != null) {
		node.right.parent = node;
	    }
	}

	return node;
    }

    public SearchStats search (int[] fp) {
	return search (fp, null);
    }

    int[] getSignature (SNode n) {
	List<Integer> sig = new ArrayList<Integer>();
	for (SNode p = n.parent; p != null; p = p.parent) {
	    sig.add((p.bit << 1) | (p.left == n ? 0 : 1));
	    if (p.left != n && p.right != n) {
		System.err.println("FATAL ERROR!");
		System.exit(1);
	    }
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

	LinkedList<SNode> stack = new LinkedList<SNode>();
	stack.push(root);

	SearchStats stats = new SearchStats ();
	while (!stack.isEmpty()) {
	    SNode n = stack.pop();
	    if (n.isLeaf()) {
		int[] sig = getSignature (n);
		stats.signatures.add(sig);

		screen (stats, n, fp, visitor);
	    }
	    else {
		stack.push(n.right);
		if ((fp[n.bit/32] & (1<<(n.bit%32))) == 0) {
		    stack.push(n.left);
		}
	    }
	}

	return stats;
    }

    void screen (SearchStats stats, SNode n, 
		 int[] fp, ScreenVisitor visitor) {
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
	FingerprintScreener2 screener = new FingerprintScreener2 (dim, 12);

	logger.info("generating "+size+" random "
		    +dim+"-bit fingerprints...");

	for (int i = 0; i < size; ++i) {
	    int[] fp = randomFp (rand, dim, rand.nextGaussian());
	    screener.add("key"+i, fp);
	}
	FingerprintScreener2.ScreenStats stats = screener.getScreenStats();
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
	FingerprintScreener2 fps = new FingerprintScreener2 (32);
	Random rand = new Random ();
	for (int i = 0; i < 20; ++i) {
	    int[] fp = new int[1];
	    fp[0] = rand.nextInt(20);
	    fps.add("key"+i, fp);
	}
	System.out.println("** Signature Tree **\n");
	fps.dump(System.out);
	ScreenStats stats = fps.getScreenStats();
	logger.info("** Signature tree stats\n"+stats);	
    }

    public static void main (String[] argv) throws Exception {
	testRandomScreen (1024, 1000000, 100);
	//test ();
    }
}
