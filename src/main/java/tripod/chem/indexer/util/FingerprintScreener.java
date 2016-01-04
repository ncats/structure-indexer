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

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * A fast fingerprint screener based on k-d tree.  For a high-level
 * description of the code, please see the blog post at
 *
 *      http://tripod.nih.gov/?p=329
 *
 * There are lots of room for improvements here, so we'd love to 
 * know of any improvements made to the implementation. 
 */
public class FingerprintScreener {

    static final Logger logger = 
	Logger.getLogger(FingerprintScreener.class.getName());

    static final int DEBUG;
    static final int MIN_LEAF_SIZE;

    static {
	int debug = 0, leaf = 100;
	try {
	    debug = Integer.getInteger("screen.debug", 0);
	    leaf = Integer.getInteger("screen.leaf", 100);
	}
	catch (Exception ex) {
	}
	DEBUG = debug;
	MIN_LEAF_SIZE = leaf;
    }


    public interface ScreenVisitor {
	void matched (Object key, int[] fp);
    }

    static class FpSorter implements Comparator<FPV> {
	int axis;
	FpSorter (int axis) {
	    this.axis = axis;
	}

	// Comparator interface
	int compare1 (FPV v1, FPV v2) {
	    int d = Integer.bitCount(v1.fp[axis]) 
		- Integer.bitCount(v2.fp[axis]);
	    if (d == 0) {
		long x1 = v1.getFp(axis);
		long x2 = v2.getFp(axis);
		if (x1 < x2) d = -1;
		else if (x1 > x2) d = 1;
	    }
	    return d;
	}
	int compare0 (FPV v1, FPV v2) {
	    long x1 = v1.getFp(axis);
	    long x2 = v2.getFp(axis);
	    if (x1 < x2) return -1;
	    else if (x1 > x2) return 1;
	    return 0;
	}
	public int compare (FPV v1, FPV v2) {
	    return compare1 (v1, v2);
	}
    }

    /*
     * A basic k-d tree structure
     */
    static class KdNode {
	/*
	 * fp value at the designated split
	 */
	int median; 

	List<FPV> values; // values for this node
	KdNode left, right; // children
	KdNode next; // next leaf node

	KdNode (int median, List<FPV> values) {
	    this.median = median;
	    this.values = values;
	}

	KdNode getChild1 (int v) {
	    int d = Integer.bitCount(median) - Integer.bitCount(v);
	    if (d < 0) return right;
	    else if (d > 0) return left;
	    else {
		long x = median & 0xffffffffl;
		long y = v & 0xffffffffl;
		if (x < y) return right;
		return left;
	    }
	}

	KdNode getChild0 (int v) {
	    long x = median & 0xffffffffl;
	    long y = v & 0xffffffffl;
	    if (x < y) return right;
	    return left;
	}

	KdNode getChild (int v) {
	    return getChild1 (v);
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

	// return fingerprint as unsigned value
	long getFp (int pos) {
	    return fp[pos] & 0xffffffffl;
	}

	public String toString () {
	    StringBuilder sb = new StringBuilder ();
	    sb.append(key.toString()+"[");
	    sb.append(getFp (0));
	    for (int i = 1; i < fp.length; ++i) {
		sb.append(","+getFp (i));
	    }
	    sb.append("]");
	    return sb.toString();
	}
    }

    final int[] dimen;
    List<FPV> values = new ArrayList<FPV>();
    KdNode root;

    public FingerprintScreener (int dimen) {
	if (dimen == 0) {
	    throw new IllegalArgumentException ("Bad dimension: "+dimen);
	}
	this.dimen = new int[dimen];
	for (int i = 0; i < dimen; ++i) {
	    this.dimen[i] = i;
	}
    }

    public FingerprintScreener (int[] dimen) {
	if (dimen == null || dimen.length == 0) {
	    throw new IllegalArgumentException ("Bad dimension array");
	}
	this.dimen = dimen;
    }

    public void add (Object key, int[] fp) {
	if (fp.length != dimen.length) {
	    throw new IllegalArgumentException
		("Invalid fingerprint length: "+fp.length);
	}

	if (root != null) {
	    throw new IllegalArgumentException
		("Can't no longer add since a tree has already "
		 +"been constructed!");
	}

	values.add(new FPV (key, fp));
    }

    public int getDim () { return dimen.length; }
    public int size () { return values.size(); }

    // construct the kd tree
    int _maxdepth, _leafsize, _leafcount;
    KdNode _lastLeaf;

    public void build () {
	build (MIN_LEAF_SIZE);
    }

    public void build (int leafSize) {
	if (root == null) {
	    _maxdepth = 0;
	    _leafsize = Math.max(1, leafSize);
	    _leafcount = 0;
	    _lastLeaf = null;
	    
	    root = build (values, 0);
	}
    }

    KdNode build (List<FPV> subset, final int depth) {
	if (subset == null || subset.isEmpty()) {
	    return null;
	}

	if (depth > _maxdepth) {
	    _maxdepth = depth;
	}

	final int axis = dimen[depth % dimen.length];
	int size = subset.size();

	// sort the array based on the unsigned fp values 
	Collections.sort(subset, new FpSorter (axis));

    	int median = subset.get(size/2).fp[axis];
	KdNode node = new KdNode (median, subset);

	if (size <= _leafsize) { // leaf node
	    ++_leafcount;
	    if (_lastLeaf != null) {
		_lastLeaf.next = node;
	    }
	    _lastLeaf = node;
	}
	else {
	    // < median
	    int length = size/2;
	    List<FPV> left = new ArrayList<FPV>(length);
	    Iterator<FPV> it = subset.iterator();
	    for (int i = 0; i < length; ++i) {
		left.add(it.next());
	    }
	    node.left = build (left, depth+1);
	    
	    // >= median
	    length = size - size/2;
	    List<FPV> right = new ArrayList<FPV>(length);
	    for (int i = 0; i < length; ++i) {
		right.add(it.next());
	    }
	    node.right = build (right, depth+1);
	}

	return node;
    }

    void debug (PrintStream ps) {
	ps.println("** maxdepth: "+_maxdepth);
	ps.println(toString (root));
	ps.println("## DepthFirst traversal order");
	for (Enumeration<FPV> en = depthFirstEnum (root); 
	     en.hasMoreElements(); ) {
	    FPV v = en.nextElement();
	    ps.println(v);
	}
    }

    public void stats () {
	logger.info("## maxdepth="+_maxdepth+" leafs="
		    +_leafcount+" leafsize="+_leafsize+" dim="+getDim());
	System.err.print("## split order: ["+dimen[0]);
	for (int i = 1; i < dimen.length; ++i) {
	    System.err.print(","+dimen[i]);
	}
	System.err.println("]");
    }

    String toString (KdNode node) {
	StringBuilder sb = new StringBuilder ();
	depthFirst (sb, 0, node);
	return sb.toString();
    }

    void depthFirst (StringBuilder sb, int depth, KdNode n) {
	if (n == null) {
	    return;
	}
	
	for (int i = 0; i <= depth; ++i) {
	    sb.append("  ");
	}
	if (n.isLeaf()) {
	    sb.append("* d="+depth+" a="+(depth%dimen.length));
	    for (FPV v : n.values) {
		sb.append(" "+v);
	    }
	    sb.append("\n");
	}
	else {
	    Iterator<FPV> it = n.values.iterator();
	    sb.append("+ m="+n.median+" d="+depth+" a="+(depth%dimen.length)
		      +" "+n.values.size() +"={"+values.indexOf(it.next()));
	    while (it.hasNext()) {
		sb.append(","+values.indexOf(it.next()));
	    }
	    sb.append("}\n");
	}
	depthFirst (sb, depth+1, n.left);
	depthFirst (sb, depth+1, n.right);
    }
    
    Enumeration<FPV> depthFirstEnum (KdNode node) {
	List<FPV> en = new ArrayList<FPV>();
	depthFirstEnum (node, en);
	return Collections.enumeration(en);
    }
    
    void depthFirstEnum (KdNode n, List<FPV> en) {
	if (n.isLeaf()) {
	    for (FPV v : n.values) {
		en.add(v);
	    }
	}
	else {
	    if (n.left != null) {
		depthFirstEnum (n.left, en);
	    }
	    if (n.right != null) {
		depthFirstEnum (n.right, en);
	    }
	}
    }
    

    public int screen (int[] fp) {
	return screen (fp, null);
    }

    public int screen (int[] fp, ScreenVisitor visitor) {
	if (fp.length != dimen.length) {
	    throw new IllegalArgumentException
		("Invalid fingerprint length: "+fp.length);
	}

	if (root == null) {
	    throw new IllegalArgumentException
		("Please construct the KdTree first!");
	}

	return _screen (fp, visitor);
    }

    KdNode getLeaf (int[] fp) {
	KdNode node = root, n = node;

	int depth = 0;
	do {
	    node = n;
	    n = n.getChild(fp[dimen[depth++ % dimen.length]]);
	}
	while (n != null);
	return node;
    }

    int _screen (int[] fp, ScreenVisitor visitor) {
	KdNode node = getLeaf (fp);

	if (DEBUG == 1) {
	    System.err.print("searching for ["+fp[0]);
	    for (int i = 1; i < fp.length; ++i) {
		System.err.print(","+fp[i]);
	    }
	    System.err.println("]");
	    System.err.println("## Candidates...");
	}

	int count = 0, i, size = 0, screened = 0;
	KdNode next;
	do {
	    next = node.next;

	    for (FPV val : node.values) {
		i = 0;
		for (; i < fp.length; ++i) {
		    /*
		     * change the following expression to match the desire
		     * behavior; e.g., val.fp[i] < fp[i]
		     */
		    if ((val.fp[i] & fp[i]) != fp[i]) {
			break;
		    }
		}
		
		boolean matched = i == fp.length;
		if (matched) {
		    if (visitor != null) {
			visitor.matched(val.key, val.fp);
		    }
		    ++count;
		}
	    
		if (DEBUG == 1) {
		    System.err.println(val+(matched ?"**":""));
		}
		++screened;
	    }
	    ++size;

	    node = next;
	}
	while (node != null);

	if (true/*DEBUG == 1*/) {
	    logger.info(size+"/"+_leafcount+" nodes scanned => "
			+screened+"/"+values.size()+" ("
			+String.format("%1$.3f", (double)screened/values.size())+") screened");
	}

	return count;
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


    static FingerprintScreener createRandomScreenDb (int dim, int size) {
	java.util.Random rand = new java.util.Random ();
	FingerprintScreener screener = new FingerprintScreener (dim);
	logger.info("generating "+size+" random "
		    +(32*dim)+"-bit fingerprints...");

	for (int i = 0; i < size; ++i) {
	    screener.add("key"+i, randomFp (rand, dim, rand.nextGaussian()));
	}

	long start = System.currentTimeMillis();
	screener.build(500);
	logger.info("Time to construct Kd tree of size "+screener.size()+" is "
		    +String.format("%1$.3fs", 
				   1e-3*(System.currentTimeMillis()-start)));
	return screener;
    }

    static int[] randomFp (java.util.Random rand, int dim) {
	return randomFp (rand, dim, rand.nextDouble());
    }

    static int[] randomFp (java.util.Random rand, int dim, double density) {
	int[] fp = new int[dim];
	int size = 32*dim;
	if (density < 0) density *= -1.;
	int nb = (int)(density*size + .5);
	for (int i = 0; i < nb; ++i) {
	    int b = rand.nextInt(size); // uniformly turn on the bits
	    fp[b/32] |= 1<<(b%32);
	}
	return fp;
    }

    public static void main (String[] argv) throws Exception {
	java.util.Random rand = new java.util.Random (1);

	int dim = 32;// simulate 1024-bit fingerprint
	int size = 1000000; // a reasonable size that fits in memory

	FingerprintScreener screener = createRandomScreenDb (dim, size);
	//screener.debug(System.err);
	screener.stats();

	for (int i = 0; i < 100; ++i) {
	    int[] fp = randomFp (rand, dim);
	    System.out.print("Query: ["+fp[0]);
	    for (int j = 1; j< fp.length; ++j) {
		System.out.print(","+fp[j]);
	    }
	    System.out.println("]");
	    long start = System.currentTimeMillis();
	    int cnt = screener.screen(fp);
	    double time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("   K-d screen took "
			       +String.format("%1$.3fs",time)
			       +"; "+cnt+" found!");
	    start = System.currentTimeMillis();
	    int cnt2 = screener.linear(fp);
	    time = 1e-3*(System.currentTimeMillis()-start);
	    System.out.println("Linear screen took "
			       +String.format("%1$.3fs",time)
			       +"; "+cnt2+" found!");
	    if (cnt != cnt2) {
		System.out.println("** FATAL: bug found for this query; "
				   +"please report this error!");
		System.exit(1);
	    }
	    System.out.println();
	}
    }
}
