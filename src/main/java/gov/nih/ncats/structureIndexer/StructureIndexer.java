package gov.nih.ncats.structureIndexer;

import static org.apache.lucene.document.Field.Store.NO;
import static org.apache.lucene.document.Field.Store.YES;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.core.KeywordAnalyzer;
import org.apache.lucene.analysis.miscellaneous.PerFieldAnalyzerWrapper;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.DoubleField;
import org.apache.lucene.document.IntField;
import org.apache.lucene.document.LongField;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.facet.FacetField;
import org.apache.lucene.facet.FacetResult;
import org.apache.lucene.facet.Facets;
import org.apache.lucene.facet.FacetsCollector;
import org.apache.lucene.facet.FacetsConfig;
import org.apache.lucene.facet.LabelAndValue;
import org.apache.lucene.facet.taxonomy.FastTaxonomyFacetCounts;
import org.apache.lucene.facet.taxonomy.TaxonomyReader;
import org.apache.lucene.facet.taxonomy.directory.DirectoryTaxonomyReader;
import org.apache.lucene.facet.taxonomy.directory.DirectoryTaxonomyWriter;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexReaderContext;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.Terms;
import org.apache.lucene.index.TermsEnum;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.DisjunctionMaxQuery;
import org.apache.lucene.search.Filter;
import org.apache.lucene.search.FilteredQuery;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.MatchAllDocsQuery;
import org.apache.lucene.search.NumericRangeQuery;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.NIOFSDirectory;
import org.apache.lucene.store.NoLockFactory;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.Version;

import gov.nih.ncats.common.io.IOUtil;
import gov.nih.ncats.common.util.CachedSupplier;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.MolwitchException;
import gov.nih.ncats.molwitch.Bond.BondType;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinters.FingerprintSpecification;
import gov.nih.ncats.molwitch.io.ChemFormat;



public class StructureIndexer {
    static final Logger logger =
        Logger.getLogger(StructureIndexer.class.getName());

   // static final Molecule POISON_PILL = new Molecule ();
    static final Version LUCENE_VERSION = Version.LATEST;
    /**
     * Fields for each document
     */
    public static final String FIELD_TEXT = "text";
    public static final String FIELD_NAME = "name";
    public static final String FIELD_ID = "_id";
    public static final String FIELD_SOURCE = "_source";
    public static final String FIELD_CODEBOOK = "_codebook";
    public static final String FIELD_FINGERPRINT_SUB = "_fingerprint_sub";
    public static final String FIELD_FINGERPRINT_SIM = "_fingerprint_sim";
    public static final String FIELD_FIELDS = "_fields";
    public static final String FIELD_FORMULA = "_formula";
    // fingerprint pop count    
    public static final String FIELD_POPCNT = "_popcount";
    public static final String FIELD_MOLFILE = "_molfile";
    public static final String FIELD_MOLWT = "_molwt";
    public static final String FIELD_NATOMS = "_natoms";
    public static final String FIELD_NBONDS = "_nbonds";

    static final String FIELD_DICT = "_dict";
    static final String FIELD_CODE = "_code";
    
    static final int FPSIZE = 16;
    static final int FPBITS = 2;
    static final int FPDEPTH = 6;
    
    static final int CODESIZE = 8; // 8-bit or 256
    static final int CODEBOOKS = 256;

    static final char[] ALPHA = {
        'Q','X','Y','Z','U','V','W'
    };
    static final String[] CODEWORDS;
    static {
        List<String> codes = new ArrayList<String>();
        // 2^8 = 256 < 7^3 (343)
        for (int i = 0; i < ALPHA.length; ++i)
            for (int j = 0; j < ALPHA.length; ++j)
                for (int k = 0; k < ALPHA.length; ++k)
                    codes.add(new String
                              (new char[]{ALPHA[i],ALPHA[j],ALPHA[k]}));
        CODEWORDS = codes.toArray(new String[codes.size()]);
    }
    
    static public class Codebook {
        final String name;
        int[] dict;
        int[][] eqv;
        int[] counts;
        
        public Codebook (int fpsize) {
            byte[] buf = new byte[4];
            Random rand = new Random ();            
            rand.nextBytes(buf);
            this.name = "CB"+toHex(buf);
            
            BitSet used = new BitSet ();
            int[] d = new int[CODESIZE];
            for (int i = 0; i < d.length; ++i) {
                int b = rand.nextInt(fpsize);
                if (!used.get(b)) {
                    used.set(b);
                    d[i] = b;
                }
            }
            setDictionary (d);
        }
        
        public Codebook (String name, int[] dict) {
            this.name = name;
            setDictionary (dict);
        }

        public int size () { return counts.length; }
        
        static Document instrument (Codebook cb) {
            Document doc = new Document ();
            doc.add(new StringField (FIELD_ID, cb.getName(), YES));
            for (int i = 0; i < cb.dict.length; ++i) {
                doc.add(new IntField (FIELD_DICT, cb.dict[i], YES));
            }
            
            for (int i = 1; i < cb.counts.length; ++i) {
                doc.add(new IntField (FIELD_CODE+"_"+i, cb.counts[i], YES));
            }
            //logger.info("++ "+load (doc));
            return doc;
        }

        static Codebook load (Document doc) {
            String name = doc.get(FIELD_ID);
            IndexableField[] fields = doc.getFields(FIELD_DICT);
            int[] dict = new int[fields.length];
            for (int i = 0; i < dict.length; ++i)
                dict[i] = fields[i].numericValue().intValue();

            Codebook cb = new Codebook (name, dict);
            for (int i = 1; i < cb.counts.length; ++i) {
                IndexableField f = doc.getField(FIELD_CODE+"_"+i);
                if (f != null) {
                    cb.counts[i] = f.numericValue().intValue();
                }
                else {
                    logger.warning
                        (name+": field "+FIELD_CODE+"_"+i+" is null!");
                }
            }
            return cb;
        }

        public String getName () { return name; }
        public int[] getDictionary () { return dict; }
        public void setDictionary (int[] dict) {
            if (dict == null || dict.length == 0)
                throw new IllegalArgumentException ("Invalid dictionary!");

            counts = new int[1<<dict.length];
            BitSet[] bsets = new BitSet[counts.length];
            for (int i = 0; i < bsets.length; ++i)
                bsets[i] = new BitSet ();

            for (int i = 1; i < bsets.length; ++i) {
                for (int j = 1; j <  bsets.length; ++j) 
                    if ((i & j) == j) 
                        bsets[j].set(i);
            }
            
            eqv = new int[bsets.length][];
            for (int i = 0; i < eqv.length; ++i) {
                BitSet bs = bsets[i];
                eqv[i] = new int[bs.cardinality()];
                for (int k = 0, j = bs.nextSetBit(0);
                     j>=0; j = bs.nextSetBit(j+1))
                    eqv[i][k++] = j;
            }
            
            this.dict = dict;
        }

       
        public int[] apply(Fingerprint fingerPrint){
        	 int code = encode (fingerPrint);
             return code == 0 ? null : eqv[code];
        }
        public int[] apply (byte[] fp) {
            int code = encode (fp);
            return code == 0 ? null : eqv[code];
        }

        public int count (int code) {
            return counts[code];
        }
        
        public synchronized int incr (int code) {
            return ++counts[code];
        }

        public synchronized int decr (int code) {
            return --counts[code];
        }

        public String encode (int code) {
            return name+CODEWORDS[code];
        }

        public int decode (String code) {
            if (code.startsWith(name)) {
                String id = code.substring(0, name.length());
                for (int i = 0; i < CODEWORDS.length; ++i)
                    if (id.equals(CODEWORDS[i]))
                        return i;
            }
            return -1;
        }

        protected synchronized void adjustCounts (IndexSearcher searcher)
            throws IOException {
            int total = searcher.getIndexReader().numDocs();
            for (int i = 1; i < counts.length; ++i) {
                TermQuery tq = new TermQuery
                    (new Term (FIELD_CODEBOOK, encode (i)));
                TopDocs hits = searcher.search(tq, total);
                counts[i] = hits.totalHits;
            }
        }

        public int encode (Fingerprint fp) {
            int code = 0;
            BitSet bs = fp.toBitSet();
            for (int i = 0; i < dict.length; ++i) {
//                if (get(fpArray, dict[i])){
                if(bs.get(dict[i])){
                    code |= 1 << i;
                }
            }

             return code & 0xff;
        }
        
        public int encode (byte[] fp) {
            int code = 0;
            for (int i = 0; i < dict.length; ++i) {
                if (get (fp, dict[i]))
                    code |= 1<<i;
            }
            return code & 0xff;
        }

        public int encode (BitSet fp) {
            int code = 0;
            for (int i = 0; i < dict.length; ++i) {
                if (fp.get(dict[i]))
                    code |= i << i;
            }
            return code & 0xff;
        }

        public String toString () {
            StringBuilder sb = new StringBuilder ("{name="+name+",bits=");
            for (int i = 0; i < dict.length; ++i) {
                sb.append(dict[i]);
                if (i+1 < dict.length) sb.append(" ");
            }
            sb.append(",counts=");
            for (int i = 0; i < counts.length; ++i) {
                if (counts[i] > 0) {
                    sb.append(i+":"+counts[i]);
                    if (i+1 < counts.length) sb.append(" ");
                }
            }
            sb.append("}");
            return sb.toString();
        }
    }

    static boolean get (int[] fp, int bit) {
        return (fp[bit/32] & ((1 << (31-(bit % 32))))) != 0;
    }
    
    static boolean get (byte[] fp, int bit) {
        return (fp[bit/8] & ((1 << (7-(bit % 8))))) != 0;
    }

    static String toHex (byte[] buf) {
        StringBuilder sb = new StringBuilder ();
        for (int i = 0; i < buf.length; ++i)
            sb.append(String.format("%1$02X", buf[i] & 0xff));
        return sb.toString();
    }
    
    static class Payload {
        final int id;   
        final Document doc;
        Chemical mol;
        Fingerprint fpSub;
        Fingerprint fpSim;

        Payload (int id, Document doc) {
            this.id = id;
            this.doc = doc;
        }
        
        Payload () {
            id = -1;
            doc = null;
        }

        public Fingerprint getFpSub () {
        	
            if (fpSub == null) {
                BytesRef ref = doc.getBinaryValue(FIELD_FINGERPRINT_SUB);
                byte[] bytes;
                if (ref.offset > 0) {
                	bytes = new byte[ref.length];
                    System.arraycopy(ref.bytes, ref.offset, fpSub, 0, ref.length);
                }
                else {
                	bytes = ref.bytes;
                }
                
                fpSub = new Fingerprint(bytes);
            }
            return fpSub;
        }
        public Fingerprint getFpSim () {
        	
            if (fpSim == null) {
                BytesRef ref = doc.getBinaryValue(FIELD_FINGERPRINT_SIM);
                byte[] bytes;
                if (ref.offset > 0) {
                	bytes = new byte[ref.length];
                    System.arraycopy(ref.bytes, ref.offset, fpSim, 0, ref.length);
                }
                else {
                	bytes = ref.bytes;
                }
                
                fpSim = new Fingerprint(bytes);
            }
            return fpSim;
        }
        public int getId () { return id; }
        public Document getDoc () { return doc; }
        public Chemical getMol () {
            if (mol == null) {
                String mol = doc.get(FIELD_MOLFILE);


                try {

                    this.mol = Chemical.parseMol(mol);
//                    try{
//                        this.mol.aromatize();
//                    }catch(Exception e){
//                        //ignore?
//                    }
                    this.mol.setName(doc.get(FIELD_ID));
                    for (IndexableField f : doc.getFields(FIELD_FIELDS)) {
                        String v = doc.get(f.stringValue());
                        this.mol.setProperty(f.stringValue(), v);
                    }
                }
                catch (Exception ex) {
                	ex.printStackTrace();
                	System.err.println("bbadmol=\n"+mol);
                    throw new RuntimeException
                        ("Document "+doc.get(FIELD_ID)+" contains bogus "
                         +"field "+FIELD_MOLFILE+"!\n" , ex);
                }
            }
            return mol;
        }
    }
    static final Payload POISON_PAYLOAD = new Payload ();

    static public class Result implements Comparable<Result> {
        final int id;
        final Document doc;
        final CachedSupplier<Chemical> mol;
        final int[] hit;
        final Double similarity;
        
        Result (Payload payload, Double similarity, int[] hit) {
            id = payload.getId();
            doc = payload.getDoc();


            //the input hit is 0-based offset must make it 1-based position
            if(hit==null){
                this.hit= null;
                mol = CachedSupplier.of(()->payload.getMol());
            }else {
                this.hit = new int[hit.length];
                for (int i = 0; i < hit.length; i++) {
                    int old = hit[i];
                    if (old >= 0) {
                        //atoms[hits[i]].setAtomMap(i+1);
//                        mol.getAtom(hit[i]).setAtomToAtomMap( i+1);
                        this.hit[i] = old + 1;
                    }
                }
                mol=CachedSupplier.of(()->{
                    Chemical mol = payload.getMol();
                    for (int i = 0; i < this.hit.length; i++) {
                        if(this.hit[i] >=0){
                            mol.getAtom(this.hit[i]-1).setAtomToAtomMap( i+1);
                        }
                    }
                    return mol;
                });

            }
            this.similarity = similarity;
        }

        public int[] getHits(){
            if(hit==null) {
                return null;
            }
            return Arrays.copyOf(hit, hit.length);
        }
        
        Result (Payload payload) {
            this (payload, null,null);
        }
        
        Result () {
            id = Integer.MAX_VALUE;
            doc = null;
            mol = null;
            hit= null;
            similarity = null;
        }

        public int compareTo (Result r) {
            int d = 0;
            if (similarity != null && r.similarity != null) {
                double dif = r.similarity - similarity;
                if (dif > 0.) d = 1;
                else if (dif < 0.) d = -1;
            }
            
            if (d == 0 && doc != null && r.doc != null) {
                // tiebreaker
                IndexableField f1 = doc.getField(FIELD_NATOMS);
                IndexableField f2 = r.doc.getField(FIELD_NATOMS);
                if (f1 != null && f2 != null) {
                    d = f1.numericValue().intValue()
                        - f2.numericValue().intValue();
                }
            }

            if (d == 0)
                d = id - r.id;
            
            //logger.info(id+" ("+similarity+") vs "+r.id+" ("+r.similarity+") => "+d);       
            return d;
        }

        public String get (String field) { return doc.get(field);}
        public String getId () { return doc.get(FIELD_ID); } 
        public String getSource () { return doc.get(FIELD_SOURCE); }
        public Double getSimilarity () { return similarity; }
        public Chemical getMol () { return mol.get(); }
        public Document getDoc () { return doc; }
        public String[] getFields () {
            List<String> fields = new ArrayList<String>();
            for (IndexableField f : doc.getFields(FIELD_FIELDS)) {
                fields.add(f.stringValue());
            }
            return fields.toArray(new String[0]);
        }
        public Number getFieldAsNumber (String field) {
            IndexableField f = doc.getField(field);
            if (f == null)
                throw new IllegalArgumentException ("Unknown field: "+field);
            return f.numericValue();
        }
        public String getFieldAsString (String field) {
            IndexableField f = doc.getField(field);
            if (f == null)
                throw new IllegalArgumentException ("Unknown field: "+field);
            return f.stringValue();
        }
    }
    static final Result POISON_RESULT = new Result ();

    public static class ResultEnumeration implements Enumeration<Result> {
        final BlockingQueue<Result> queue;
        Result next;
        
        ResultEnumeration (BlockingQueue<Result> queue) {
            this.queue = queue;
//            next ();
        }

        void next () {
            try {
                next = queue.take();
            }
            catch (Exception ex) {
                ex.printStackTrace();
                next = POISON_RESULT; // terminate
            }
        }           

        public boolean hasMoreElements () {
        	 if (next == null) next ();
            return next != POISON_RESULT;
        }
        
        public Result nextElement () {
        	if (next == null) next ();
            Result current = next;
            next ();
            return current;
        }
    }
    
    // it would have been faster to use a simple lookup table here?
    static int popcnt (byte[] b) {
        int c = 0;
        for (int i = 0; i < b.length; ++i) {
            c += Integer.bitCount(b[i] & 0xff);
        }
        return c;
    }

    static class ResultBlockingQueue<T>
        extends PriorityBlockingQueue<T> {

        final ReentrantLock lock = new ReentrantLock ();
        final Condition filled = lock.newCondition();

        private final int bufsiz;
        private boolean ready = false;

        public ResultBlockingQueue () {
            this (10000);
        }

        public ResultBlockingQueue (int bufsiz) {
            super (bufsiz);
            this.bufsiz = bufsiz;
            //System.err.println("#### "+getClass()+": bufsiz="+bufsiz);
        }

        @Override
        public void put (T r) {
            super.put(r);

            lock.lock();
            try {
                if (!ready && (size() >= bufsiz || r == POISON_RESULT)) {
                    /*
                    System.err.println("************ QUEUE IS NOW READY "
                                       +"FOR CONSUMPTION ("+size()
                                       +")! ************");
                    */
                    ready = true;
                    filled.signal();
                }
            }
            finally {
                lock.unlock();
            }
        }

        @Override
        public T take () throws InterruptedException {
            lock.lock();
            try {
                /*
                 * here we wait enough there are sufficient elements
                 * in the priority queue to have a meaningful ordering
                 * before the first take().
                 */
                while (!ready)
                    filled.await();
            }
            finally {
                lock.unlock();
            }

            return super.take();
        }
    }

    static class Tanimoto implements Callable<Integer> {
        final BlockingQueue<Payload> in;
        final BlockingQueue<Result> out;
        final int max;
        final double threshold;
        final Fingerprint query;

        Tanimoto (BlockingQueue<Payload> in,
                  BlockingQueue<Result> out,
                  Fingerprint query, int max,
                  double threshold) {
            this.in = in;
            this.out = out;
            this.max = max;
            this.threshold = threshold;
            this.query = query;
        }

        public Integer call () throws Exception {
            int count = 0;
            for (Payload p; (p = in.take()) != POISON_PAYLOAD
                     && (max <= 0 || (max > 0 && out.size() < max));) {
            
            	Fingerprint fp = p.getFpSim();
            	double similarity = fp.tanimotoSimilarity(query);
            	 if (similarity >= threshold) {
                     ++count;
                     out.put(new Result (p, similarity,null));
                 }
            	
            }
            return count;
        }
    }
    
    static class GraphIso implements Callable<Integer> {
        final BlockingQueue<Payload> in;
        final BlockingQueue<Result> out;
        final IsoMorphismSearcher isomorphismSearcher;
        
        final int max;
        final BitSet fpBits;
        final byte[] fp;
        final byte[] fpS;

        GraphIso (BlockingQueue<Payload> in,
                  BlockingQueue<Result> out,
                  IsoMorphismSearcher isomorphismSearcher, Fingerprint fp, int max, Fingerprint fpSim) {
            this.in = in;
            this.out = out;
            this.max = max;
            this.fp = fp.toByteArray();
            this.fpBits = fp.toBitSet();
            this.fpS = fpSim.toByteArray();
            this.isomorphismSearcher = isomorphismSearcher;
        }
        
        public Integer call () throws Exception {
            int count = 0;
            for (Payload p; (p = in.take()) != POISON_PAYLOAD
                     && (max <= 0 || (max > 0 && out.size() < max));) {
                byte[] pfp = p.getFpSub().toByteArray();
//               p.getFpSub()
                int i = 0;
                for (; i < fp.length; ++i) {
                    if ((pfp[i] & fp[i]) != fp[i]) {
                        break;
                    }
                }
                
                if (i == fp.length) {
                	byte[] pfpSim = p.getFpSim().toByteArray();
                	int a=0;
                	int b=0;
//                	p.getFpSim().tanimotoSimilarity()
                	for (i=0; i < pfpSim.length && i < fpS.length; ++i) {
                         a += Integer.bitCount(pfpSim[i] & fpS[i]);
                         b += Integer.bitCount(pfpSim[i] | fpS[i]);
                    }
                    if(pfpSim.length > fpS.length){
                	    for(; i < pfpSim.length; i++){
                	        b+= Integer.bitCount(pfpSim[i]);
                        }
                    }else{
                        for(; i < fpS.length; i++){
                            b+= Integer.bitCount(fpS[i]);
                        }
                    }
                	
                	int[] hits = isomorphismSearcher.findMax(p.getMol());
                	if(hits.length !=0){
                		 out.put(new Result (p, (double)a/b, hits));
                         ++count;
                	}
                	/*
                    if (hits != null) {
                        MolAtom[] atoms = mol.getAtomArray();
                        for (i = 0; i < hits.length; ++i) {
                        	
                            if (hits[i] >= 0){
                            	//by setting the atom map to something >0
                            	//we can later highlight that atom
                            	//when we render the molecules that hit
                            	//so that's all this does
                                atoms[hits[i]].setAtomMap(i+1);
                            }
                        }
                        out.put(new Result (p, (double)a/b, hits));
                        ++count;
                    }
                    */
                }
            }
            return count;
        }
    }

    static class Output implements Callable<Integer> {
        final BlockingQueue<Payload> in;
        final BlockingQueue<Result> out;
        final int max;
        Output (BlockingQueue<Payload> in,
                BlockingQueue<Result> out,
                int max) {
            this.in = in;
            this.out = out;
            this.max = max;
        }

        public Integer call () throws Exception {
            int count = 0;
            for (Payload p; (p = in.take()) != POISON_PAYLOAD
                     && (max <= 0 || (max > 0 && out.size() < max));) {
                out.put(new Result (p, null,null));
            }
            return count;
        }
    }

    private File baseDir;
    private Directory indexDir;
    private Directory facetDir;
    private Directory metaDir;
    private IndexWriter metaWriter;
    private IndexWriter indexWriter;
    private DirectoryReader indexReader;
    private DirectoryTaxonomyWriter facetWriter;
    private Analyzer indexAnalyzer;
    private FacetsConfig facetsConfig;
    private Codebook[] codebooks;
    
    private ExecutorService threadPool;
    private boolean localThreadPool = false;
    private final ConcurrentMap<IndexReader, Long>
        readerBuffer = new ConcurrentHashMap<IndexReader, Long>();

    private ScheduledExecutorService scheduledPool =
        Executors.newSingleThreadScheduledExecutor();
    private AtomicLong updatesSinceSaved = new AtomicLong (0);
    private AtomicLong lastModified = new AtomicLong (0);
    
    private  Fingerprinter fingerPrinterSub = Fingerprinters.getFingerprinter(FingerprintSpecification.PATH_BASED.create()
    																										.setLength(512)
    												);
    
    private Fingerprinter fingerPrinterSim = new ConcatFingerprinter()
    		 								.addFP(new ECFingerprinter(512,3,1),512)
    		 								.addFP(fingerPrinterSub,512)
    		 								.folded(512);
    
    
    

//    private Fingerprinter fingerPrinterSim = Fingerprinters.getFingerprinter(FingerprintSpecification.PATH_BASED.create()
//			.setLength(512)
//);
    

    public static StructureIndexer openReadOnly (File dir) throws IOException {
        return new StructureIndexer (dir);
    }
    
    public static StructureIndexer open (File dir) throws IOException {
        return new StructureIndexer (dir, false);
    }
    
    public StructureIndexer (File dir) throws IOException {
        this (dir, true);
    }
    
    public StructureIndexer (File dir, boolean readOnly) throws IOException {
        this (dir, readOnly, Executors.newCachedThreadPool());
        localThreadPool = true;
    }
    
    public StructureIndexer (File dir, boolean readOnly,
                             ExecutorService threadPool) throws IOException {
        if(!dir.exists()){
            dir.mkdirs();
        }
        if (!dir.isDirectory())
            throw new IllegalArgumentException ("Not a directory: "+dir);
        if (threadPool == null)
            throw new IllegalArgumentException ("Thread pool is null");
        
        this.threadPool = threadPool;
        this.baseDir = dir;

        indexAnalyzer = createIndexAnalyzer ();
        
        File index = new File (dir, "index");
        if (!index.exists())
            index.mkdirs();
        File meta = new File (dir, "codebook");
        if (!meta.exists())
            meta.mkdirs();
        File facet = new File (dir, "facet");
        if (!facet.exists())
            facet.mkdirs();

        indexDir = new NIOFSDirectory(index.toPath(), NoLockFactory.INSTANCE);
        metaDir = new NIOFSDirectory (meta.toPath(), NoLockFactory.INSTANCE);
        facetDir = new NIOFSDirectory
            (facet.toPath(), NoLockFactory.INSTANCE);

        if (!readOnly) {
            metaWriter = new IndexWriter
                (metaDir, new IndexWriterConfig
                 ( indexAnalyzer));
            indexWriter = new IndexWriter (indexDir, new IndexWriterConfig 
                                           ( indexAnalyzer));

            if (metaWriter.numDocs() == 0) {
                /*
                logger.info("No meta documents found; "
                            +"configuring a new set of codebooks...");
                */
                codebooks = new Codebook[CODEBOOKS];
                for (int i = 0; i < codebooks.length; ++i) {
                    codebooks[i] = new Codebook (32*FPSIZE);
                }
            }
            else {
                codebooks = load (DirectoryReader.open(metaWriter, true));
            }
            indexReader = DirectoryReader.open(indexWriter, true);
            facetWriter = new DirectoryTaxonomyWriter (facetDir);

            scheduledPool.scheduleAtFixedRate(new Runnable () {
                    public void run () {
                        flush ();
                        // release the IndexReader if it's been in the buffer
                        // for more than 5s
                        releaseReaders (5000l);
                    }
                }, 5, 2, TimeUnit.SECONDS);
        }
        else {
            codebooks = load (DirectoryReader.open(metaDir));
            indexReader = DirectoryReader.open(indexDir);           
        }
        /*
        for (Codebook cb : codebooks) {
            logger.info("Codebook "+cb);
        }
        */
        
        facetsConfig = new FacetsConfig ();
        facetsConfig.setMultiValued(FIELD_SOURCE, true);
        facetsConfig.setRequireDimCount(FIELD_SOURCE, true);            
    }

    static Analyzer createIndexAnalyzer () {
        Map<String, Analyzer> fields = new HashMap<String, Analyzer>();
        fields.put(FIELD_ID, new KeywordAnalyzer ());
        fields.put(FIELD_CODEBOOK, new KeywordAnalyzer ());
        fields.put(FIELD_FIELDS, new KeywordAnalyzer ());
        return  new PerFieldAnalyzerWrapper 
            (new StandardAnalyzer (), fields);
    }

    protected synchronized DirectoryReader getReader (boolean applyDeletes)
        throws IOException {
        DirectoryReader reader = indexWriter != null
            ? DirectoryReader.openIfChanged(indexReader, indexWriter, applyDeletes)
            : DirectoryReader.openIfChanged(indexReader);
        
        if (reader != null) {
            readerBuffer.putIfAbsent(indexReader, System.currentTimeMillis());
            indexReader = reader;
        }
        return indexReader;
    }

    protected IndexSearcher getIndexSearcher () throws IOException {
        return getIndexSearcher (true);
    }
    
    protected IndexSearcher getIndexSearcher (boolean applyDeletes)
        throws IOException {
        return new IndexSearcher (getReader (applyDeletes), threadPool);
    }

    protected void releaseReaders (long elapsed) {
        List<IndexReader> remove = new ArrayList<IndexReader>();
        for (Map.Entry<IndexReader, Long> me : readerBuffer.entrySet()) {
            if (System.currentTimeMillis() - me.getValue() > elapsed) {
                IndexReader reader = me.getKey();
                try {
                    reader.close();
                }
                catch (IOException ex) {
                    logger.warning("Can't close reader "+reader);
                }
                remove.add(reader);
            }
        }
        
        for (IndexReader r : remove)
            readerBuffer.remove(r);
    }

    public File getBasePath () { return baseDir; }
    public String[] getFields () {
        Set<String> fields = new TreeSet<String>();
        try {
            getFields (getReader (true), fields);
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        return fields.toArray(new String[0]);
    }

    protected void getFields (IndexReader reader, Set<String> fields)
        throws IOException {
        List<IndexReaderContext> children = reader.getContext().children();
        if (children == null) {
            for (LeafReaderContext ctx : reader.leaves()) {
                Terms terms = ctx.reader().terms(FIELD_FIELDS);
                //logger.info("terms "+terms);
                if (terms != null) {
                    TermsEnum en = terms.iterator();
                    for (BytesRef ref; (ref = en.next()) != null; ) {
                        fields.add(new String
                                   (ref.bytes, ref.offset, ref.length));
                    }
                }
            }
        }
        else {
            for (IndexReaderContext child : children)
                getFields (child.reader(), fields);
        }
    }
    
    public void shutdown () {
            scheduledPool.shutdown();
            if (metaWriter != null) {
                flush ();
                IOUtil.closeQuietly(metaWriter);
            }
            
            //closeQuietly does null check so we don't have to
            IOUtil.closeQuietly(indexReader);
            
            
            
            if(((NIOFSDirectory)facetDir).getDirectory().toFile().exists()) {
                IOUtil.closeQuietly(facetWriter);
            }
            
            for (IndexReader reader : readerBuffer.keySet()) {
                IOUtil.closeQuietly(reader);
            }

            IOUtil.closeQuietly(indexDir);
            IOUtil.closeQuietly(metaDir);

            if(((NIOFSDirectory)facetDir).getDirectory().toFile().exists()) {
                IOUtil.closeQuietly(facetDir);
            }
            if (localThreadPool) {
                threadPool.shutdown();
            }

    }

    protected void update (Codebook cb) throws IOException {
        Term term = new Term (FIELD_ID, cb.getName());
        //logger.info("Persist codebook "+cb);    
        metaWriter.updateDocument(term, Codebook.instrument(cb));
    }

    synchronized void flush () {
        if (updatesSinceSaved.get() > 0) {
            /*
            logger.info("### flushing "
                        +updatesSinceSaved.get()+" updates...");
            */
            try {
                for (Codebook cb : codebooks) {
                    update (cb);
                }
                metaWriter.commit();
                indexWriter.commit();
                facetWriter.commit();
                updatesSinceSaved.set(0);
            }catch (IOException ex) {
//                ex.printStackTrace();
            }
        }
    }

    protected Codebook[] load (DirectoryReader reader) throws IOException {
        try {
            Codebook[] cbooks = new Codebook[reader.numDocs()];
            for (int i = 0; i < reader.numDocs(); ++i) {
                Document doc = reader.document(i);
                cbooks[i] = Codebook.load(doc);
                //logger.info("Loading doc "+i+" "+cbooks[i]);
            }
            return cbooks;
        }
        finally {
            reader.close();
        }
    }

    public int size () {
        int size;
        if (indexWriter != null) {
            if (indexWriter.hasUncommittedChanges()) {
                try {
                    indexWriter.commit();
                }
                catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
            size = indexWriter.numDocs();
        }
        else {
            size = indexReader.numDocs();
        }
        
        return size;
    }

    public Map<String, Integer> getSources () throws IOException {
        IndexSearcher searcher = getIndexSearcher ();
        return getSources (searcher);
    }
    
    protected Map<String, Integer> getSources (IndexSearcher searcher)
        throws IOException {
        Map<String, Integer> map = new TreeMap<String, Integer>();
        try {
            FacetsCollector fc = new FacetsCollector ();
            TaxonomyReader taxon = new DirectoryTaxonomyReader (facetWriter);
            TopDocs hits = FacetsCollector.search
                (searcher, new MatchAllDocsQuery (),
                 searcher.getIndexReader().numDocs(), fc);
            Facets facets = new FastTaxonomyFacetCounts
                    (taxon, facetsConfig, fc);
            List<FacetResult> facetResults = facets.getAllDims(20);
            for (FacetResult result : facetResults) {
                if (result.dim.equals(FIELD_SOURCE)) {
                    for (int i = 0; i < result.labelValues.length; ++i) {
                        LabelAndValue lv = result.labelValues[i];
                        map.put(lv.label, lv.value.intValue());
                    }
                }
            }
            taxon.close();
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        return map;
    }

    public void add (String id, String struc) throws IOException {
        add (null, id, struc);
    }
    
    public void add (String source, String id, String struc)
        throws IOException {
        try {
           
			add (source, id,  Chemical.parse(struc));
        }
        catch (Exception ex) {
            throw new IllegalArgumentException ("Bogus molecule format", ex);
        }
    }
    
    public void add (String id, Chemical struc) throws IOException {
        add (null, id, struc);
    }
    
    public void add (String source, String id, Chemical struc)
        throws IOException {
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");
        
        Document doc = new Document ();
        doc.add(new StringField (FIELD_ID, id, YES));
        if (source != null) {
            doc.add(new FacetField (FIELD_SOURCE, source));
            doc.add(new StringField (FIELD_SOURCE, source, YES));
            doc.add(new TextField (FIELD_TEXT, source, NO));
        }
        doc.add(new TextField (FIELD_TEXT, id, NO));
        instrument (doc, struc);
        doc = facetsConfig.build(facetWriter, doc);
        indexWriter.addDocument(doc);
        updated ();
    }

    protected void updated () throws IOException {
        updatesSinceSaved.incrementAndGet();
        lastModified.set(System.currentTimeMillis());
    }

    protected void instrument (Document doc, Chemical orig)
        throws IOException {
       
       Chemical chemical = orig.copy();
       
       chemical.makeHydrogensImplicit();
       
       
        try{
            processQuery(chemical);
        }catch(Exception e){
//            e.printStackTrace();
            //ignore?
        }
        
		Fingerprint fingerprintSub = fingerPrinterSub.computeFingerprint(chemical);
		byte[] fp =  fingerprintSub.toByteArray();
		Fingerprint fingerprintSim = fingerPrinterSim.computeFingerprint(chemical);
		byte[] fpSim =  fingerprintSim.toByteArray();
		
		chemical.makeHydrogensExplicit();
		String indexMolHExp = chemical.toMol(new ChemFormat.MolFormatSpecification()
                .setKekulization(ChemFormat.KekulizationEncoding.FORCE_AROMATIC));
                
		
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            int code = cb.encode(fingerprintSub);
//            System.out.println("code book " + i + " = " + code);
            if (code != 0) {
                cb.incr(code); // this must be in-sync with the document count!
                String encodedCB = cb.encode(code);
                doc.add(new StringField
                        (FIELD_CODEBOOK,encodedCB , NO));
            }
        }
        for(Entry<String, String> entry : chemical.getProperties().entrySet()){
            String prop = entry.getKey();
            String value = entry.getValue();
            if (value != null) {
                doc.add(new TextField (FIELD_TEXT, value, NO));
                doc.add(new TextField (FIELD_TEXT, prop, NO));
                doc.add(new TextField (FIELD_FIELDS, prop, YES));
                try {
                    double dv = Double.parseDouble(value);
                    doc.add(new DoubleField (prop, dv, NO));
                }
                catch (NumberFormatException ex) {
                }
                try {
                    long lv = Long.parseLong(value);
                    doc.add(new LongField (prop, lv, NO));
                }
                catch (NumberFormatException ex) {
                }
                try {
                    int iv = Integer.parseInt(value);
                    doc.add(new IntField (prop, iv, NO));
                }
                catch (NumberFormatException ex) {
                }
                doc.add(new TextField (prop, value, YES));
            }
        }
        String name = chemical.getName();
        if (name != null && name.length() > 0) { 
            doc.add(new TextField (FIELD_NAME, name, NO));
            doc.add(new TextField (FIELD_TEXT, name, NO));
        }
        doc.add(new StoredField (FIELD_FINGERPRINT_SUB, fp));
        doc.add(new StoredField (FIELD_FINGERPRINT_SIM, fpSim));
        
        doc.add(new IntField (FIELD_POPCNT, popcnt (fpSim), NO));
        
        
       
        doc.add(new StoredField(FIELD_MOLFILE, indexMolHExp));
        
       
       doc.add(new StringField (FIELD_FORMULA, chemical.getFormula(), YES));
        doc.add(new IntField (FIELD_NATOMS, chemical.getAtomCount(), NO));
        doc.add(new IntField (FIELD_NBONDS, chemical.getBondCount(), NO));
        doc.add(new DoubleField (FIELD_MOLWT,chemical.getMass(), NO));
    }

    private static int getLineOfDefline(String mol) throws IOException{
        try(BufferedReader reader = new BufferedReader(new StringReader(mol))){
            int i=1;
            String line=null;
            while((line = reader.readLine()) !=null && !line.endsWith("V2000")){
                i++;
            }
            return i;
        }
    }
    
    public void remove (String source, String id) throws IOException {
        remove (getIndexSearcher (false), source, id);
    }
    
    protected void remove (IndexSearcher searcher,
                           String source, String id) throws IOException {
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");
        
        /* TODO: this method should only be one line but because we have 
         * to keep the document count in sync with the codebooks, we need
         * to do this..
         */
        Query q = null;
        if (source != null)
            q = new TermQuery (new Term (FIELD_SOURCE, source));

        if (id != null) {
            TermQuery tq = new TermQuery (new Term (FIELD_ID, id));
            if (q != null) {
                BooleanQuery bq = new BooleanQuery ();
                bq.add(q, BooleanClause.Occur.MUST);
                bq.add(tq, BooleanClause.Occur.MUST);
                q = bq;
            }
            else
                q = tq;
        }

        if (q == null)
            throw new IllegalArgumentException
                ("Either source or id must be specified!");
            
        int total = searcher.getIndexReader().numDocs();
        TopDocs hits = searcher.search(q, total);
        
        logger.info("Deleting "+id+" ["+source+"].."+hits.totalHits
                    +"/"+total);
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            for (String code : doc.getValues(FIELD_CODEBOOK))
                for (Codebook cb : codebooks) {
                    int c = cb.decode(code);
                    if (c >= 0)
                        cb.decr(c);
                }
        }
        
        indexWriter.deleteDocuments(q);
        updated ();
    }

    public void remove (String source) throws IOException {
        // this potentially can dramatically alter the counts, so we need
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");

        indexWriter.deleteDocuments(new Term (FIELD_SOURCE, source));
        updated ();
        
        // do update in background
        threadPool.submit(new Runnable () {
                public void run () {
                    try {
                        IndexSearcher searcher = getIndexSearcher ();
                        for (Codebook cb : codebooks) {
                            cb.adjustCounts(searcher);
                        }
                    }
                    catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            });
    }

    public Codebook[] getCodebooks () { return codebooks; }
    public long lastModified () { return lastModified.get(); }
    
    public ResultEnumeration substructure (String query, Filter... filters)
        throws Exception {
        return substructure (query, -1, 2, filters);
    }

    public ResultEnumeration substructure
        (String query, int max) throws Exception {
        return substructure (query, max, 2, null);
    }

    public ResultEnumeration substructure
        (String query, int max, int nthreads, Filter... filters)
        throws Exception {
    	//query string could be a mol or a smiles use reader
       Chemical chemical = Chemical.parse(query);
        //GSRS-1095 check for aromatized bonds
        //katzelda Sept 21 don't need to aromatize when using CDK
//        boolean present = chemical.bonds().filter(Bond::isAromatic).findAny().isPresent();
////
//        if(!present) {
//            boolean hasQueryBonds = chemical.bonds().filter(Bond::isQueryBond).findAny().isPresent();
//            if(!hasQueryBonds) {
////            System.out.println("AROMATIZING INPUT ");
//                chemical.aromatize();
//            }
//        }
//        chemical.removeNonDescriptHydrogens();
       
       //I don't think this is necessary
//        if(!chemical.hasCoordinates()){
//            chemical.generateCoordinates();
//        }
        return substructure (chemical, max, nthreads, filters);
    }

    public ResultEnumeration substructure
        (Chemical query) throws Exception {
        return substructure (query, -1, 2, null);
    }
    
    public ResultEnumeration substructure
        (Chemical query, Filter... filters) throws Exception {
        return substructure (query, -1, 2, filters);
    }

    public ResultEnumeration substructure
        (Chemical query, final int max, int nthreads) throws Exception {
        return substructure (getIndexSearcher (), query, max, nthreads, null);
    }
    
    public ResultEnumeration substructure
        (Chemical query, final int max, int nthreads, Filter... filters)
        throws Exception {
        return substructure (getIndexSearcher (),
                             query, max, nthreads, filters);
    }
    
    protected void processQuery(Chemical query) {
        //We do want to aromatize the query, but we DO NOT want to DE-aromatize parts that aren't aromatized
        Set<Bond> abs=query.bonds()
                           .filter(bb->bb.isAromatic())
                           .collect(Collectors.toSet());
        
        query.aromatize();
        
        abs.stream()
           .filter(bb->!bb.isAromatic())
           .forEach(bb->bb.setBondType(BondType.AROMATIC));
        
    }
    protected Chemical processQueryForFP(Chemical query) {
        Chemical copyr=query.copy();
        copyr.makeHydrogensImplicit();
//        copyr.atoms().filter(aa->"H".equals(aa.getSymbol())).collect(Collectors.toList())
//        .forEach(aa->{
//            copyr.removeAtom(aa);
//        });
        return copyr;
    }
    
    protected ResultEnumeration substructure
        (IndexSearcher searcher, Chemical query,
         final int max, int nthreads, Filter... filters) throws Exception {

        processQuery(query);        
        Chemical copyr=processQueryForFP(query);
        
        Fingerprint qfp = fingerPrinterSub.computeFingerprint(copyr);
        Fingerprint qfpSim = fingerPrinterSim.computeFingerprint(copyr);
        
        
        Codebook bestCb = null;
        int bestHits = Integer.MAX_VALUE;
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            int[] eqv = cb.apply(qfp);
            if (eqv != null) {
                int hits = 0;
                for (int j = 0; j < eqv.length; ++j) {
                    hits += cb.count(eqv[j]);
                }
//                if(hits !=0){
//                    System.out.println("code book "+ i + " hits = " + hits);
//                }
                if (/*hits !=0 &&*/ hits < bestHits) {
                    bestHits = hits;
                    bestCb = cb;
                }
            }
        }
        
        Query q = null;
        if (bestCb == null) {
            // iterate over all documents
            q = new MatchAllDocsQuery ();
        }
        else {
            DisjunctionMaxQuery mq = new DisjunctionMaxQuery (1.f);
            int[] eqv = bestCb.apply(qfp);
//            System.out.println("eqv = " + Arrays.toString(eqv));
            for (int j = 0; j < eqv.length; ++j) {
                TermQuery tq = new TermQuery
                    (new Term (FIELD_CODEBOOK, bestCb.encode(eqv[j])));
                mq.add(tq);
            }
            q = mq;
        }
        
        int total = searcher.getIndexReader().numDocs();
        long start = System.currentTimeMillis();
        if (filters != null) {
            for (Filter f : filters) {
                q = new FilteredQuery
                    (q, f, FilteredQuery.QUERY_FIRST_FILTER_STRATEGY);
            }
        }
        TopDocs hits = searcher.search(q, total);
        logger.info("## total screened: "+(total-hits.totalHits)+"/"+total
                    +" screen efficiency: "
                    +String.format("%1$.2f", 1.-(double)hits.totalHits/total)
                    +" ellapsed: "
                    +String.format("%1$.2fs",
                                   (System.currentTimeMillis()-start)*1e-3));
        
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final BlockingQueue<Result> out = new ResultBlockingQueue<>();
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>();
        for (int i = 0; i < nthreads; ++i){
        	IsoMorphismSearcher chemSearcher = new IsoMorphismSearcher(query);
            threads.add(threadPool.submit
                        (new GraphIso (in, out, chemSearcher, qfp, max,qfpSim)));
        }
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (hits.scoreDocs[i].doc, doc));
        }
        for (int i = 0; i < nthreads; ++i)
            in.put(POISON_PAYLOAD);

        threadPool.submit(new Runnable () {
                public void run () {
                    try {                   
                        int total = 0;
                        for (Future<Integer> f : threads) {
                        	try{
                        		total += f.get();
                        	}catch(Throwable ex){
                        		//catch exception inside loop so
                        		//not only do we keep getting the others
                        		//but we eventually add the POISON_RESULT
                        		//to the queue.
                        		//otherwise, an exception will make this loop forever
                        		ex.printStackTrace();
                        	}
                        }
                        out.put(POISON_RESULT);
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    public ResultEnumeration similarity
        (String query, double threshold, int max) throws Exception {
    	Chemical chemical = Chemical.parse(query);
        return similarity (chemical, threshold, max, 2);
    }
    
    public ResultEnumeration similarity (String query, double threshold)
        throws Exception {
    	Chemical chemical = Chemical.parse(query);
        return similarity (chemical, threshold, -1, 2);
    }

    public ResultEnumeration similarity
        (String query, double threshold, Filter... filters) throws Exception {
        Chemical chemical = Chemical.parse(query);
        return similarity (chemical, threshold, filters);
    }

    public ResultEnumeration similarity
        (Chemical query, double threshold, Filter... filters) throws Exception {
        return similarity (query, threshold, -1, 2, filters);
    }
    
    public ResultEnumeration similarity
        (Chemical query, final double threshold, 
         final int max, final int nthreads) throws Exception {
        return similarity (query, threshold, max, nthreads, null);
    }
    
    public ResultEnumeration similarity
        (Chemical query, final double threshold, 
         final int max, final int nthreads, Filter... filters)
        throws Exception {
        return similarity
            (getIndexSearcher (), query, threshold, max, nthreads, filters);
    }

    protected ResultEnumeration similarity
        (IndexSearcher searcher, Chemical query, final double threshold,
         final int max, final int nthreads, Filter... filters)
        throws Exception {
        /*
         * first calculate the minimum popcnt needed to satisfy the
         * cutoff:
         *   (q & t)/(q | t) <= threshold <= min(|q|,|t|)/max(|q|,|t|)
         * where |q| and |t| are the popcnt of the query and target, 
         * respectively. This in turn provides the following lower bound:
         *   q*threshold <= t
         *   
         * TP: That's just a lower bound. There's an upper bound too.
         * 
         *  To find that, you just consider the case where there's a FP with every bit
         *  that's in the query, and then a lot of other ones too. That is, the target is
         *  a feature-based superstructure of the query. That would give the
         *  following simple formula:
         *  
         *  tanimoto <= |q|/|t| 
         *  
         *  Which means:
         *  
         *  |t| <= |q|/tanimoto
         *  
         *  So, the upper bound pop count for target is pop(q)/cutoff. Note that this assumes that
         *  the target is a super structure. If, instead, the target has some missing stuff from
         *  the query too:
         *  
         *  tanimoto <= (|q|-miss)/(|t|+miss)
         *  
         *  |t|+miss <= (|q|-miss)/tanimoto
         *  |t| <= (|q|-miss)/tanimoto - miss <= |q|/tanimoto 
         *  
         *  Therefore |q|/tanimoto is a good pessimistic upperbound.
         *  
         */
        processQuery(query);        
        Chemical copyr=processQueryForFP(query);
        Fingerprint q = fingerPrinterSim.computeFingerprint(copyr);
        int popcnt = q.populationCount();

        int minpop = (int)(popcnt*threshold+0.5);
        int maxpop = (int)(popcnt*(1.0/threshold)+0.5);
        Query range = NumericRangeQuery.newIntRange
            (FIELD_POPCNT, minpop, maxpop, true, true);
        if (filters != null) {
            for (Filter f : filters)
                range = new FilteredQuery
                    (range, f, FilteredQuery.QUERY_FIRST_FILTER_STRATEGY);
        }
        long start = System.currentTimeMillis();
        TopDocs hits = searcher.search
            (range, searcher.getIndexReader().numDocs());
        logger.info("## range query >="+minpop
                    +": "+hits.totalHits+" ellapsed: "
                    +String.format("%1$.2fs",
                                   (System.currentTimeMillis()-start)*1e-3));
        
        final BlockingQueue<Result> out = new ResultBlockingQueue<>();
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>(); 
        for (int i = 0; i < nthreads; ++i)
            threads.add(threadPool.submit
                        (new Tanimoto (in, out, q, max, threshold)));
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (hits.scoreDocs[i].doc, doc));
        }

        for (int i = 0; i < nthreads; ++i)
            in.put(POISON_PAYLOAD);
        threadPool.submit(new Runnable () {
                public void run () {
                    try {
                        int total = 0;
                        for (Future<Integer> f : threads) {
                            total += f.get();
                        }

                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }finally{
                        //issue 3 put poison payload in finally block
                        // in case try block above errors out
                        try {
                            out.put(POISON_RESULT);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    public ResultEnumeration search (Query query) throws Exception {
        return search (query, 0, null);
    }

    public ResultEnumeration search (Filter... filters) throws Exception {
        return search (new MatchAllDocsQuery (), 0, filters);
    }

    public ResultEnumeration search (String query, int max, Filter... filters)
        throws Exception {
        QueryParser parser = new QueryParser (FIELD_TEXT, indexAnalyzer);
        return search (parser.parse(query), max, filters);
    }
    
    public ResultEnumeration search (Query query, int max, Filter... filters)
        throws Exception {
        return search (getIndexSearcher (), query, max, filters);
    }

    protected ResultEnumeration search
        (IndexSearcher searcher, Query query, int max, Filter... filters)
        throws Exception {

        if (query == null && filters == null)
            throw new IllegalArgumentException
                ("Both query and filter are null!");
        

        final BlockingQueue<Result> out = new LinkedBlockingQueue<Result>();
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final Future<Integer> future =
            threadPool.submit(new Output(in, out, max));
        
        if (max <= 0)
            max = searcher.getIndexReader().numDocs();

        if (filters != null) {
            for (Filter f : filters)
                query = new FilteredQuery (query, f);
        }
        
        TopDocs hits = searcher.search(query, max);
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (hits.scoreDocs[i].doc, doc));
        }
        in.put(POISON_PAYLOAD);
        
        threadPool.submit(new Runnable () {
                public void run () {
                    try {
                        int total = future.get();
                        out.put(POISON_RESULT);
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    public ResultEnumeration search (String query) throws Exception {
        return search (query, 0);
    }
    
    public ResultEnumeration search (String query, final int max)
        throws Exception {
        return search (query, max, 1);
    }
    
    public ResultEnumeration search (String query, final int max,
                                     final int nthreads) throws Exception {
        return search (getIndexSearcher (), query, max, nthreads);
    }

    protected ResultEnumeration search
        (final IndexSearcher searcher, String query,
         final int max, final int nthreads) throws Exception {
        QueryParser parser = new QueryParser (FIELD_TEXT, indexAnalyzer);
        Query q = parser.parse(query);
        long start = System.currentTimeMillis();
		TopDocs hits = searcher.search(q, searcher.getIndexReader().numDocs());
        logger.info("## query "+q
                    +" yields "+hits.totalHits+" ellapsed: "
                    +String.format("%1$.2fs",
                                   (System.currentTimeMillis()-start)*1e-3));
        
        final BlockingQueue<Result> out = new LinkedBlockingQueue<Result>();
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>();
        for (int i = 0; i < nthreads; ++i)
            threads.add(threadPool.submit(new Output(in, out, max)));
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (hits.scoreDocs[i].doc, doc));
        }

        for (int i = 0; i < nthreads; ++i)
            in.put(POISON_PAYLOAD);
        
        threadPool.submit(new Runnable () {
                public void run () {
                    try {
                        int total = 0;
                        for (Future<Integer> f : threads) {
                            total += f.get();
                        }
                        out.put(POISON_RESULT);
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    protected int[] histogram (IndexSearcher searcher,
                               String field, int[] range) throws IOException {
        if (range == null)
            throw new IllegalArgumentException ("Range is null");
        
        int[] hist = new int[range.length+2];
        int numDocs = searcher.getIndexReader().numDocs();
        if (numDocs > 0) {
            Query query = NumericRangeQuery.newIntRange
                (field, null, range[0], false, false);
            TopDocs hits = searcher.search(query, numDocs);
            hist[0] = hits.totalHits;
            for (int i = 1; i < range.length; ++i) {
                query = NumericRangeQuery.newIntRange
                    (field, range[i-1], range[i], true, false);
                hits = searcher.search(query, numDocs);
                hist[i] = hits.totalHits;
            }
            query = NumericRangeQuery.newIntRange
                (field, range[range.length-1], null, true, false);
            hits = searcher.search(query, numDocs);
            hist[range.length] = hits.totalHits;
        }
       
        
        return hist;
    }

    protected int[] histogram (IndexSearcher searcher,
                               String field, long[] range) throws IOException {
        if (range == null)
            throw new IllegalArgumentException ("Range is null");
        
        int[] hist = new int[range.length+2];
        int numDocs = searcher.getIndexReader().numDocs();
        if (numDocs > 0) {
            Query query = NumericRangeQuery.newLongRange
                (field, null, range[0], false, false);
            TopDocs hits = searcher.search(query, numDocs);
            hist[0] = hits.totalHits;
            for (int i = 1; i < range.length; ++i) {
                query = NumericRangeQuery.newLongRange
                    (field, range[i-1], range[i], true, false);
                hits = searcher.search(query, numDocs);
                hist[i] = hits.totalHits;
            }
            query = NumericRangeQuery.newLongRange
                (field, range[range.length-1], null, true, false);
            hits = searcher.search(query, numDocs);
            hist[range.length] = hits.totalHits;
        }
       
        
        return hist;
    }
    
    protected int[] histogram (IndexSearcher searcher,
                               String field, double[] range)
        throws IOException {
        if (range == null)
            throw new IllegalArgumentException ("Range is null");
        
        int[] hist = new int[range.length+2];
        int numDocs = searcher.getIndexReader().numDocs();
        if (numDocs > 0) {
            Query query = NumericRangeQuery.newDoubleRange
                (field, null, range[0], false, false);
            TopDocs hits = searcher.search(query, numDocs);
            hist[0] = hits.totalHits;
            for (int i = 1; i < range.length; ++i) {
                query = NumericRangeQuery.newDoubleRange
                    (field, range[i-1], range[i], true, false);
                hits = searcher.search(query, numDocs);
                hist[i] = hits.totalHits;
            }
            query = NumericRangeQuery.newDoubleRange
                (field, range[range.length-1], null, true, false);
            hits = searcher.search(query, numDocs);
            hist[range.length] = hits.totalHits;
        }
        
        return hist;
    }

    public int[] histogram (String field, double[] range) throws IOException {
        return histogram (getIndexSearcher (), field, range);
    }
    public int[] histogram (String field, int[] range) throws IOException {
        return histogram (getIndexSearcher (), field, range);
    }
    public int[] histogram (String field, long[] range) throws IOException {
        return histogram (getIndexSearcher (), field, range);
    }
    
    protected void statPopcnt (IndexSearcher searcher, PrintStream ps)
        throws IOException {
        int[] range = new int[]{0, 50, 100, 150, 200, 250, 300, 350};
        int[] hist = histogram (searcher, FIELD_POPCNT, range);
        ps.println("Popcnt histogram:");
        for (int i = 1; i < range.length; ++i) {
            ps.println(String.format("[%1$3d,%2$3d): %3$d",
                                     range[i-1], range[i], hist[i]));
        }
        ps.println(String.format("    >%1$3d:", range[range.length-1])
                       +"  "+hist[range.length]);
    }

    protected void statAtoms (IndexSearcher searcher, PrintStream ps)
        throws IOException {
        int[] range = new int[]{0, 10, 20, 30, 40, 50, 100};
        int[] hist = histogram (searcher, FIELD_NATOMS, range);
        ps.println("Atom count histogram:");
        for (int i = 1; i < range.length; ++i) {
            ps.println(String.format("[%1$3d,%2$3d): %3$d",
                                     range[i-1], range[i], hist[i]));
        }
        ps.println(String.format("    >%1$3d:", range[range.length-1])
                       +"  "+hist[range.length]);
    }
    
    protected void statMolwt (IndexSearcher searcher, PrintStream ps)
        throws IOException {
        double[] range = new double[]{100., 200., 350., 500., 700.};
        int[] hist = histogram (searcher, FIELD_MOLWT, range);
        ps.println("Molwt histogram:");
        ps.println(String.format("    <%1$3.0f :", range[0])+" "+hist[0]);
        for (int i = 1; i < range.length; ++i) {
            ps.println(String.format("[%1$3.0f,%2$3.0f): %3$d",
                                     range[i-1], range[i], hist[i]));
        }
        ps.println(String.format("    >%1$3.0f :", range[range.length-1])
                       +" "+hist[range.length]);
    }
    
    protected void stats (IndexSearcher searcher, PrintStream ps)
        throws IOException {
        ps.println("[*** stats for "+baseDir+" ***]");
        statPopcnt (searcher, ps);
        statMolwt (searcher, ps);
        statAtoms (searcher, ps);
    }
    
    public void stats (PrintStream ps) throws IOException {
        stats (getIndexSearcher (), ps);
    }
}
