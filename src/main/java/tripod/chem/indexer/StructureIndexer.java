package tripod.chem.indexer;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;
import java.util.concurrent.locks.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.apache.lucene.document.Document;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.LongField;
import org.apache.lucene.document.IntField;
import org.apache.lucene.document.FloatField;
import org.apache.lucene.document.DoubleField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.FieldType;
import org.apache.lucene.document.IntDocValuesField;
import org.apache.lucene.document.NumericDocValuesField;
import org.apache.lucene.document.FloatDocValuesField;
import org.apache.lucene.document.DoubleDocValuesField;
import static org.apache.lucene.document.Field.Store.*;

import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.SerialMergeScheduler;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.FieldInfo.IndexOptions;

import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.store.NIOFSDirectory;
import org.apache.lucene.store.NoLockFactory;
import org.apache.lucene.util.Version;
import org.apache.lucene.util.BytesRef;

import org.apache.lucene.analysis.core.KeywordAnalyzer;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;
import org.apache.lucene.analysis.miscellaneous.PerFieldAnalyzerWrapper;

import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.DisjunctionMaxQuery;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.RegexpQuery;
import org.apache.lucene.search.SortField;
import org.apache.lucene.search.Sort;
import org.apache.lucene.search.Filter;
import org.apache.lucene.search.MatchAllDocsQuery;
import org.apache.lucene.search.FieldCacheTermsFilter;
import org.apache.lucene.search.NumericRangeQuery;

import org.apache.lucene.queryparser.flexible.standard.StandardQueryParser;
import org.apache.lucene.queryparser.flexible.core.QueryNodeException;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.queryparser.classic.MultiFieldQueryParser;
import org.apache.lucene.queryparser.classic.ParseException;

import org.apache.lucene.facet.*;
import org.apache.lucene.facet.range.*;
import org.apache.lucene.facet.taxonomy.*;
import org.apache.lucene.facet.taxonomy.directory.*;
import org.apache.lucene.facet.sortedset.*;

import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.sss.search.MolSearch;

public class StructureIndexer {
    static final Logger logger =
        Logger.getLogger(StructureIndexer.class.getName());

    static final Molecule POISON_PILL = new Molecule ();
    static final Version LUCENE_VERSION = Version.LATEST;
    /**
     * Fields for each document
     */
    static final String FIELD_ID = "_id";
    static final String FIELD_SOURCE = "_source";
    static final String FIELD_CODEBOOK = "_codebook";
    static final String FIELD_FINGERPRINT = "_fingerprint";
    static final String FIELD_POPCNT = "_popcount"; // fingerprint pop count
    static final String FIELD_MOLFILE = "_molfile";
    static final String FIELD_MOLWT = "_molwt";
    static final String FIELD_NATOMS = "_natoms";
    static final String FIELD_NBONDS = "_nbonds";
    
    static final int FPSIZE = 16;
    static final int FPBITS = 2;
    static final int FPDEPTH = 6;
    
    static final int CODESIZE = 8; // 8-bit or 256
    static final int CODEBOOKS = 256;

    static final String CONFIG_FILE = "indexer.xml";

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
        CODEWORDS = codes.toArray(new String[0]);
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

        public static Codebook create (String name, Properties props) {
            String value = props.getProperty(name+".dict");
            if (value == null)
                throw new IllegalArgumentException
                    ("Codebook \""+name+"\" has no dict property!");
            String[] d = value.split(",");
            int[] dict = new int[d.length];
            for (int i = 0; i < d.length; ++i)
                dict[i] = Integer.parseInt(d[i]);
            Codebook cb = new Codebook (name, dict);
            value = props.getProperty(name+".counts");
            if (value != null) {
                d = value.split(",");
                if (d.length != cb.counts.length)
                    throw new IllegalArgumentException
                        ("Mismatch dictionary length; expecting "
                         +cb.counts.length+" but got "+d.length+"!");
                for (int i = 0; i < d.length; ++i)
                    cb.counts[i] = Integer.parseInt(d[i]);
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

        public int[] apply (int[] fp) {
            int code = encode (fp);
            return code == 0 ? null : eqv[code];
        }
        
        public int[] apply (byte[] fp) {
            int code = encode (fp);
            return code == 0 ? null : eqv[code];
        }

        public int count (int code) {
            return counts[code];
        }
        
        public int incr (int code) {
            return ++counts[code];
        }

        public int decr (int code) {
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
        
        /**
         * an arbitrary encoding.. 
         */
        public int encode (int[] fp) { // fingerprint bits stored as ints
            int code = 0;
            for (int i = 0; i < dict.length; ++i) {
                if (get (fp, dict[i]))
                    code |= 1<<i;
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
            StringBuilder sb = new StringBuilder ("{");
            for (int i = 0; i < dict.length; ++i) {
                if (i > 0) sb.append(",");
                sb.append(dict[i]);
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
        final Document doc;
        Molecule mol;
        byte[] fp;

        Payload (Document doc) {
            this.doc = doc;
        }
        
        Payload () {
            doc = null;
            mol = null;
            fp = null;
        }

        byte[] getFp () {
            if (fp == null) {
                BytesRef ref = doc.getBinaryValue(FIELD_FINGERPRINT);
                fp = ref.bytes;
            }
            return fp;
        }
        Molecule getMol () {
            if (mol == null) {
                BytesRef molref = doc.getBinaryValue(FIELD_MOLFILE);
                try {
                    MolHandler mh = new MolHandler
                        (new String (molref.bytes,
                                     molref.offset, molref.length));
                    mol = mh.getMolecule();
                    mol.setName(doc.get(FIELD_ID));
                }
                catch (Exception ex) {
                    throw new RuntimeException
                        ("Document "+doc.get(FIELD_ID)+" contains bogus "
                         +"field "+FIELD_MOLFILE+"!");
                }
            }
            return mol;
        }
    }
    static final Payload POISON_PAYLOAD = new Payload ();

    static public class Result {
        final Document doc;
        final Molecule mol;
        
        Result (Document doc, Molecule mol) {
            this.doc = doc;
            this.mol = mol;
        }
        Result () {
            this (null, null);
        }

        public String getSource () {
            return doc.get(FIELD_SOURCE);
        }
        public Molecule getMol () { return mol; }
        public Document getDoc () { return doc; }
    }
    static final Result POISON_RESULT = new Result ();

    public static class ResultEnumeration implements Enumeration<Result> {
        final BlockingQueue<Result> queue;
        Result next;
        
        ResultEnumeration (BlockingQueue<Result> queue) {
            this.queue = queue;
        }

        public boolean hasMoreElements () {
            try {
                next = queue.take();
                return next != POISON_RESULT;
            }
            catch (Exception ex) {
                ex.printStackTrace();
            }
            return false;
        }
        
        public Result nextElement () { return next; }
    }
    
    static String encodeDocId (String source, String id) {
        return source+"::"+id;
    }
    
    static String[] decodeDocId (String id) {
        String[] toks = id.split("::");
        if (toks.length != 2) {
            throw new IllegalArgumentException ("Bogus doc id: "+id);
        }
        return toks;
    }

    // it would have been faster to use a simple lookup table here?
    static int popcnt (byte[] b) {
        int c = 0;
        for (int i = 0; i < b.length; ++i)
            c += Integer.bitCount(b[i] & 0xff);
        return c;
    }

    static class Tanimoto implements Callable<Integer> {
        final BlockingQueue<Payload> in;
        final BlockingQueue<Result> out;
        final int max;
        final double threshold;
        final byte[] query;

        Tanimoto (BlockingQueue<Payload> in,
                  BlockingQueue<Result> out,
                  byte[] query, int max,
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
                     && (max < 0 || (max > 0 && out.size() < max));) {
                int a = 0, b = 0;
                byte[] fp = p.getFp();
                for (int j = 0; j < query.length; ++j) {
                    a += Integer.bitCount(query[j] & fp[j]);
                    b += Integer.bitCount(query[j] | fp[j]);
                }
            
                double tan = (double)a/b;
                if (tan >= threshold) {
                    ++count;
                    Molecule mol = p.getMol();
                    mol.setProperty("TANIMOTO", String.format("%1$.4f", tan));
                    out.put(new Result (p.doc, mol));
                }
            }
            logger.info(Thread.currentThread().getName()+" "
                        +count+" passed tanimoto cutoff "+threshold+"!");
            return count;
        }
    }
    
    static class GraphIso implements Callable<Integer> {
        final BlockingQueue<Payload> in;
        final BlockingQueue<Result> out;
        final MolSearch msearch;
        final int max;
        final byte[] fp;

        GraphIso (BlockingQueue<Payload> in,
                  BlockingQueue<Result> out,
                  Molecule query, byte[] fp, int max) {
            this.in = in;
            this.out = out;
            this.max = max;
            this.fp = fp;
            msearch = new MolSearch ();
            msearch.setQuery(query);
        }
        
        public Integer call () throws Exception {
            int count = 0;
            for (Payload p; (p = in.take()) != POISON_PAYLOAD
                     && (max < 0 || (max > 0 && out.size() < max));) {
                byte[] pfp = p.getFp();
                int i = 0;
                for (; i < fp.length; ++i) {
                    if ((pfp[i] & fp[i]) != fp[i])
                        break;
                }
                
                if (i == fp.length) {
                    Molecule mol = p.getMol();
                    mol.aromatize(); // it should already be aromatized
                    msearch.setTarget(mol);
                    int[] hits = msearch.findFirst();
                    if (hits != null) {
                        MolAtom[] atoms = mol.getAtomArray();
                        for (i = 0; i < hits.length; ++i) {
                            if (hits[i] >= 0)
                                atoms[hits[i]].setAtomMap(i+1);
                        }
                        out.put(new Result (p.doc, mol));
                        ++count;
                    }
                }
            }
            /*
            logger.info(Thread.currentThread().getName()+" "
                        +count+" subgraph isomorphisms performed!");
            */
            return count;
        }
    }

    private File baseDir;
    private Directory indexDir;
    private Directory facetDir;
    private IndexWriter indexWriter;
    private DirectoryTaxonomyWriter facetWriter;
    private Analyzer indexAnalyzer;
    private FacetsConfig facetsConfig;
    private Codebook[] codebooks;
    
    private ExecutorService threadPool;
    private boolean localThreadPool = false;
    
    public StructureIndexer (File dir) throws IOException {
        this (dir, Executors.newCachedThreadPool());
        localThreadPool = true;
    }
    
    public StructureIndexer (File dir, ExecutorService threadPool)
        throws IOException {
        if (!dir.isDirectory())
            throw new IllegalArgumentException ("Not a directory: "+dir);
        if (threadPool == null)
            throw new IllegalArgumentException ("Thread pool is null");
        
        this.threadPool = threadPool;
        this.baseDir = dir;
        
        File index = new File (dir, "index");
        if (!index.exists())
            index.mkdirs();
        indexDir = new NIOFSDirectory 
            (index, NoLockFactory.getNoLockFactory());

        File facet = new File (dir, "facet");
        if (!facet.exists())
            facet.mkdirs();
        facetDir = new NIOFSDirectory
            (facet, NoLockFactory.getNoLockFactory());

        indexAnalyzer = createIndexAnalyzer ();
        indexWriter = new IndexWriter (indexDir, new IndexWriterConfig 
                                       (LUCENE_VERSION, indexAnalyzer));
        
        File conf = new File (dir, CONFIG_FILE);
        int numDocs = indexWriter.numDocs();
        if (!conf.exists()) {
            if (numDocs > 0) {
                logger.warning("Indexer configuration "+conf+" doesn't exist "
                               +"and yet there exist "+numDocs+" documents!");
            }

            logger.info("Setting up new indexer with default values...");
            codebooks = new Codebook[CODEBOOKS];
            for (int i = 0; i < codebooks.length; ++i) {
                codebooks[i] = new Codebook (32*FPSIZE);
            }
            save (conf, codebooks);
        }
        else {
            logger.info("Loading configuration "+conf+"...");
            codebooks = load (conf);
        }
        logger.info("Codebook size: "+codebooks.length);
        logger.info("Number of documents: "+numDocs);

        facetWriter = new DirectoryTaxonomyWriter (facetDir);
        facetsConfig = new FacetsConfig ();
        facetsConfig.setMultiValued(FIELD_SOURCE, true);
        facetsConfig.setRequireDimCount(FIELD_SOURCE, true);
    }

    Analyzer createIndexAnalyzer () {
        Map<String, Analyzer> fields = new HashMap<String, Analyzer>();
        fields.put(FIELD_ID, new KeywordAnalyzer ());
        fields.put(FIELD_CODEBOOK, new KeywordAnalyzer ());
        return  new PerFieldAnalyzerWrapper 
            (new StandardAnalyzer (LUCENE_VERSION), fields);
    }
    
    public void shutdown () {
        try {
            save (new File (baseDir, CONFIG_FILE), codebooks);
            if (indexWriter != null)
                indexWriter.close();
            if (facetWriter != null)
                facetWriter.close();
            indexDir.close();
            facetDir.close();
            
            if (localThreadPool)
                threadPool.shutdown();
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    protected static Codebook[] load (File conf) throws IOException {
        Properties props = new Properties ();
        props.loadFromXML(new FileInputStream (conf));
        String param = props.getProperty("CODEBOOKS");
        if (param == null) {
            throw new IllegalArgumentException
                ("Invalid configuration; no property \"codebooks\" defined!");
        }
        List<Codebook> codebooks = new ArrayList<Codebook>();
        for (String name : param.split(",")) {
            Codebook cb = Codebook.create(name, props);
            codebooks.add(cb);
        }
        
        return codebooks.toArray(new Codebook[0]);
    }

    protected static void save (File conf, Codebook[] codebooks)
        throws IOException {
        Properties props = new Properties ();
        StringBuilder sb = new StringBuilder ();
        sb.append(codebooks[0].getName());
        for (int i = 1; i < codebooks.length; ++i)
            sb.append(","+codebooks[i].getName());
        props.setProperty("CODEBOOKS", sb.toString());
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            sb = new StringBuilder ();
            sb.append(cb.dict[0]);
            for (int j = 1; j < cb.dict.length; ++j)
                sb.append(","+cb.dict[j]);
            props.setProperty(cb.getName()+".dict", sb.toString());
            sb = new StringBuilder ();
            sb.append(cb.counts[0]);
            for (int j = 1; j < cb.counts.length; ++j)
                sb.append(","+cb.counts[j]);
            props.setProperty(cb.getName()+".counts", sb.toString());
        }
        props.storeToXML(new FileOutputStream (conf),
                         "Generated by "+StructureIndexer.class.getName()
                         +" on "+(new java.util.Date()));
    }

    public void add (String source, Molecule struc) throws IOException {
        add (source, struc.getName(), struc);
    }
    
    public void add (String source, String id, Molecule struc)
        throws IOException {
        Document doc = new Document ();
        doc.add(new StoredField (FIELD_ID, encodeDocId (source, id)));
        doc.add(new FacetField (FIELD_SOURCE, source));
        instrument (doc, struc);
        doc = facetsConfig.build(facetWriter, doc);
        indexWriter.addDocument(doc);
    }

    protected void instrument (Document doc, Molecule struc)
        throws IOException {
        MolHandler mh = new MolHandler (struc.cloneMolecule());
        mh.aromatize();
        Molecule mol = mh.getMolecule();
        mol.hydrogenize(false);
        
        byte[] fp = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            int code = cb.encode(fp);
            if (code != 0) {
                cb.incr(code); // this must be in-sync with the document count!
                doc.add(new StringField
                        (FIELD_CODEBOOK, cb.encode(code), NO));
            }
        }
        
        doc.add(new StoredField (FIELD_FINGERPRINT, fp));
        doc.add(new IntField (FIELD_POPCNT, popcnt (fp), NO));
        doc.add(new StoredField
                (FIELD_MOLFILE, mol.toFormat("csmol").getBytes("utf8")));
        doc.add(new IntField (FIELD_NATOMS, mol.getAtomCount(), NO));
        doc.add(new IntField (FIELD_NBONDS, mol.getBondCount(), NO));
        doc.add(new DoubleField (FIELD_MOLWT, mol.getMass(), NO));
    }

    public ResultEnumeration substructure (String query)
        throws Exception {
        return substructure (query, -1, 2);
    }

    public ResultEnumeration substructure
        (String query, int max, int nthreads) throws Exception {
        MolHandler mh = new MolHandler (query, true);
        return substructure (mh.getMolecule(), max, nthreads);
    }

    public ResultEnumeration substructure
        (Molecule query) throws Exception {
        return substructure (query, -1, 2);
    }
    
    public ResultEnumeration substructure
        (Molecule query, final int max, int nthreads) throws Exception {
        
        MolHandler mh = new MolHandler (query);
        mh.aromatize();
        byte[] q = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);

        IndexSearcher searcher = new IndexSearcher
            (DirectoryReader.open(indexWriter, true));

        Codebook bestCb = null;
        int bestHits = Integer.MAX_VALUE;
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            int[] eqv = cb.apply(q);
            if (eqv != null) {
                int hits = 0;
                for (int j = 0; j < eqv.length; ++j) {
                    hits += cb.count(eqv[j]);
                }
                
                if (hits < bestHits) {
                    bestHits = hits;
                    bestCb = cb;
                }
            }
        }
        
        if (bestCb == null) {
            // this means that one should iterate over ALL documents!
            logger.warning("No matches found!"); 
            return null;
        }
        
        long start = System.currentTimeMillis();
        DisjunctionMaxQuery mq = new DisjunctionMaxQuery (1.f);
        int[] eqv = bestCb.apply(q);
        for (int j = 0; j < eqv.length; ++j) {
            TermQuery tq = new TermQuery
                (new Term (FIELD_CODEBOOK, bestCb.encode(eqv[j])));
            mq.add(tq);
        }
        int total = indexWriter.numDocs();
        TopDocs hits = searcher.search(mq, total);
        logger.info("## total screened: "+hits.totalHits+"/"+bestHits
                    +" screen efficiency: "
                    +String.format("%1$.2f", 1.-(double)hits.totalHits/total)
                    +" ellapsed: "
                    +String.format("%1$.2fs",
                                   (System.currentTimeMillis()-start)*1e-3));

        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final BlockingQueue<Result> out = new LinkedBlockingQueue<Result>();
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>();
        for (int i = 0; i < nthreads; ++i)
            threads.add(threadPool.submit
                        (new GraphIso (in, out, query, q, max)));
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (doc));
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
                        /*
                        logger.info
                            ("Background threads finished; total="+total);
                        */
                        out.put(POISON_RESULT);
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    public ResultEnumeration similarity (String query, double threshold)
        throws Exception {
        MolHandler mh = new MolHandler (query, true);
        return similarity (mh.getMolecule(), threshold, -1, 2);
    }
    
    public ResultEnumeration similarity
        (Molecule query, final double threshold, // minimum tanimoto cutoff
         final int max, final int nthreads) throws Exception {

        /*
         * first calculate the minimum popcnt needed to satisfy the
         * cutoff:
         *   (a & b)/(a | b) <= threshold <= min(a,b)/max(a,b)
         * where a and b are the popcnt of the query and target, 
         * respectively. This in turn provides the following bound:
         *   a*threshold <= b <= a/threshold
         */
        MolHandler mh = new MolHandler (query);
        mh.aromatize();
        byte[] q = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);
        int popcnt = popcnt (q);

        IndexSearcher searcher = new IndexSearcher
            (DirectoryReader.open(indexWriter, true));
        int minpop = (int)(popcnt*threshold+0.5);
        Query range = NumericRangeQuery.newIntRange
            (FIELD_POPCNT, minpop, null, true, true);
        long start = System.currentTimeMillis();
        TopDocs hits = searcher.search(range, indexWriter.numDocs());
        logger.info("## range query >="+minpop
                    +": "+hits.totalHits+" ellapsed: "
                    +String.format("%1$.2fs",
                                   (System.currentTimeMillis()-start)*1e-3));
        
        final BlockingQueue<Result> out = new LinkedBlockingQueue<Result>();
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>(); 
        for (int i = 0; i < nthreads; ++i)
            threads.add(threadPool.submit
                        (new Tanimoto (in, out, q, max, threshold)));
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            in.put(new Payload (doc));
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

    public void remove (String source, String id) throws IOException {
        /* TODO: this method should only be one line but because we have 
         * to keep the document count in sync with the codebooks, we need
         * to do this..
         */
        TermQuery tq = new TermQuery
            (new Term (FIELD_ID, encodeDocId (source, id)));
        IndexSearcher searcher = new IndexSearcher
            (DirectoryReader.open(indexWriter, true));
        TopDocs hits = searcher.search(tq, indexWriter.numDocs());
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            for (String code : doc.getValues(FIELD_CODEBOOK))
                for (Codebook cb : codebooks) {
                    int c = cb.decode(code);
                    if (c >= 0)
                        cb.decr(c);
                }
        }
        indexWriter.deleteDocuments(tq);
    }

    public void stats (PrintStream ps) throws IOException {
        IndexSearcher searcher = new IndexSearcher
            (DirectoryReader.open(indexWriter, true));
        ps.println("[*** stats for "+baseDir+" ***]");
        ps.println("Popcnt histogram:");
        int[] ranges = new int[]{0, 50, 100, 150, 200, 250, 300, 350};
        for (int i = 1; i < ranges.length; ++i) {
            Query range = NumericRangeQuery.newIntRange
                (FIELD_POPCNT, ranges[i-1], ranges[i], true, false);
            TopDocs hits = searcher.search(range, indexWriter.numDocs());
            ps.println(String.format("  [%1$3d,%2$3d)", ranges[i-1], ranges[i])
                       +" "+hits.totalHits);
        }
        Query range = NumericRangeQuery.newIntRange
            (FIELD_POPCNT, ranges[ranges.length-1], null, true, false);
        TopDocs hits = searcher.search(range, indexWriter.numDocs());
        ps.println(String.format("      >%1$3d", ranges[ranges.length-1])
                       +"  "+hits.totalHits);
    }
}
