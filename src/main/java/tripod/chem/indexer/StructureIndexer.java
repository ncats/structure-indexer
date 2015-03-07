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

import org.apache.lucene.index.AtomicReaderContext;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.SerialMergeScheduler;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.Terms;
import org.apache.lucene.index.TermsEnum;
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
import org.apache.lucene.search.FilteredQuery;

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
    public static final String FIELD_TEXT = "text";
    public static final String FIELD_ID = "_id";
    public static final String FIELD_SOURCE = "_source";
    public static final String FIELD_CODEBOOK = "_codebook";
    public static final String FIELD_FINGERPRINT = "_fingerprint";
    public static final String FIELD_FIELDS = "_fields";
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
        Molecule mol;
        byte[] fp;

        Payload (int id, Document doc) {
            this.id = id;
            this.doc = doc;
        }
        
        Payload () {
            id = -1;
            doc = null;
        }

        public byte[] getFp () {
            if (fp == null) {
                BytesRef ref = doc.getBinaryValue(FIELD_FINGERPRINT);
                if (ref.offset > 0) {
                    fp = new byte[ref.length];
                    System.arraycopy(ref.bytes, ref.offset, fp, 0, ref.length);
                }
                else {
                    fp = ref.bytes;
                }
            }
            return fp;
        }
        public int getId () { return id; }
        public Document getDoc () { return doc; }
        public Molecule getMol () {
            if (mol == null) {
                BytesRef molref = doc.getBinaryValue(FIELD_MOLFILE);
                try {
                    MolHandler mh = new MolHandler
                        (new String (molref.bytes,
                                     molref.offset, molref.length));
                    mol = mh.getMolecule();
                    mol.setName(doc.get(FIELD_ID));
                    for (IndexableField f : doc.getFields(FIELD_FIELDS)) {
                        String v = doc.get(f.stringValue());
                        mol.setProperty(f.stringValue(), v);
                    }
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

    static public class Result implements Comparable<Result> {
        final int id;
        final Document doc;
        final Molecule mol;
        final Double similarity;
        
        Result (Payload payload, Double similarity) {
            id = payload.getId();
            doc = payload.getDoc();
            mol = payload.getMol();
            this.similarity = similarity;
        }
        
        Result (Payload payload) {
            this (payload, null);
        }
        
        Result () {
            id = Integer.MAX_VALUE;
            doc = null;
            mol = null;
            similarity = null;
        }

        public int compareTo (Result r) {
            int d = 0;
            if (similarity != null && r.similarity != null) {
                double dif = r.similarity - similarity;
                if (dif > 0.) d = 1;
                else if (dif < 0.) d = -1;
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
        public Molecule getMol () { return mol; }
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
        for (int i = 0; i < b.length; ++i) {
            c += Integer.bitCount(b[i] & 0xff);
        }
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
            int count = 0, total = 0;
            for (Payload p; (p = in.take()) != POISON_PAYLOAD
                     && (max <= 0 || (max > 0 && out.size() < max)); ++total) {
                int a = 0, b = 0;
                byte[] fp = p.getFp();
                for (int j = 0; j < query.length; ++j) {
                    a += Integer.bitCount(query[j] & fp[j]);
                    b += Integer.bitCount(query[j] | fp[j]);
                }
            
                double tan = (double)a/b;
                if (tan >= threshold) {
                    ++count;
                    out.put(new Result (p, tan));
                }
            }
            /*
            logger.info(Thread.currentThread().getName()+" "
                        +count+"/"+total
                        +" passed tanimoto cutoff "+threshold+"!");
            */
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
                     && (max <= 0 || (max > 0 && out.size() < max));) {
                byte[] pfp = p.getFp();
                assert pfp.length == fp.length;
                
                int i = 0, a = 0, b = 0;
                for (; i < fp.length; ++i) {
                    if ((pfp[i] & fp[i]) != fp[i]) {
                        break;
                    }
                    a += Integer.bitCount(pfp[i] & fp[i]);
                    b += Integer.bitCount(pfp[i] | fp[i]);
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
                        out.put(new Result (p, (double)a/b));
                        ++count;
                    }
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
                out.put(new Result (p, null));
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

        indexDir = new NIOFSDirectory(index, NoLockFactory.getNoLockFactory());
        metaDir = new NIOFSDirectory (meta, NoLockFactory.getNoLockFactory());
        facetDir = new NIOFSDirectory
            (facet, NoLockFactory.getNoLockFactory());

        if (!readOnly) {
            metaWriter = new IndexWriter
                (metaDir, new IndexWriterConfig
                 (LUCENE_VERSION, indexAnalyzer));
            indexWriter = new IndexWriter (indexDir, new IndexWriterConfig 
                                           (LUCENE_VERSION, indexAnalyzer));

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
        }
        else {
            codebooks = load (DirectoryReader.open(metaDir));
            indexReader = DirectoryReader.open(indexDir);           
        }
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
            (new StandardAnalyzer (LUCENE_VERSION), fields);
    }

    protected synchronized DirectoryReader getReader () throws IOException {
        DirectoryReader reader = DirectoryReader.openIfChanged(indexReader);
        if (reader != null) {
            indexReader.close();
            indexReader = reader;
        }
        return indexReader;
    }

    protected IndexSearcher getIndexSearcher () throws IOException {
        return new IndexSearcher (getReader (), threadPool);
    }

    public File getBasePath () { return baseDir; }
    public String[] getFields () {
        List<String> fields = new ArrayList<String>();
        try {
            for (AtomicReaderContext ctx : getReader().leaves()) {
                Terms terms = ctx.reader().terms(FIELD_FIELDS);
                //logger.info("terms "+terms);
                if (terms != null) {
                    TermsEnum en = terms.iterator(null);
                    for (BytesRef ref; (ref = en.next()) != null; ) {
                        fields.add(new String (ref.bytes, ref.offset, ref.length));
                    }
                }
            }
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        return fields.toArray(new String[0]);
    }
    
    public void shutdown () {
        try {
            if (indexReader != null)
                indexReader.close();
            if (indexWriter != null)
                indexWriter.close();
            if (facetWriter != null)
                facetWriter.close();
            
            if (metaWriter != null) {
                for (Codebook cb : codebooks)
                    update (cb);
                metaWriter.close();
            }
            
            indexDir.close();
            facetDir.close();
            metaDir.close();
            if (localThreadPool)
                threadPool.shutdown();
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    protected void update (Codebook cb) throws IOException {
        Term term = new Term (FIELD_ID, cb.getName());
        //logger.info("Persist codebook "+cb);    
        metaWriter.updateDocument(term, Codebook.instrument(cb));
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
        try {
            return getReader().numDocs();
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        return -1;
    }

    public Map<String, Integer> getSources () throws IOException {
        Map<String, Integer> map = new TreeMap<String, Integer>();
        try {
            IndexSearcher searcher = getIndexSearcher ();
            FacetsCollector fc = new FacetsCollector ();
            TaxonomyReader taxon = new DirectoryTaxonomyReader (facetWriter);
            TopDocs hits = FacetsCollector.search
                (searcher, new MatchAllDocsQuery (),
                 null, searcher.getIndexReader().numDocs(), fc);
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
            MolHandler mh = new MolHandler (struc);
            add (source, id, mh.getMolecule());
        }
        catch (Exception ex) {
            throw new IllegalArgumentException ("Bogus molecule format", ex);
        }
    }
    
    public void add (String id, Molecule struc) throws IOException {
        add (null, id, struc);
    }
    
    public void add (String source, String id, Molecule struc)
        throws IOException {
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");
        
        Document doc = new Document ();
        doc.add(new StringField (FIELD_ID, id, YES));
        if (source != null) {
            doc.add(new FacetField (FIELD_SOURCE, source));
            doc.add(new StringField (FIELD_SOURCE, source, NO));
            doc.add(new TextField (FIELD_TEXT, source, NO));
        }
        doc.add(new TextField (FIELD_TEXT, id, NO));
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
        
        for (int i = 0; i < struc.getPropertyCount(); ++i) {
            String prop = struc.getPropertyKey(i);
            String value = struc.getProperty(prop);
            if (value != null) {
                doc.add(new TextField (FIELD_TEXT, value, NO));
                doc.add(new TextField (FIELD_TEXT, prop, NO));
                doc.add(new StoredField (FIELD_FIELDS, prop));
                try {
                    double dv = Double.parseDouble(value);
                    doc.add(new DoubleField (prop, dv, NO));
                }
                catch (NumberFormatException ex) {
                }
                try {
                    int iv = Integer.parseInt(value);
                    doc.add(new IntField (prop, iv, NO));
                }
                catch (NumberFormatException ex) {
                }
                doc.add(new StoredField (prop, value));
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
    
    public void remove (String source, String id) throws IOException {
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");
        
        /* TODO: this method should only be one line but because we have 
         * to keep the document count in sync with the codebooks, we need
         * to do this..
         */
        TermQuery tq = new TermQuery (new Term (FIELD_ID, id));
        IndexSearcher searcher = getIndexSearcher ();
        int total = searcher.getIndexReader().numDocs();
        Filter filter = null;
        if (source != null) {
            filter = new FieldCacheTermsFilter (FIELD_SOURCE, source);
        }
        TopDocs hits = searcher.search(tq, filter, total);
        
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
        indexWriter.deleteDocuments(tq);
    }

    public void remove (String source) throws IOException {
        // this potentially can dramatically alter the counts, so we need
        if (indexWriter == null)
            throw new RuntimeException ("Index is read-only!");

        indexWriter.deleteDocuments(new Term (FIELD_SOURCE, source));
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
        MolHandler mh = new MolHandler (query, true);
        return substructure (mh.getMolecule(), max, nthreads, filters);
    }

    public ResultEnumeration substructure
        (Molecule query) throws Exception {
        return substructure (query, -1, 2, null);
    }
    
    public ResultEnumeration substructure
        (Molecule query, Filter... filters) throws Exception {
        return substructure (query, -1, 2, filters);
    }

    public ResultEnumeration substructure
        (Molecule query, final int max, int nthreads) throws Exception {
        return substructure (getIndexSearcher (), query, max, nthreads, null);
    }
    
    public ResultEnumeration substructure
        (Molecule query, final int max, int nthreads, Filter... filters)
        throws Exception {
        return substructure (getIndexSearcher (),
                             query, max, nthreads, filters);
    }
    
    protected ResultEnumeration substructure
        (IndexSearcher searcher, Molecule query,
         final int max, int nthreads, Filter... filters) throws Exception {
        MolHandler mh = new MolHandler (query);
        mh.aromatize();
        byte[] qfp = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);
        
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
                
                if (hits < bestHits) {
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
        final BlockingQueue<Result> out = new PriorityBlockingQueue<Result>();
        final List<Future<Integer>> threads = new ArrayList<Future<Integer>>();
        for (int i = 0; i < nthreads; ++i)
            threads.add(threadPool.submit
                        (new GraphIso (in, out, query, qfp, max)));
        
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

    public ResultEnumeration similarity
        (String query, double threshold, int max) throws Exception {
        MolHandler mh = new MolHandler (query, true);   
        return similarity (mh.getMolecule(), threshold, max, 2);
    }
    
    public ResultEnumeration similarity (String query, double threshold)
        throws Exception {
        MolHandler mh = new MolHandler (query, true);
        return similarity (mh.getMolecule(), threshold, -1, 2);
    }

    public ResultEnumeration similarity
        (String query, double threshold, Filter... filters) throws Exception {
        MolHandler mh = new MolHandler (query, true);
        return similarity (mh.getMolecule(), threshold, filters);
    }

    public ResultEnumeration similarity
        (Molecule query, double threshold, Filter... filters) throws Exception {
        return similarity (query, threshold, -1, 2, filters);
    }
    
    public ResultEnumeration similarity
        (Molecule query, final double threshold, 
         final int max, final int nthreads) throws Exception {
        return similarity (query, threshold, max, nthreads, null);
    }
    
    public ResultEnumeration similarity
        (Molecule query, final double threshold, 
         final int max, final int nthreads, Filter... filters)
        throws Exception {
        return similarity
            (getIndexSearcher (), query, threshold, max, nthreads, filters);
    }

    protected ResultEnumeration similarity
        (IndexSearcher searcher, Molecule query, final double threshold,
         final int max, final int nthreads, Filter... filters)
        throws Exception {
        /*
         * first calculate the minimum popcnt needed to satisfy the
         * cutoff:
         *   (a & b)/(a | b) <= threshold <= min(a,b)/max(a,b)
         * where a and b are the popcnt of the query and target, 
         * respectively. This in turn provides the following bound:
         *   a*threshold <= b
         */
        MolHandler mh = new MolHandler (query);
        mh.aromatize();
        byte[] q = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);
        int popcnt = popcnt (q);

        int minpop = (int)(popcnt*threshold+0.5);
        Query range = NumericRangeQuery.newIntRange
            (FIELD_POPCNT, minpop, null, true, true);
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
        
        final BlockingQueue<Result> out = new PriorityBlockingQueue<Result>();
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
                        out.put(POISON_RESULT);
                    }
                    catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        
        return new ResultEnumeration (out);
    }

    public ResultEnumeration search (Query query) throws Exception {
        return search (query, null, 0);
    }

    public ResultEnumeration search (Filter filter) throws Exception {
        return search (new MatchAllDocsQuery (), filter, 0);
    }
    
    public ResultEnumeration search (Query query, Filter filter, int max)
        throws Exception {
        if (query == null && filter == null)
            throw new IllegalArgumentException
                ("Both query and filter are null!");
        
        IndexSearcher searcher = getIndexSearcher ();   
        final BlockingQueue<Result> out = new LinkedBlockingQueue<Result>();
        final BlockingQueue<Payload> in = new LinkedBlockingQueue<Payload>();
        final Future<Integer> future =
            threadPool.submit(new Output(in, out, max));
        
        if (max <= 0)
            max = searcher.getIndexReader().numDocs();
        TopDocs hits = searcher.search(query, filter, max);
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
        QueryParser parser = new QueryParser ("text", indexAnalyzer);
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
        int[] hist = new int[range.length+2];
        
        int numDocs = searcher.getIndexReader().numDocs();
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
        
        return hist;
    }

    protected int[] histogram (IndexSearcher searcher,
                               String field, double[] range)
        throws IOException {
        int[] hist = new int[range.length+2];
        
        int numDocs = searcher.getIndexReader().numDocs();
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
        
        return hist;
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
