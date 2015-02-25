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

    public static final Molecule DONE = new Molecule ();
    static final Version LUCENE_VERSION = Version.LATEST;
    /**
     * Fields for each document
     */
    static final String FIELD_ID = "_id";
    static final String FIELD_SOURCE = "_source";
    static final String FIELD_CODEBOOK = "_codebook";
    static final String FIELD_FINGERPRINT = "_fingerprint";
    static final String FIELD_MOLFILE = "_molfile";
    static final String FIELD_MOLWT = "_molwt";
    static final String FIELD_NATOMS = "_natoms";
    static final String FIELD_NBONDS = "_nbonds";
    
    static final int FPSIZE = 16;
    static final int FPBITS = 2;
    static final int FPDEPTH = 6;
    
    static final int CODESIZE = 8; // 8-bit or 256
    static final int CODEBOOKS = 256;

    public static class Result {
        public final Molecule query;
        public List<Molecule> matches = new ArrayList<Molecule>();
        
        Result (Molecule query) {
            this.query = query;
        }
    }

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
        Random rand = new Random ();
        
        public Codebook (int fpsize) {
            byte[] buf = new byte[4];
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

    static class GraphIso implements Callable<Integer> {
        final BlockingQueue<Molecule> in;
        final BlockingQueue<Molecule> out;
        final MolSearch msearch;
        final Molecule poisonPill;

        GraphIso (BlockingQueue<Molecule> in,
                  BlockingQueue<Molecule> out,
                  Molecule query, Molecule poisonPill) {
            this.in = in;
            this.out = out;
            msearch = new MolSearch ();
            msearch.setQuery(query);
            this.poisonPill = poisonPill;
        }
        
        public Integer call () throws Exception {
            int count = 0;
            for (Molecule m; (m = in.take()) != poisonPill;) {
                msearch.setTarget(m);
                int[] hits = msearch.findFirst();
                if (hits != null) {
                    MolAtom[] atoms = m.getAtomArray();
                    for (int i = 0; i < hits.length; ++i) {
                        if (hits[i] >= 0)
                            atoms[hits[i]].setAtomMap(i+1);
                    }
                    out.put(m);
                    ++count;
                }
            }
            logger.info(Thread.currentThread().getName()+" "
                        +count+" subgraph isomorphisms performed!");
            
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
    
    public StructureIndexer (File dir) throws IOException {
        this (dir, Executors.newFixedThreadPool(2));
    }
    
    public StructureIndexer (File dir, ExecutorService threadPool)
        throws IOException {
        if (!dir.isDirectory())
            throw new IllegalArgumentException ("Not a directory: "+dir);
        if (threadPool == null)
            throw new IllegalArgumentException ("Thread pool is null");
        
        this.threadPool = threadPool;

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
        IndexWriterConfig conf = new IndexWriterConfig 
            (LUCENE_VERSION, indexAnalyzer);
        indexWriter = new IndexWriter (indexDir, conf);

        facetWriter = new DirectoryTaxonomyWriter (facetDir);
        facetsConfig = new FacetsConfig ();
        facetsConfig.setMultiValued(FIELD_SOURCE, true);
        facetsConfig.setRequireDimCount(FIELD_SOURCE, true);
            
        codebooks = new Codebook[CODEBOOKS];
        for (int i = 0; i < codebooks.length; ++i) {
            codebooks[i] = new Codebook (32*FPSIZE);
        }
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
            if (indexWriter != null)
                indexWriter.close();
            if (facetWriter != null)
                facetWriter.close();
            indexDir.close();
            facetDir.close();       
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public void add (String source, Molecule struc) throws IOException {
        add (source, struc.getName(), struc);
    }
    
    public void add (String source, String id, Molecule struc)
        throws IOException {
        Document doc = new Document ();
        doc.add(new StoredField (FIELD_ID, source+"::"+id));
        if (source != null) {
            doc.add(new FacetField (FIELD_SOURCE, source));
        }
        instrument (doc, struc);
        doc = facetsConfig.build(facetWriter, doc);
        indexWriter.addDocument(doc);
    }

    protected void instrument (Document doc, Molecule struc)
        throws IOException {
        MolHandler mh = new MolHandler (struc.cloneMolecule());
        mh.aromatize();
        Molecule mol = mh.getMolecule();
        byte[] fp = mh.generateFingerprintInBytes(FPSIZE, FPBITS, FPDEPTH);
        for (int i = 0; i < codebooks.length; ++i) {
            Codebook cb = codebooks[i];
            int code = cb.encode(fp);
            if (code != 0) {
                cb.incr(code); // this must be insync with the document count!
                doc.add(new StringField
                        (FIELD_CODEBOOK, cb.encode(code), NO));
            }
        }
        
        doc.add(new StoredField (FIELD_FINGERPRINT, fp));
        doc.add(new StoredField
                (FIELD_MOLFILE, mol.toFormat("csmol").getBytes("utf8")));
        doc.add(new IntField (FIELD_NATOMS, mol.getAtomCount(), NO));
        doc.add(new IntField (FIELD_NBONDS, mol.getBondCount(), NO));
        doc.add(new DoubleField (FIELD_MOLWT, mol.getMass(), NO));
    }

    public BlockingQueue<Molecule> substructure (String query)
        throws Exception {
        return substructure (query, 2);
    }
    
    public BlockingQueue<Molecule> substructure (String query, int nthreads)
        throws Exception {
        MolHandler mh = new MolHandler (query, true);
        return substructure (mh.getMolecule(), nthreads);
    }

    public BlockingQueue<Molecule> substructure
        (Molecule query) throws Exception {
        return substructure (query, 2);
    }
    
    public BlockingQueue<Molecule> substructure
        (Molecule query, int nthreads) throws Exception {
        
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

        Molecule poisonPill = new Molecule ();
        BlockingQueue<Molecule> in = new LinkedBlockingQueue<Molecule>();
        BlockingQueue<Molecule> out = new LinkedBlockingQueue<Molecule>();
        for (int i = 0; i < nthreads; ++i)
            threadPool.submit(new GraphIso (in, out, query, poisonPill));
        
        for (int i = 0; i < hits.totalHits; ++i) {
            Document doc = searcher.doc(hits.scoreDocs[i].doc);
            BytesRef ref = doc.getBinaryValue(FIELD_FINGERPRINT);
            if (ref.length != q.length) {
                logger.log(Level.SEVERE,
                           "Bogus fingerprint stored in Document "
                           +doc.get(FIELD_ID));
                continue;
            }
            
            int j = 0;
            for (; j < q.length; ++j) {
                if ((ref.bytes[j] & q[j]) != q[j])
                    break;
            }
            
            if (j == q.length) {
                // now do (sub-) graph isomorphism here!
                BytesRef molref = doc.getBinaryValue(FIELD_MOLFILE);
                mh.setMolecule(new String (molref.bytes,
                                           molref.offset, molref.length));
                mh.aromatize();
                in.put(mh.getMolecule());
            }
        }
        for (int i = 0; i < nthreads; ++i)
            in.put(poisonPill);
        out.put(DONE);
        
        return out;
    }

    public void remove (String source, String id) throws IOException {
        TermQuery tq = new TermQuery (new Term (FIELD_ID, source+"::"+id));
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

    public static void main (String[] argv) throws Exception {
        if (argv.length < 2) {
            System.err.println("Usage: StructureIndexer INDEXDIR FILES...");
            System.exit(1);
        }

        File dir = new File (argv[0]);
        if (!dir.exists())
            dir.mkdirs();
        ExecutorService threadPool = Executors.newCachedThreadPool();
        StructureIndexer indexer = new StructureIndexer (dir, threadPool);
        logger.info("Indexing structures...");
        long start = System.currentTimeMillis(), total = 0;
        for (int i = 1; i < argv.length; ++i) {
            File file = new File (argv[i]);
            MolImporter mi = new MolImporter (argv[i]);
            int count = 0;
            for (Molecule m = new Molecule (); mi.read(m); ++count) {
                indexer.add(file.getName(), m);
            }
            logger.info(file.getName()+": "+count);
            total += count;
        }
        logger.info("Indexing time "+String.format
                    ("%1$.2fs",1e-3*(System.currentTimeMillis()-start))
                    +" for "+total+" structures!");

        start = System.currentTimeMillis();
        BlockingQueue<Molecule> q = indexer.substructure("c1ccncc1", 3);
        if (q != null) {
            double ellapsed = (System.currentTimeMillis()-start)*1e-3;
            int count = 0;
            for (Molecule m; (m = q.take()) != DONE; ++count) {
                //System.out.println(m.toFormat("smiles:q"));
            }
            logger.info(count+" matches found in "
                        +String.format("%1$.2fs", ellapsed));
        }
        threadPool.shutdown();
    }
}
