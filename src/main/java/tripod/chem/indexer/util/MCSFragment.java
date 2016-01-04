
package tripod.chem.indexer.util;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.sql.*;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

public class MCSFragment {
    static class FragmentGroup {
	Map<String, String> group;
	long fragid;

	FragmentGroup (long fragid, Map<String, String> g) {
	    this.group = new HashMap<String, String>(g);
	    this.fragid = fragid;
	}
	FragmentGroup () {}
    }

    static final FragmentGroup EMPTY_GROUP = new FragmentGroup ();

    static class Worker extends Thread {
	BlockingQueue<FragmentGroup> queue;

	public Worker (BlockingQueue<FragmentGroup> queue) {
	    this.queue = queue;
	}

	public void run () {
	    try {
		while (true) {
		    FragmentGroup fg = queue.take();
		    if (fg == EMPTY_GROUP) {
			break;
		    }
		    try {
			doMCS (fg.fragid, fg.group);
		    }
		    catch (Exception ex) {
			ex.printStackTrace();
		    }
		}
	    }
	    catch (InterruptedException ex) {
	    }
	}
    }

    static void doMCS (long fragid, Map<String, String> group) 
	throws Exception {
	MolHandler mh = new MolHandler ();
	MolStandardizer molstd = new MolStandardizer ();

	long id = Thread.currentThread().getId();

	System.out.println
	    ("# Thread " + id + " processing group " + fragid + "...");
	Vector<Molecule> mols = new Vector<Molecule>();

	for (Map.Entry<String, String> e : group.entrySet()) {
	    //System.out.println(e.getKey());
	    mh.setMolecule(e.getValue());
	    Molecule m = mh.getMolecule();
	    if (m.getAtomCount() < 100) {
		// don't process peptide...
		molstd.standardize(m);
		m.calcHybridization();
		m.setName(e.getKey());
		mols.add(m);
	    }
	}

	System.out.println
	    ("# Thread " + id + " " + mols.size() + " structure(s)");

	PrintStream ps = new PrintStream 
	    (new FileOutputStream ("frag_"+fragid+"_mcs.smi"), true);

	Molecule []mary = mols.toArray(new Molecule[0]);
	AlmostMCS mcs = new AlmostMCS ();

	int count = 0;
	for (int i = 0; i < mary.length; ++i) {
	    Molecule mi = mary[i];
	    mcs.setQuery(mi);
	    for (int j = i + 1; j < mary.length; ++j) {
		Molecule mj = mary[j];
		mcs.setTarget(mj);
		if (mcs.search()) {
		    Molecule core = mcs.getResultAsMolecule(false);
		    int c = ChemUtil.complexity(core);
		    if (c > 110) {
			String smi = core.toFormat("smiles:ua_bas0");
			ps.println(smi + "\t"  + mary[i].getName() 
				   + ":" + mary[j].getName() + "\t" + c);
			if (++count % 5000 == 0) {
			    System.out.print(".");
			}
		    }
		}
	    }
	}

	ps.close();
	System.out.println
	    ("# Thread " + id + ": " + count + " MCS(s) generated!");
	System.out.println("# Thread " + id + " is done; current active " 
			   + Thread.activeCount());
    }

    public static void main (String []argv) throws Exception {
	if (argv.length == 0) {
	    System.out.println("Usage: MCSFragment TYPE");
	    System.exit(1);
	}

	int type = Integer.parseInt(argv[0]);

	Class.forName("oracle.jdbc.driver.OracleDriver");
	Connection con = DriverManager.getConnection
	    ("jdbc:oracle:thin:spotlite/spotlite@ncgcprobe.nhgri.nih.gov:1521:probedb");

	String sql = "SELECT c.smiles_iso,"
	    +" c.sample_id, "
	    +" a.frag_id "
	    +"FROM fragment a,"
	    +"  fragment_group b,"
	    +"  registry.ncgc_sample c "
	    +"WHERE a.type_id = ? "
	    +"AND a.flags is null "
	    +"AND a.frag_id = b.frag_id "
	    +"AND b.sample_id = c.sample_id "
	    +"AND c.smiles_can is not null "
	    +"order by a.frag_id"
	    ;

	PreparedStatement pstm = con.prepareStatement(sql);
	pstm.setInt(1, type);

	long pid = 0;
	Map<String, String> group = new HashMap<String, String>();

	// use only four threads...
	int maxthreads = Integer.getInteger("thread.count", 4);
	
	BlockingQueue<FragmentGroup> queue = 
	    new LinkedBlockingQueue<FragmentGroup> ();

	Worker workers[] = new Worker[maxthreads];
	for (int i = 0; i < workers.length; ++i) {
	    workers[i] = new Worker (queue);
	    workers[i].start();
	}

	try {
	    ResultSet rset = pstm.executeQuery();
	    while (rset.next()) {
		String smi = rset.getString(1);
		String sample = rset.getString(2);
		long id = rset.getLong(3);
		if (id != pid) {
		    if (group.size() > 1) {
			queue.put(new FragmentGroup (pid, group));
		    }
		    group.clear();
		}
		group.put(sample, smi);
		pid = id;
	    }
	    
	    if (group.size() > 1) {
		queue.put(new FragmentGroup (pid, group));
	    }
	    
	    queue.put(EMPTY_GROUP);
	    rset.close();
	}
	catch (InterruptedException ex) {
	}

	pstm.close();
	con.close();
    }
}
