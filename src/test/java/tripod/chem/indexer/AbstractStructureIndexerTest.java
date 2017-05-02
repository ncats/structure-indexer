package tripod.chem.indexer;

import java.io.IOException;

import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

public class AbstractStructureIndexerTest {

	 	@Rule
	    public TemporaryFolder tmpDir = new TemporaryFolder();

	    protected StructureIndexer indexer;

	   
	    
	    @Before
	    public void createIndexer() throws IOException {
	        indexer = StructureIndexer.open(tmpDir.getRoot());
	    }

	    @After
	    public void shutdownIndexer() {
	        if (indexer != null) {
	            indexer.shutdown();
	        }
	    }
}
