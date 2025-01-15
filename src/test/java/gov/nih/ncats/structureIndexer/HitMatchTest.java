package gov.nih.ncats.structureIndexer;

import gov.nih.ncats.molwitch.Chemical;
import org.apache.lucene.document.Document;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

public class HitMatchTest extends AbstractStructureIndexerTest {

    @Test
    public void TestInstrument() throws IOException, NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        String localDir = System.getProperty("user.dir");
        System.out.printf("localDir = %s\n", localDir);
        String molfilePath = localDir + "/src/test/resources/mols/85740c44-7680-4025-bd22-98d2f764469e.mol";
        File molFile = new File(molfilePath);
        Chemical chem = Chemical.parseMol(molFile);
        Integer chemicalId =1;
        Document doc = new Document();
        StructureIndexer.Payload payload = new StructureIndexer.Payload(chemicalId, doc);

        StructureIndexer.Result result =new StructureIndexer.Result(payload);

        String tmpdir = System.getProperty("java.io.tmpdir") + "/index";
        File tempDir = new File(tmpdir);
        tempDir.mkdir();

        Method method = indexer.getClass().getDeclaredMethod("instrument", Document.class, Chemical.class);
        method.setAccessible(true);
        method.invoke(indexer, doc, chem);

        Assert.assertNotNull(result);
    }
}
