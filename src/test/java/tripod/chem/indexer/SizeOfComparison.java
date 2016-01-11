package tripod.chem.indexer;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import chemaxon.formats.MolFormatException;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import net.sourceforge.sizeof.SizeOf;

public class SizeOfComparison {

	public static void main(String[] args) throws InvalidSmilesException, MolFormatException {
		
		String smilesString = "CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5";
		
		System.out.println("for smiles string " + smilesString);
		SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		
		List<Object> list = new ArrayList<>();
		
		for(int i= 0; i< 1000; i++){
			list.add(smilesParser.parseSmiles(smilesString));
		}
		System.out.println("1000 CDK containers size = " + SizeOf.humanReadable(SizeOf.deepSizeOf(list)));

		List<Object> list2 = new ArrayList<>();
		
		for(int i= 0; i< 1000; i++){
			Molecule molecule = new MolHandler(smilesString).getMolecule();
		//	molecule.clean(0, "");
		//	molecule.clean(1, "");
			molecule.clean(0, null);
		//	molecule.clean(3, "");
			list2.add(molecule);
		}
		
		
		System.out.println("1000 JChem mol size = " + SizeOf.humanReadable(SizeOf.deepSizeOf(list2)));

	}

	
}
