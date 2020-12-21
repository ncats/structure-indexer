package gov.nih.ncats.structureIndexer;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;

public class ConcatFingerprinter implements Fingerprinter{
	private List<Tuple<Fingerprinter,Integer>> fps=new ArrayList<>();
	private int totLength=0;
	
	public ConcatFingerprinter(){}
	
	public ConcatFingerprinter addFP(Fingerprinter fp, int length){
		int offset=totLength;
		fps.add(Tuple.of(fp,offset));
		totLength=offset+length;
		return this;
	}
	
	@Override
	public Fingerprint computeFingerprint(Chemical chemical) {
		
		BitSet bs1= new BitSet();
		
		
		fps.stream()
		   .map(Tuple.kmap(fp->fp.computeFingerprint(chemical)))
		   .forEach(t->{
			   int off=t.v();
			   Fingerprint fp1 = t.k();
			   fp1.toBitSet().stream()
				   .forEach(i->{
					   bs1.set(i+off);
				   });			   
		   });
		return new Fingerprint(bs1);
	}
	
	
	public Fingerprinter folded(int nlength){
		Fingerprinter _this=this;
		
		return new Fingerprinter(){

			@Override
			public Fingerprint computeFingerprint(Chemical chemical) {
				Fingerprint fp1 = _this.computeFingerprint(chemical);
				BitSet bs1= new BitSet();
				
				fp1.toBitSet().stream()
				   .forEach(i->{
					   bs1.set((i*101)%nlength);
				   });
				
				return new Fingerprint(bs1);
			}
			
		};
	}
	
	
}
