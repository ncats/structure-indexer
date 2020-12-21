package gov.nih.ncats.structureIndexer;


import java.util.BitSet;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;


public class ECFingerprinter implements Fingerprinter{
	int nBits;
	int MAX_LENGTH=8;
	int BITS_PER_STRING=1;
	
	
	public ECFingerprinter(int nBits, int MAX_LENGTH, int BITS_PER_STRING){
		this.nBits=nBits;
		this.MAX_LENGTH=MAX_LENGTH;
		this.BITS_PER_STRING=BITS_PER_STRING;
	}
	
	@Override
	public Fingerprint computeFingerprint(Chemical chemical) {
		ECFingerprint ecf = new ECFingerprint(nBits, MAX_LENGTH, BITS_PER_STRING);
		
		long[] fp1=ecf.getFingerprint(chemical);
		return new Fingerprint(BitSet.valueOf(fp1));
	}
}
