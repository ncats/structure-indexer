package gov.nih.ncats.structureIndexer;


import java.util.BitSet;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.fingerprint.Fingerprint;
import gov.nih.ncats.molwitch.fingerprint.Fingerprinter;


public class ECFingerprinter implements Fingerprinter{
    private int nBits;
    private int maxLength=8;
    private int bitsPerString=1;
    private boolean encodeWhole = true;

    public ECFingerprinter(int nBits, int maxLength, int bitsPerString, boolean encodeWhole){
        this.nBits=nBits;
        this.maxLength=maxLength;
        this.bitsPerString=bitsPerString;
        this.encodeWhole=encodeWhole;
    }

    @Override
    public Fingerprint computeFingerprint(Chemical chemical) {
        ECFingerprint ecf = new ECFingerprint(nBits, maxLength, bitsPerString, encodeWhole);
        long[] fp1=ecf.getFingerprint(chemical);
        return new Fingerprint(BitSet.valueOf(fp1), nBits);
    }
}
