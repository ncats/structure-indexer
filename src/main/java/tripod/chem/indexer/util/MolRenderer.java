
package tripod.chem.indexer.util;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import chemaxon.marvin.MolPrinter;
import chemaxon.marvin.common.UserSettings;
import chemaxon.marvin.paint.DispOptConsts;
import chemaxon.struc.Molecule;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MDocument;
import chemaxon.util.MolHandler;
import chemaxon.struc.StereoConstants;

import java.awt.Image;
import java.awt.AlphaComposite;
import java.awt.Rectangle;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Color;
import java.awt.Stroke;
import java.awt.BasicStroke;
import java.awt.color.ColorSpace;
import java.awt.geom.RoundRectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.BufferedImageOp;
import java.awt.image.ColorConvertOp;
import java.awt.image.ConvolveOp;
import java.awt.image.Kernel;

import java.util.logging.Level;
import java.util.logging.Logger;

public class MolRenderer implements StereoConstants {
    static final private Logger logger = 
	Logger.getLogger(MolRenderer.class.getName());

    static final ConcurrentMap<Integer, BufferedImageOp> OperCache = 
        new ConcurrentHashMap<Integer, BufferedImageOp>();
    

    static final ColorSpace GREYCS = ColorSpace.getInstance(ColorSpace.CS_GRAY);
	static final ColorConvertOp GREYOP = new ColorConvertOp(GREYCS, null);  

    static final Color DEFAULT_BACKGROUND = Color.white;

    private MolPrinter renderer;
    private Color background;
    private Color borderColor = new Color (0, 0, 0, 220);
    private boolean border = false;
    private boolean shadow; // show shadow
    private boolean monochrome = false;
    private int shadowOffset = 3;
    private float shadowRadius = 0.01f; // portional to the size
    private float shadowTranslucent = .25f;

    public MolRenderer () {
        this (true, DEFAULT_BACKGROUND);
    }
    public MolRenderer (boolean shadow) {
        this (shadow, DEFAULT_BACKGROUND);
    }
    public MolRenderer (Color background) {
        this (true, background);
    }
    public MolRenderer (boolean shadow, Color background) {
        this.renderer = createMolPrinter ();
        this.shadow = shadow;
        setBackground (background);
    }

    static public MolPrinter createMolPrinter () {
        MolPrinter renderer = new MolPrinter();
        renderer.setAtomsize(.6);
	renderer.setBondSpacing(.3);
        renderer.setDispopts
            (UserSettings.IMPLICITH_HETERO, UserSettings.IMPLICITH_MASK);
        renderer.setDispopts
            (renderer.getDispopts() & ~DispOptConsts.ATMAP_FLAG);
        return renderer;
    }

    public void setBackground (Color color) {
        if (color == null) {
            color = DEFAULT_BACKGROUND;
        }
	renderer.setBackgroundColor(color);
	this.background = color;
    }

    public void setMonochrome (boolean monochrome) {
        this.monochrome = monochrome;
    }
    public boolean isMonochrome () { return monochrome; }

    public void setBorderVisible (boolean border) { this.border = border; }
    public boolean getBorderVisible () { return border; }
    public void setBorderColor (Color color) { this.borderColor = color; }
    public Color getBorderColor () { return borderColor; }

    public void setShadowVisible (boolean shadow) {
        this.shadow = shadow;
    }
    public boolean getShadowVisible () { return shadow; }

    public void setShadowParams (int offset, float radius, float translucent){
        shadowOffset = offset;
        shadowRadius = radius;
        shadowTranslucent = translucent;
    }
    public int getShadowOffset () { return shadowOffset; }
    public void setShadowOffset (int offset) { this.shadowOffset = offset; }
    public float getShadowRadius () { return shadowRadius; }
    public void setShadowRadius (float radius) { shadowRadius = radius; }
    public float getShadowTranslucent () { return shadowTranslucent; }
    public void setShadowTranslucent (float translucent) {
        shadowTranslucent = translucent;
    }

    public MolPrinter getMolPrinter () { return renderer; }
    public void setMolPrinter (MolPrinter renderer) {
        this.renderer = renderer;
    }

    public void setAtomMapVisible (boolean visible) {
	if (visible) {
	    renderer.setDispopts
                (renderer.getDispopts() | DispOptConsts.ATMAP_FLAG);
	}
	else {
	    renderer.setDispopts
                (renderer.getDispopts() & ~DispOptConsts.ATMAP_FLAG);
	}
    }
    public boolean getAtomMapVisible () {
	return (renderer.getDispopts() & DispOptConsts.ATMAP_FLAG) 
	    == DispOptConsts.ATMAP_FLAG;
    }

    public void setDisplayOptions (int flags) {
	renderer.setDispopts(flags);
    }
    public int getDisplayOptions () { return renderer.getDispopts(); }

    public void setStereoCentersVisible (boolean visible) {
	renderer.setChiralitySupport(visible? CHIRALITYSUPPORT_ALL
                                     : CHIRALITYSUPPORT_NONE);
    }
    public boolean getStereoCentersVisible () {
	return renderer.getChiralitySupport() != CHIRALITYSUPPORT_ALL;
    }

    public void setEZVisible (boolean visible) {
	renderer.setEzVisible(visible);
    }
    public boolean getEZVisible () { return renderer.isEzVisible(); }


    public void renderMol (Graphics2D g2, Molecule m, 
			   int width, int height, boolean round) {
	renderMol (g2, m, 0, 0, width, height, round);
    }

    public void renderMol (Graphics2D g2, Molecule m, int x, int y,
			   int width, int height, boolean round) {
	Molecule mol = m.cloneMolecule();
	mol.dearomatize();
        if (!getAtomMapVisible ()) {
            // only highlight the atoms if atom map isn't visible
            for (MolAtom a : mol.getAtomArray()) {
                int map = a.getAtomMap();
                if (map > 0) {
                    a.setAtomMap(0);
                    a.setSetSeq(3); // blue
                }
            }
        }

        renderer.setColorScheme
            (monochrome || mol.hasAtomSet() || mol.hasBondSet()
             ? DispOptConsts.MONO_SCHEME_S : DispOptConsts.CPK_SCHEME_S);

	MDocument doc = new MDocument (mol);
	doc.setBondSetThickness(0, .125);
	renderer.setDoc(doc);

        g2.setRenderingHint(RenderingHints.KEY_RENDERING,
                            RenderingHints.VALUE_RENDER_QUALITY);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                            RenderingHints.VALUE_ANTIALIAS_ON);

	Shape shape;
	if (round) {
	    shape = new RoundRectangle2D.Double
		(x, y, (double)width, (double)height, 
		 width/4., height/4.);
	}
	else {
	    shape = new Rectangle (x, y, width, height);
	}

	// Clear the image so all pixels have zero alpha
	g2.setComposite(AlphaComposite.Clear);
	g2.fillRect(x, y, width, height);

	g2.setComposite(AlphaComposite.Src);
	g2.fill(shape);
	g2.setPaint(background);
	
	g2.setComposite(AlphaComposite.SrcAtop);
	g2.fill(shape);

        //g2.setComposite(AlphaComposite.SrcOver);
	Rectangle r = new Rectangle (x, y, width, height);
        renderer.setScale(.92 * renderer.maxScale(r));
        renderer.paint(g2, r);

        // MolPrinter reset these hints
        g2.setRenderingHint(RenderingHints.KEY_RENDERING,
                            RenderingHints.VALUE_RENDER_QUALITY);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                            RenderingHints.VALUE_ANTIALIAS_ON);
    }

    public BufferedImage createImage (String mol, int size) {
	return createImage (mol, size, true);
    }

    public BufferedImage createImage (String mol, int size, boolean round) {
	try {
	    return createImage 
		(new MolHandler (mol).getMolecule(), size, round);
	}
	catch (Exception ex) {
	    throw new IllegalArgumentException (ex);
	}
    }

    public BufferedImage createImage (Molecule mol, int size) {
	return createImage (mol, size, true);
    }

    public BufferedImage createImageShadow 
        (Molecule mol, int size, boolean round, 
         int radius, int offset, float translucent) {

        BufferedImage shadow = new BufferedImage 
	    (size, size, BufferedImage.TYPE_INT_ARGB);

        Graphics2D g2 = shadow.createGraphics();
        if (mol != null)
            renderMol (g2, mol, size, size, false);

        BufferedImage blur = createBlur (shadow, radius, translucent);

        g2.drawImage(blur, null, offset, offset);
        /*
        try {
            javax.imageio.ImageIO.write
                (shadow, "png", new java.io.File("_copy.png"));
            javax.imageio.ImageIO.write
                (blur, "png", new java.io.File("_blur.png"));
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        */

        g2.drawImage(shadow, null, 0, 0);
        if (round) {
            BufferedImage im = g2.getDeviceConfiguration()
                .createCompatibleImage
                (size, size, java.awt.Transparency.TRANSLUCENT);
            Graphics2D g = im.createGraphics();
            g.setRenderingHints(g2.getRenderingHints());
            
            Shape shape = new RoundRectangle2D.Double
                (0, 0, (double)size, (double)size, size/4., size/4.);
            g.setComposite(AlphaComposite.Clear);
            g.fillRect(0, 0, size, size);
            
            g.setComposite(AlphaComposite.Src);
            g.fill(shape);
            g.setPaint(background);
            
            g.setComposite(AlphaComposite.SrcAtop);
            g.fill(shape);
            
            g.drawImage(shadow, 0, 0, null);
            if (border) {
                g.setPaint(borderColor);
                g.setStroke(new BasicStroke (4.f));
                double arc = size/4. - 3.;
                shape = new RoundRectangle2D.Double
                    (1, 1, (double)size-3, (double)size-3, arc, arc);
                g.draw(shape);
            }
            g.dispose();

            shadow = im;
        }
        g2.dispose();

        return shadow;
    }

    public BufferedImage createImage (Molecule mol, int size, boolean round) {
        BufferedImage img;

        int radius = (int)(shadowRadius*size + .5f);

        if (shadow && radius > 0) {
            img = createImageShadow 
                (mol, size, round, radius, shadowOffset, shadowTranslucent);
        }
        else {
            img = new BufferedImage (size, size, BufferedImage.TYPE_INT_ARGB);
            Graphics2D g2 = img.createGraphics();
            renderMol (g2, mol, size, size, round);
            g2.dispose();
        }

	return img;
    }

    /*
     * This just makes a drop shadow from the given src image
     * 
     * Box blur is applied of given radius, 3 times to approximate gaussian.
     * 
     * Then the contrast is fixed to be completely black
     *
     * Then the image is painted at a certain semi-transparency.
     * 
     * TODO: make a better OpFilter that does all of this at once at the
     * pixel level. That means, do the box blur more specifically (possibly 4
     * times faster or more). Include the desired radius to box blur at, and iteration
     * number. 
     * 
     */
    static private BufferedImage createBlur
        (BufferedImage src, int rad, float opac){
        BufferedImage img = new BufferedImage
            (src.getWidth(), src.getHeight(), BufferedImage.TYPE_INT_ARGB);
        Graphics2D gf =  img.createGraphics();

        BufferedImageOp blur = OperCache.get(rad);
        if (blur == null) {
            int fsize=(rad*2+1)*(rad*2+1);
            int fwid=(rad*2+1);
            float[] blurKernel= new float[fsize];
            float dist=1/(float)fsize;
            for(int i=0;i<fsize;i++){
                blurKernel[i]=dist;
            }

            blur = new ConvolveOp
                (new Kernel(fwid, fwid, blurKernel), 
                 ConvolveOp.EDGE_NO_OP, gf.getRenderingHints());

            OperCache.putIfAbsent(rad, blur);
        }
				
        //theoretically, three box blurs are approximately
        //one gaussian.
        BufferedImage dst = GREYOP.filter(
    			blur.filter(
    			blur.filter(
    			blur.filter(src, null),
    			null), 
    			null),
    			null);

        gf.setComposite(AlphaComposite.getInstance
                        (AlphaComposite.SRC_OVER, opac));

        gf.drawImage(dst, 0, 0, null);
        gf.dispose();

        return img;
    }

    public static void main (String[] argv) throws Exception {
	System.setProperty("java.awt.headless", "true");

	MolRenderer renderer = new MolRenderer ();
	// some light background
	renderer.setBackground(new Color (0xff, 0xf2, 0xe0, 0xff));
        renderer.setBorderVisible(true);
	BufferedImage img = 
	    renderer.createImage("N[C@@H](CC1=CC(O)=C(O)C=C1)C(O)=O", 100);
	javax.imageio.ImageIO.write(img, "png", new java.io.File("test1.png"));

	renderer.setBackground(null);
	img = renderer.createImage("[CH3:6][C:5]1=[C:7](C)[CH:8]=[C:2]([NH:1]S(=O)(=O)[C:3]2=[C:2]([CH:8]=[C:7]3NC(=O)C[CH2:6][C:5]3=[CH:4]2)[N:1]4CCC4)[CH:3]=[CH:4]1", 300);
	javax.imageio.ImageIO.write(img, "png", new java.io.File("test2.png"));
    }
}
