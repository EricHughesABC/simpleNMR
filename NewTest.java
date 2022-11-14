import org.openscience.cdk.interfaces.IAtom;
import java.util.Locale;
import org.openscience.nmrshiftdb.PredictionTool;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.nmrshiftdb.util.AtomUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import java.io.Reader;
import org.openscience.cdk.io.MDLReader;
import java.io.FileReader;

// 
// Decompiled by Procyon v0.5.36
// 

public class NewTest
{
    public static void main(final String[] args) throws Exception {
        if (args.length == 0) {
            System.exit(1);
        }
        String solvent = "Unreported";
        if (args.length == 2 || args.length == 3) {
            solvent = args[1];
            if (!"Chloroform-D1 (CDCl3)Methanol-D4 (CD3OD)Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)".contains(solvent)) {
                System.out.println("solvent must be \"Chloroform-D1 (CDCl3)\", \"Methanol-D4 (CD3OD)\", or \"Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)\" (with quotation marks)!");
                System.exit(1);
            }
        }
        boolean use3d = true;
        if (args.length == 3 && args[2].equals("no3d")) {
            use3d = false;
        }
        final MDLReader mdlreader = new MDLReader((Reader)new FileReader(args[0]));
        final IAtomContainer mol = (IAtomContainer)mdlreader.read((IChemObject)DefaultChemObjectBuilder.getInstance().newInstance((Class)IAtomContainer.class, new Object[0]));
        AtomUtils.addAndPlaceHydrogens(mol);
        CDKHueckelAromaticityDetector.detectAromaticity(mol);
        final PredictionTool predictor = new PredictionTool();
        System.out.println("C,min,mean,max");

        for (int i = 0; i < mol.getAtomCount(); ++i) {
            final IAtom curAtom = mol.getAtom(i);
            if (curAtom.getAtomicNumber() == 6) {
                final float[] result = predictor.predict(mol, curAtom, use3d, solvent);
                System.out.format(Locale.US, "%3d,%8.2f,%8.2f,%8.2f\n", i + 1, result[0], result[1], result[2]);
            }
        }
        System.exit(0);
    }
}