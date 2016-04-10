package utils;

import matrix.SparseMatrix;

import java.io.*;

/**
 * Created by Admin on 03.04.2016.
 */
public class MTXWriter {

    public void write(String filename, SparseMatrix A) throws IOException{
        OutputStream s = new FileOutputStream(filename);
        BufferedWriter bwr = new BufferedWriter(new OutputStreamWriter(s));
        bwr.write("%%MatrixMarket matrix coordinate real general\n");
        int nRows = A.getRowDimension();
        int nColumns = A.getColumnDimension();
        int nNonZeros = A.getNNZ();

        String line = new String();
        line=String.valueOf(nRows)+" "+String.valueOf(nColumns) + " "+ String.valueOf(nNonZeros)+"\n";
        bwr.write(line);

        int[] ir = A.getIr();
        int[] jc = A.getJc();
        double[] pr = A.getPr();

        for (int i = 0; i < nColumns; i++) {
//            line = "";
            for (int j = jc[i]; j < jc[i+1]; j++) {
                line=String.valueOf(ir[j]+1)+" "+String.valueOf(i+1)+" "+String.valueOf(pr[j])+"\n";

                bwr.write(line);
            }
        }
        bwr.close();
    }
}
