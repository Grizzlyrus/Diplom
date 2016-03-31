package utils;

import matrix.SparseMatrix;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by Admin on 09.03.2016.
 */
public class MTXReader {
    private String typecode;
    private SparseMatrix matrix;

    public void read(String filename) throws java.io.IOException {
//        InputStream s = SparseMatrixEx2.class.getResourceAsStream(filename);
        InputStream s = new FileInputStream(filename);
        BufferedReader br = new BufferedReader(new InputStreamReader(s));

        // read type code initial line
        String line = br.readLine();
        typecode = line;

        // read comment lines if any
        boolean comment = true;

        boolean issymmetric = false;

        while (comment) {
            String headerLine[] = line.split(" ");
            if(headerLine[4].equals("symmetric")){
                issymmetric=true;
            }
            line = br.readLine();
            comment = line.startsWith("%");

        }

        // line now contains the size information which needs to be parsed
        String[] str = line.split("( )+");
        int nRows = (Integer.valueOf(str[0].trim())).intValue();
        int nColumns = (Integer.valueOf(str[1].trim())).intValue();
        int nNonZeros = (Integer.valueOf(str[2].trim())).intValue();

        // now we're into the data section
        matrix = new SparseMatrix(nRows, nColumns);
        while (true) {
            line = br.readLine();
            if (line == null)  break;
            str = line.split("( )+");
            int i = (Integer.valueOf(str[0].trim())).intValue();
            int j = (Integer.valueOf(str[1].trim())).intValue();
            double x = (Double.valueOf(str[2].trim())).doubleValue();

            matrix.setEntry(i - 1, j - 1, x);

            if(issymmetric){
                matrix.setEntry(j - 1, i - 1, x);
            }
        }
        br.close();
    }

    public String getTypeCode() {
        return this.typecode;
    }

    public SparseMatrix getMatrix(){
        return matrix;
    }
}
