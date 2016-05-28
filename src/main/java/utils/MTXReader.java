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
        InputStream s = new FileInputStream(filename);
        BufferedReader br = new BufferedReader(new InputStreamReader(s));

        String line = br.readLine();
        typecode = line;

        boolean comment = true;

        boolean issymmetric = false;
        boolean pattern = false;

        String headerLine[] = line.split(" ");
        if(headerLine[4].equals("symmetric")){
            issymmetric=true;
        }
        if(headerLine[3].equals("pattern")){
            pattern=true;
        }

        while (comment) {
            line = br.readLine();
            comment = line.startsWith("%");

        }


        String[] str = line.split("( )+");
        int nRows = (Integer.valueOf(str[0].trim())).intValue();
        int nColumns = (Integer.valueOf(str[1].trim())).intValue();
        int nNonZeros = (Integer.valueOf(str[2].trim())).intValue();



        matrix = new SparseMatrix(nRows, nColumns);
        while (true) {
            double x;
            line = br.readLine();
            if (line == null)  break;
            str = line.split("( )+");
            int i = (Integer.valueOf(str[0].trim())).intValue();
            int j = (Integer.valueOf(str[1].trim())).intValue();
            if(!pattern) {
                x = (Double.valueOf(str[2].trim())).doubleValue();
            }else{
//                x = 2.0*Math.random()*10E9-10E9;
                if(Math.abs(i-j)<10) {
                    x = Math.random() * 10E-4;
                }else{
                    x = Math.random()*2.0*10E2-10E2;
                }
            }
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
