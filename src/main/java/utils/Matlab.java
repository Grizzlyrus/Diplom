package utils;

import decomposition.LUDecomposition;
import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by Admin on 08.03.2016.
 */
public class Matlab {

    public static DenseVector full(Vector V) {
        if(!(V instanceof SparseVector)) {
            return (DenseVector)V;
        } else {
            int dim = V.getDim();
            int[] ir = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            double[] values = ArrayOperator.allocateVector(dim, 0.0D);

            for(int k = 0; k < ((SparseVector)V).getNNZ(); ++k) {
                values[ir[k]] = pr[k];
            }

            return new DenseVector(values);
        }
    }

    public static Matrix hilb(int m, int n) {
        DenseMatrix A = new DenseMatrix(m, n);
        double[][] data = A.getData();
        double[] A_i = (double[])null;

        for(int i = 0; i < m; ++i) {
            A_i = data[i];

            for(int j = 0; j < n; ++j) {
                A_i[j] = 1.0D / (double)(i + j + 1);
            }
        }

        return A;
    }

    public static Matrix inv(Matrix A) {
        int M = A.getRowDimension();
        int N = A.getColumnDimension();
        if(M != N) {
            System.err.println("Input should be a square matrix.");
        }

        LUDecomposition LUDecomp = new LUDecomposition(A);
        if(LUDecomp.det() == 0.0D) {
            System.err.println("The input matrix is not invertible.");
        }

        return LUDecomp.inverse();
    }

    public static SparseMatrix sparse(Matrix A) {
        if(!(A instanceof DenseMatrix)) {
            return (SparseMatrix)A;
        } else {
            boolean rIdx = false;
            boolean cIdx = false;
            int nzmax = 0;
            double value = 0.0D;
            TreeMap map = new TreeMap();
            int numRows = A.getRowDimension();
            int numColumns = A.getColumnDimension();
            double[][] data = ((DenseMatrix)A).getData();

            int var18;
            for(int ir = 0; ir < numColumns; ++ir) {
                var18 = ir;

                for(int jc = 0; jc < numRows; ++jc) {
                    value = data[jc][ir];
                    if(value != 0.0D) {
                        map.put(Pair.of(Integer.valueOf(var18), Integer.valueOf(jc)), Double.valueOf(value));
                        ++nzmax;
                    }
                }
            }

            int[] var19 = new int[nzmax];
            int[] var20 = new int[numColumns + 1];
            double[] pr = new double[nzmax];
            int k = 0;
            var20[0] = 0;
            int currentColumn = 0;

            for(Iterator var16 = map.entrySet().iterator(); var16.hasNext(); ++k) {
                Map.Entry entry = (Map.Entry)var16.next();
                int var17 = ((Integer)((Pair)entry.getKey()).second).intValue();
                var18 = ((Integer)((Pair)entry.getKey()).first).intValue();
                pr[k] = ((Double)entry.getValue()).doubleValue();
                var19[k] = var17;
                if(currentColumn < var18) {
                    var20[currentColumn + 1] = k;
                    ++currentColumn;
                }
            }

            while(currentColumn < numColumns) {
                var20[currentColumn + 1] = k;
                ++currentColumn;
            }

            var20[numColumns] = k;
            return SparseMatrix.createSparseMatrixByCSCArrays(var19, var20, pr, numRows, numColumns, nzmax);
        }
    }

    public static DenseMatrix full(Matrix S) {
        if(!(S instanceof SparseMatrix)) {
            return (DenseMatrix)S;
        } else {
            int M = S.getRowDimension();
            int N = S.getColumnDimension();
            double[][] data = new double[M][];
            int[] ic = ((SparseMatrix)S).getIc();
            int[] jr = ((SparseMatrix)S).getJr();
            int[] valCSRIndices = ((SparseMatrix)S).getValCSRIndices();
            double[] pr = ((SparseMatrix)S).getPr();

            for(int i = 0; i < M; ++i) {
                double[] rowData = ArrayOperator.allocateVector(N, 0.0D);

                for(int k = jr[i]; k < jr[i + 1]; ++k) {
                    rowData[ic[k]] = pr[valCSRIndices[k]];
                }

                data[i] = rowData;
            }

            return new DenseMatrix(data);
        }
    }

    public static Vector[] sparseMatrix2SparseRowVectors(Matrix S) {
        if(!(S instanceof SparseMatrix)) {
            System.err.println("SparseMatrix input is expected.");
            System.exit(1);
        }

        int M = S.getRowDimension();
        int N = S.getColumnDimension();
        Vector[] Vs = new Vector[M];
        int[] ic = ((SparseMatrix)S).getIc();
        int[] jr = ((SparseMatrix)S).getJr();
        double[] pr = ((SparseMatrix)S).getPr();
        int[] valCSRIndices = ((SparseMatrix)S).getValCSRIndices();
        int[] indices = (int[])null;
        double[] values = (double[])null;
        boolean nnz = false;
        int dim = N;

        for(int r = 0; r < M; ++r) {
            int var15 = jr[r + 1] - jr[r];
            indices = new int[var15];
            values = new double[var15];
            int idx = 0;

            for(int k = jr[r]; k < jr[r + 1]; ++k) {
                indices[idx] = ic[k];
                values[idx] = pr[valCSRIndices[k]];
                ++idx;
            }

            Vs[r] = new SparseVector(indices, values, var15, dim);
        }

        return Vs;
    }

    public static Vector[] sparseMatrix2SparseColumnVectors(Matrix S) {
        if(!(S instanceof SparseMatrix)) {
            System.err.println("SparseMatrix input is expected.");
            System.exit(1);
        }

        int M = S.getRowDimension();
        int N = S.getColumnDimension();
        Vector[] Vs = new Vector[N];
        int[] ir = ((SparseMatrix)S).getIr();
        int[] jc = ((SparseMatrix)S).getJc();
        double[] pr = ((SparseMatrix)S).getPr();
        int[] indices = (int[])null;
        double[] values = (double[])null;
        boolean nnz = false;
        int dim = M;

        for(int c = 0; c < N; ++c) {
            int var14 = jc[c + 1] - jc[c];
            indices = new int[var14];
            values = new double[var14];
            int idx = 0;

            for(int k = jc[c]; k < jc[c + 1]; ++k) {
                indices[idx] = ir[k];
                values[idx] = pr[k];
                ++idx;
            }

            Vs[c] = new SparseVector(indices, values, var14, dim);
        }

        return Vs;
    }

    public static Matrix sparseRowVectors2SparseMatrix(Vector[] Vs) {
        int nnz = 0;
        int numColumns = Vs.length > 0?Vs[0].getDim():0;

        int numRows;
        for(numRows = 0; numRows < Vs.length; ++numRows) {
            if(!(Vs[numRows] instanceof SparseVector)) {
                Printer.fprintf("Vs[%d] should be a sparse vector.%n", new Object[]{Integer.valueOf(numRows)});
                System.exit(1);
            }
//TODO
            nnz += ((SparseVector)Vs[numRows]).getNNZ();
            if(numColumns != Vs[numRows].getDim()) {
                Printer.fprintf("Vs[%d]\'s dimension doesn\'t match.%n", new Object[]{Integer.valueOf(numRows)});
                System.exit(1);
            }
        }

        numRows = Vs.length;
        int[] ic = new int[nnz];
        int[] jr = new int[numRows + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        boolean cIdx = true;
        int cnt = 0;
        jr[0] = 0;
        int currentRow = 0;
//TODO
        for(int i = 0; i < numRows; ++i) {
            int[] indices = ((SparseVector)Vs[i]).getIr();
            double[] values = ((SparseVector)Vs[i]).getPr();
            int nnz1 = ((SparseVector)Vs[i]).getNNZ();
//            nnz = ((SparseVector)Vs[i]).getNNZ();

            for(int k = 0; k < nnz1; ++k) {
                int var17 = indices[k];
                int var16 = i;
                pr[cnt] = values[k];

                for(ic[cnt] = var17; currentRow < var16; ++currentRow) {
                    jr[currentRow + 1] = cnt;
                }

                ++cnt;
            }
        }

        while(currentRow < numRows) {
            jr[currentRow + 1] = cnt;
            ++currentRow;
        }

        return SparseMatrix.createSparseMatrixByCSRArrays(ic, jr, pr, numRows, numColumns, nnz);
    }

    public static Matrix sparseColumnVectors2SparseMatrix(Vector[] Vs) {
        int nnz = 0;
        int numRows = Vs.length > 0?Vs[0].getDim():0;

        int numColumns;
        for(numColumns = 0; numColumns < Vs.length; ++numColumns) {
            if(!(Vs[numColumns] instanceof SparseVector)) {
                Printer.fprintf("Vs[%d] should be a sparse vector.%n", new Object[]{Integer.valueOf(numColumns)});
                System.exit(1);
            }

            nnz += ((SparseVector)Vs[numColumns]).getNNZ();
            if(numRows != Vs[numColumns].getDim()) {
                Printer.fprintf("Vs[%d]\'s dimension doesn\'t match.%n", new Object[]{Integer.valueOf(numColumns)});
                System.exit(1);
            }
        }

        numColumns = Vs.length;
        int[] ir = new int[nnz];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        boolean cIdx = true;
        int k = 0;
        jc[0] = 0;
        int currentColumn = 0;

        for(int c = 0; c < numColumns; ++c) {
            int[] indices = ((SparseVector)Vs[c]).getIr();
            double[] values = ((SparseVector)Vs[c]).getPr();
            nnz = ((SparseVector)Vs[c]).getNNZ();

            for(int r = 0; r < nnz; ++r) {
                int var16 = indices[r];
                int var17 = c;
                pr[k] = values[r];

                for(ir[k] = var16; currentColumn < var17; ++currentColumn) {
                    jc[currentColumn + 1] = k;
                }

                ++k;
            }
        }

        while(currentColumn < numColumns) {
            jc[currentColumn + 1] = k;
            ++currentColumn;
        }

        return SparseMatrix.createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nnz);
    }

}
