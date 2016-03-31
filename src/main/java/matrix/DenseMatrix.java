package matrix;

import decomposition.QRDecomposition;
import utils.ArrayOperator;
import utils.Matlab;
import utils.Printer;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.io.Serializable;
import java.util.Arrays;

/**
 * Created by Admin on 08.03.2016.
 */

public class DenseMatrix implements Matrix, Serializable {
    private static final long serialVersionUID = 6821454132254344419L;
    private int M;
    private int N;
    private double[][] data;

    public DenseMatrix() {
        this.M = 0;
        this.N = 0;
        this.data = null;
    }

    public DenseMatrix(double v) {
        this.M = 1;
        this.N = 1;
        this.data = new double[][]{{v}};
    }

    public DenseMatrix(int M, int N) {
        this.data = new double[M][];

        for(int i = 0; i < M; ++i) {
            this.data[i] = new double[N];

            for(int j = 0; j < N; ++j) {
                this.data[i][j] = 0.0D;
            }
        }

        this.M = M;
        this.N = N;
    }

    public DenseMatrix(int[] size) {
        if(size.length != 2) {
            System.err.println("The input integer array should have exactly two entries!");
            System.exit(1);
        }

        int M = size[0];
        int N = size[1];
        this.data = new double[M][];

        for(int i = 0; i < M; ++i) {
            this.data[i] = new double[N];

            for(int j = 0; j < N; ++j) {
                this.data[i][j] = 0.0D;
            }
        }

        this.M = M;
        this.N = N;
    }

    public DenseMatrix(double[][] data) {
        this.data = data;
        this.M = data.length;
        this.N = this.M > 0?data[0].length:0;
    }

    public DenseMatrix(double[] data, int dim) {
        if(dim == 1) {
            this.M = data.length;
            this.N = 1;
            this.data = new double[this.M][];

            for(int i = 0; i < this.M; ++i) {
                this.data[i] = new double[this.N];
                this.data[i][0] = data[i];
            }
        } else if(dim == 2) {
            this.M = 1;
            this.N = data.length;
            this.data = new double[this.M][];
            this.data[0] = data;
        }

    }

    public DenseMatrix(int M, int N, double v) {
        this.data = new double[M][];

        for(int i = 0; i < M; ++i) {
            this.data[i] = new double[N];

            for(int j = 0; j < N; ++j) {
                this.data[i][j] = v;
            }
        }

        this.M = M;
        this.N = N;
    }

    public DenseMatrix(int[] size, double v) {
        if(size.length != 2) {
            System.err.println("The input integer array should have exactly two entries!");
            System.exit(1);
        }

        int M = size[0];
        int N = size[1];
        this.data = new double[M][];

        for(int i = 0; i < M; ++i) {
            this.data[i] = new double[N];

            for(int j = 0; j < N; ++j) {
                this.data[i][j] = v;
            }
        }

        this.M = M;
        this.N = N;
    }

    public double[][] getData() {
        return this.data;
    }

    public int getRowDimension() {
        return this.M;
    }

    public int getColumnDimension() {
        return this.N;
    }

    public Matrix mtimes(Matrix A) {
        DenseMatrix res = null;
        double[][] resData = new double[this.M][];
        int NA = A.getColumnDimension();

        for(int rowData = 0; rowData < this.M; ++rowData) {
            resData[rowData] = new double[NA];
        }

        double[] var15 = (double[])null;
        int i1;
        if(A instanceof DenseMatrix) {
            double[][] ir = ((DenseMatrix)A).getData();
            double[] jc = new double[A.getRowDimension()];
            double pr = 0.0D;

            for(int s = 0; s < NA; ++s) {
                int i;
                for(i = 0; i < A.getRowDimension(); ++i) {
                    jc[i] = ir[i][s];
                }

                for(i = 0; i < this.M; ++i) {
                    var15 = this.data[i];
                    pr = 0.0D;

                    for(i1 = 0; i1 < this.N; ++i1) {
                        pr += var15[i1] * jc[i1];
                    }

                    resData[i][s] = pr;
                }
            }
        } else if(A instanceof SparseMatrix) {
            int[] var16 = (int[])null;
            int[] var17 = (int[])null;
            double[] var18 = (double[])null;
            var16 = ((SparseMatrix)A).getIr();
            var17 = ((SparseMatrix)A).getJc();
            var18 = ((SparseMatrix)A).getPr();
            boolean r = true;
            double var20 = 0.0D;

            for(i1 = 0; i1 < this.M; ++i1) {
                var15 = this.data[i1];

                for(int j = 0; j < NA; ++j) {
                    var20 = 0.0D;

                    for(int k = var17[j]; k < var17[j + 1]; ++k) {
                        int var19 = var16[k];
                        var20 += var15[var19] * var18[k];
                    }

                    resData[i1][j] = var20;
                }
            }
        }

        res = new DenseMatrix(resData);
        return res;
    }


    public double getEntry(int r, int c) {
        return this.data[r][c];
    }

    public void setEntry(int r, int c, double v) {
        this.data[r][c] = v;
    }

    public Matrix transpose() {
        double[][] resData = new double[this.N][];

        for(int i = 0; i < this.N; ++i) {
            resData[i] = new double[this.M];

            for(int j = 0; j < this.M; ++j) {
                resData[i][j] = this.data[j][i];
            }
        }

        return new DenseMatrix(resData);
    }

    public Matrix plus(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            DenseMatrix res = (DenseMatrix)this.copy();
            double[] pr;
            int r;
            int j;
            if(A instanceof DenseMatrix) {
                double[][] ir = ((DenseMatrix)A).getData();
                double[] jc = (double[])null;
                pr = (double[])null;

                for(r = 0; r < this.M; ++r) {
                    pr = res.data[r];
                    jc = ir[r];

                    for(j = 0; j < this.N; ++j) {
                        pr[j] += jc[j];
                    }
                }
            } else if(A instanceof SparseMatrix) {
                int[] var9 = (int[])null;
                int[] var10 = (int[])null;
                pr = (double[])null;
                var9 = ((SparseMatrix)A).getIr();
                var10 = ((SparseMatrix)A).getJc();
                pr = ((SparseMatrix)A).getPr();
                boolean var11 = true;

                for(j = 0; j < A.getColumnDimension(); ++j) {
                    for(int k = var10[j]; k < var10[j + 1]; ++k) {
                        r = var9[k];
                        res.data[r][j] += pr[k];
                    }
                }
            }

            return res;
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix minus(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            DenseMatrix res = (DenseMatrix)this.copy();
            if(A instanceof DenseMatrix) {
                for(int ir = 0; ir < this.M; ++ir) {
                    for(int jc = 0; jc < this.N; ++jc) {
                        res.data[ir][jc] -= ((DenseMatrix)A).data[ir][jc];
                    }
                }
            } else if(A instanceof SparseMatrix) {
                int[] var9 = (int[])null;
                int[] var10 = (int[])null;
                double[] pr = (double[])null;
                var9 = ((SparseMatrix)A).getIr();
                var10 = ((SparseMatrix)A).getJc();
                pr = ((SparseMatrix)A).getPr();
                boolean r = true;

                for(int j = 0; j < A.getColumnDimension(); ++j) {
                    for(int k = var10[j]; k < var10[j + 1]; ++k) {
                        int var11 = var9[k];
                        res.data[var11][j] -= pr[k];
                    }
                }
            }

            return res;
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix times(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            double[][] resData = ArrayOperator.allocate2DArray(this.M, this.N, 0.0D);
            double[] pr;
            int j;
            int k;
            if(A instanceof DenseMatrix) {
                double[][] ir = ((DenseMatrix)A).getData();
                double[] jc = (double[])null;
                pr = (double[])null;
                double[] r = (double[])null;

                for(j = 0; j < this.M; ++j) {
                    double[] var10000 = this.data[j];
                    jc = ir[j];
                    r = resData[j];
                    pr = this.data[j];

                    for(k = 0; k < this.N; ++k) {
                        r[k] = pr[k] * jc[k];
                    }
                }
            } else if(A instanceof SparseMatrix) {
                int[] var9 = (int[])null;
                int[] var10 = (int[])null;
                pr = (double[])null;
                var9 = ((SparseMatrix)A).getIr();
                var10 = ((SparseMatrix)A).getJc();
                pr = ((SparseMatrix)A).getPr();
                boolean var11 = true;

                for(j = 0; j < A.getColumnDimension(); ++j) {
                    for(k = var10[j]; k < var10[j + 1]; ++k) {
                        int var12 = var9[k];
                        resData[var12][j] = this.data[var12][j] * pr[k];
                    }
                }
            }

            return new DenseMatrix(resData);
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix times(double v) {
        DenseMatrix res = (DenseMatrix)this.copy();

        for(int i = 0; i < this.M; ++i) {
            for(int j = 0; j < this.N; ++j) {
                res.data[i][j] *= v;
            }
        }

        return res;
    }

    public Matrix copy() {
        DenseMatrix res = new DenseMatrix();
        res.M = this.M;
        res.N = this.N;
        res.data = new double[this.M][];

        for(int i = 0; i < this.M; ++i) {
            res.data[i] = (double[])this.data[i].clone();
        }

        return res;
    }

    public Matrix clone() {
        return this.copy();
    }

    public Matrix plus(double v) {
        DenseMatrix res = (DenseMatrix)this.copy();

        for(int i = 0; i < this.M; ++i) {
            for(int j = 0; j < this.N; ++j) {
                res.data[i][j] += v;
            }
        }

        return res;
    }

    public Matrix minus(double v) {
        DenseMatrix res = (DenseMatrix)this.copy();

        for(int i = 0; i < this.M; ++i) {
            for(int j = 0; j < this.N; ++j) {
                res.data[i][j] -= v;
            }
        }

        return res;
    }

    public Vector operate(Vector b) {
        if(this.N != b.getDim()) {
            System.err.println("Dimension does not match.");
            System.exit(1);
        }

        double[] V = new double[this.M];
        if(b instanceof DenseVector) {
            ArrayOperator.operate(V, this.data, ((DenseVector)b).getPr());
        } else if(b instanceof SparseVector) {
            int[] ir = ((SparseVector)b).getIr();
            double[] pr = ((SparseVector)b).getPr();
            int nnz = ((SparseVector)b).getNNZ();
            boolean idx = false;
            double[] row_i = (double[])null;

            for(int i = 0; i < this.M; ++i) {
                row_i = this.data[i];
                double s = 0.0D;

                for(int k = 0; k < nnz; ++k) {
                    int var12 = ir[k];
                    s += row_i[var12] * pr[k];
                }

                V[i] = s;
            }
        }

        return new DenseVector(V);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        if(this.data == null) {
            sb.append("Empty matrix." + System.lineSeparator());
            return sb.toString();
        } else {
            byte p = 4;

            for(int i = 0; i < this.getRowDimension(); ++i) {
                sb.append("  ");

                for(int j = 0; j < this.getColumnDimension(); ++j) {
                    String valueString = "";
                    double v = this.getEntry(i, j);
                    int rv = (int)Math.round(v);
                    if(v != (double)rv) {
                        valueString = Printer.sprintf(Printer.sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                    } else {
                        valueString = Printer.sprintf("%d", new Object[]{Integer.valueOf(rv)});
                    }

                    sb.append(Printer.sprintf(Printer.sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
                    sb.append("  ");
                }

                sb.append(System.lineSeparator());
            }

            return sb.toString();
        }
    }

    public void clear() {
        ArrayOperator.clearMatrix(this.data);
    }

    public Matrix getSubMatrix(int startRow, int endRow, int startColumn, int endColumn) {
        int nRow = endRow - startRow + 1;
        int nCol = endColumn - startColumn + 1;
        double[][] resData = new double[nRow][];
        double[] resRow = (double[])null;
        double[] thisRow = (double[])null;
        int r = 0;

        for(int i = startRow; r < nRow; ++i) {
            resRow = new double[nCol];
            thisRow = this.data[i];
            System.arraycopy(thisRow, startColumn, resRow, 0, nCol);
            resData[r] = resRow;
            ++r;
        }

        return new DenseMatrix(resData);
    }

    public Matrix getSubMatrix(int[] selectedRows, int[] selectedColumns) {
        int nRow = selectedRows.length;
        int nCol = selectedColumns.length;
        double[][] resData = new double[nRow][];
        double[] resRow = (double[])null;
        double[] thisRow = (double[])null;

        for(int r = 0; r < nRow; ++r) {
            resRow = new double[nCol];
            thisRow = this.data[selectedRows[r]];

            for(int c = 0; c < nCol; ++c) {
                resRow[c] = thisRow[selectedColumns[c]];
            }

            resData[r] = resRow;
        }

        return new DenseMatrix(resData);
    }

    public Matrix getColumnMatrix(int c) {
        DenseMatrix res = new DenseMatrix(this.M, 1);
        double[][] resData = res.data;

        for(int i = 0; i < this.M; ++i) {
            resData[i][0] = this.data[i][c];
        }

        return res;
    }

    public Vector getColumnVector(int c) {
        DenseVector res = new DenseVector(this.M);
        double[] pr = res.getPr();

        for(int i = 0; i < this.M; ++i) {
            pr[i] = this.data[i][c];
        }

        return res;
    }

    public Matrix getRowMatrix(int r) {
        return new DenseMatrix(this.data[r], 2);
    }

    public Vector getRowVector(int r) {
        return new DenseVector(this.data[r]);
    }

    public Matrix getRows(int startRow, int endRow) {
        int numRows = endRow - startRow + 1;
        double[][] resData = new double[numRows][];
        int r = startRow;

        for(int i = 0; r <= endRow; ++i) {
            resData[i] = (double[])this.data[r].clone();
            ++r;
        }

        return new DenseMatrix(resData);
    }

    public Matrix getRows(int... selectedRows) {
        int numRows = selectedRows.length;
        double[][] resData = new double[numRows][];

        for(int i = 0; i < numRows; ++i) {
            resData[i] = (double[])this.data[selectedRows[i]].clone();
        }

        return new DenseMatrix(resData);
    }

    public Vector[] getRowVectors(int startRow, int endRow) {
        int numRows = endRow - startRow + 1;
        DenseVector[] res = new DenseVector[numRows];
        int r = startRow;

        for(int i = 0; r <= endRow; ++i) {
            res[i] = new DenseVector(this.data[r]);
            ++r;
        }

        return res;
    }

    public Vector[] getRowVectors(int... selectedRows) {
        int numRows = selectedRows.length;
        DenseVector[] res = new DenseVector[numRows];

        for(int i = 0; i < numRows; ++i) {
            res[i] = new DenseVector(this.data[selectedRows[i]]);
        }

        return res;
    }

    public Matrix getColumns(int startColumn, int endColumn) {
        int nRow = this.M;
        int nCol = endColumn - startColumn + 1;
        double[][] resData = new double[nRow][];
        double[] resRow = (double[])null;
        double[] thisRow = (double[])null;

        for(int r = 0; r < nRow; ++r) {
            resRow = new double[nCol];
            thisRow = this.data[r];
            System.arraycopy(thisRow, startColumn, resRow, 0, nCol);
            resData[r] = resRow;
        }

        return new DenseMatrix(resData);
    }

    public Matrix getColumns(int... selectedColumns) {
        int nRow = this.M;
        int nCol = selectedColumns.length;
        double[][] resData = new double[nRow][];
        double[] resRow = (double[])null;
        double[] thisRow = (double[])null;

        for(int r = 0; r < nRow; ++r) {
            resRow = new double[nCol];
            thisRow = this.data[r];

            for(int c = 0; c < nCol; ++c) {
                resRow[c] = thisRow[selectedColumns[c]];
            }

            resData[r] = resRow;
        }

        return new DenseMatrix(resData);
    }

    public Vector[] getColumnVectors(int startColumn, int endColumn) {
        int numColumns = endColumn - startColumn + 1;
        DenseVector[] res = new DenseVector[numColumns];
        int c = startColumn;

        for(int i = 0; c <= endColumn; ++i) {
            res[i] = (DenseVector) this.getColumnVector(c);
            ++c;
        }

        return res;
    }

    public Vector[] getColumnVectors(int... selectedColumns) {
        int numColumns = selectedColumns.length;
        DenseVector[] res = new DenseVector[numColumns];

        for(int j = 0; j < numColumns; ++j) {
            res[j] = (DenseVector) this.getColumnVector(selectedColumns[j]);
        }

        return res;
    }

    public void setRowMatrix(int r, Matrix A) {
        if(A.getRowDimension() != 1) {
            Printer.err("Input matrix should be a row matrix.");
//            Utility.exit(1);
        }

        double[] thisRow = this.data[r];
        if(A instanceof DenseMatrix) {
            double[] jc = ((DenseMatrix)A).data[0];
            System.arraycopy(jc, 0, thisRow, 0, this.N);
        } else if(A instanceof SparseMatrix) {
            int[] var7 = ((SparseMatrix)A).getJc();
            double[] pr = ((SparseMatrix)A).getPr();

            for(int j = 0; j < this.N; ++j) {
                if(var7[j + 1] == var7[j]) {
                    thisRow[j] = 0.0D;
                } else {
                    thisRow[j] = pr[var7[j]];
                }
            }
        }

    }

    public void setRowVector(int r, Vector V) {
        double[] thisRow = this.data[r];
        if(V instanceof DenseVector) {
            double[] ir = ((DenseVector)V).getPr();
            System.arraycopy(ir, 0, thisRow, 0, this.N);
        } else if(V instanceof SparseVector) {
            int[] var11 = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            int nnz = ((SparseVector)V).getNNZ();
            int lastIdx = -1;
            boolean currentIdx = false;

            int j;
            for(j = 0; j < nnz; ++j) {
                int var12 = var11[j];

                for(int j1 = lastIdx + 1; j1 < var12; ++j1) {
                    thisRow[j1] = 0.0D;
                }

                thisRow[var12] = pr[j];
                lastIdx = var12;
            }

            for(j = lastIdx + 1; j < this.N; ++j) {
                thisRow[j] = 0.0D;
            }
        }

    }

    public void setColumnMatrix(int c, Matrix A) {
        if(A.getColumnDimension() != 1) {
            Printer.err("Input matrix should be a column matrix.");
//            Utility.exit(1);
        }

        int ir;
        if(A instanceof DenseMatrix) {
            double[][] jc = ((DenseMatrix)A).data;

            for(ir = 0; ir < this.M; ++ir) {
                this.data[ir][c] = jc[ir][0];
            }
        } else if(A instanceof SparseMatrix) {
            int[] var10 = ((SparseMatrix)A).getJc();
            if(var10[1] == 0) {
                for(ir = 0; ir < this.M; ++ir) {
                    this.data[ir][c] = 0.0D;
                }

                return;
            }

            int[] var11 = ((SparseMatrix)A).getIr();
            double[] pr = ((SparseMatrix)A).getPr();
            int lastIdx = -1;
            boolean currentIdx = false;

            int i;
            for(i = 0; i < var10[1]; ++i) {
                int var12 = var11[i];

                for(int i1 = lastIdx + 1; i1 < var12; ++i1) {
                    this.data[i1][c] = 0.0D;
                }

                this.data[var12][c] = pr[i];
                lastIdx = var12;
            }

            for(i = lastIdx + 1; i < this.M; ++i) {
                this.data[i][c] = 0.0D;
            }
        }

    }

    public void setColumnVector(int c, Vector V) {
        if(V instanceof DenseVector) {
            double[] ir = ((DenseVector)V).getPr();

            for(int pr = 0; pr < this.M; ++pr) {
                this.data[pr][c] = ir[pr];
            }
        } else if(V instanceof SparseVector) {
            int[] var10 = ((SparseVector)V).getIr();
            double[] var11 = ((SparseVector)V).getPr();
            int nnz = ((SparseVector)V).getNNZ();
            int lastIdx = -1;
            boolean currentIdx = false;

            int i;
            for(i = 0; i < nnz; ++i) {
                int var12 = var10[i];

                for(int i1 = lastIdx + 1; i1 < var12; ++i1) {
                    this.data[i1][c] = 0.0D;
                }

                this.data[var12][c] = var11[i];
                lastIdx = var12;
            }

            for(i = lastIdx + 1; i < this.M; ++i) {
                this.data[i][c] = 0.0D;
            }
        }

    }

    public double norm1(){
        double cur;
        double max = 0;
        for (int i = 0; i < this.N; i++) {
            cur = 0;
            for (int j = 0; j < this.M; j++) {
                cur+=Math.abs(data[i][j]);
            }
            max=(cur>max)?cur:max;
        }
        return max;
    }

    public double norm2(){
        double epsilon = 10E-18;
        Matrix AtA = this.copy().transpose().mtimes(this.copy());
        int m = AtA.getRowDimension();
        double temp;
        double diag[]= new double[m];
        Matrix QR[] = new Matrix[3];
        do{
            temp=0;
            QR=QRDecomposition.decompose(AtA,false);
            AtA = QR[1].mtimes(QR[0]);
            for (int i = 0; i < m; i++) {
                for (int j = m-1; j > i  ; j--) {
                    temp+=AtA.getEntry(j,i);
                }
            }

        }while (temp>epsilon);

        for (int i = 0; i < m; i++) {
            diag[i] = AtA.getEntry(i,i);
        }
        Arrays.sort(diag);
        return Math.sqrt(diag[m-1]);
    }

    public double normInf(){
        double cur;
        double max = 0;
        for (int i = 0; i < this.M; i++) {
            cur = 0;
            for (int j = 0; j < this.N; j++) {
                cur+=Math.abs(data[i][j]);
            }
            max=(cur>max)?cur:max;
        }
        return max;
    }

    public double normE(){
        double temp=0;
        for (int i = 0; i < this.getRowDimension(); i++) {
            for (int j = 0; j < this.getColumnDimension(); j++) {
                temp+=data[i][j]*data[i][j];
            }
        }
        return Math.sqrt(temp);
    }

    public double cond1(){
        return this.norm1()*Matlab.inv(this).norm1();
    }

    public double condInf(){
        return this.normInf()*Matlab.inv(this).normInf();
    }

    public double condE(){
        return this.normE()*Matlab.inv(this).normE();
    }

    public double cond2(){
        return this.norm2()*Matlab.inv(this).norm2();
    }

}