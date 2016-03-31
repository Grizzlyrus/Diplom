package matrix;

import decomposition.QRDecomposition;
import utils.ArrayOperator;
import utils.Matlab;
import utils.Pair;
import utils.Printer;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.io.Serializable;
import java.util.*;

/**
 * Created by Admin on 08.03.2016.
 */
public class SparseMatrix implements Matrix, Serializable {
    private static final long serialVersionUID = 404718895052720649L;
    private int[] ir;
    private int[] jc;
    private int[] ic;
    private int[] jr;
    private double[] pr;
    private int[] valCSRIndices;
    private int nnz;
    private int nzmax;
    private int M;
    private int N;

    private SparseMatrix() {
        this.M = 0;
        this.N = 0;
        this.nzmax = 0;
        this.nnz = 0;
    }

    public SparseMatrix(int M, int N) {
        this.M = M;
        this.N = N;
        this.nzmax = 0;
        this.nnz = 0;
        this.jc = new int[N + 1];

        int i;
        for(i = 0; i < N + 1; ++i) {
            this.jc[i] = 0;
        }

        this.jr = new int[M + 1];

        for(i = 0; i < M + 1; ++i) {
            this.jr[i] = 0;
        }

        this.ir = new int[0];
        this.pr = new double[0];
        this.ic = new int[0];
        this.valCSRIndices = new int[0];
    }

    public SparseMatrix(SparseMatrix A) {
        this.ir = A.ir;
        this.jc = A.jc;
        this.pr = A.pr;
        this.ic = A.ic;
        this.jr = A.jr;
        this.valCSRIndices = A.valCSRIndices;
        this.M = A.M;
        this.N = A.N;
        this.nzmax = A.nzmax;
        this.nnz = A.nnz;
    }

    public SparseMatrix(int[] rIndices, int[] cIndices, double[] values, int numRows, int numColumns, int nzmax) {
        SparseMatrix temp = createSparseMatrix(rIndices, cIndices, values, numRows, numColumns, nzmax);
        this.assignSparseMatrix(temp);
    }

    public void assignSparseMatrix(SparseMatrix A) {
        this.ir = (int[])A.ir.clone();
        this.jc = (int[])A.jc.clone();
        this.pr = (double[])A.pr.clone();
        this.ic = (int[])A.ic.clone();
        this.jr = (int[])A.jr.clone();
        this.valCSRIndices = (int[])A.valCSRIndices.clone();
        this.M = A.M;
        this.N = A.N;
        this.nzmax = A.nzmax;
        this.nnz = A.nnz;
    }

    public static SparseMatrix createSparseMatrix(TreeMap<Pair<Integer, Integer>, Double> inputMap, int numRows, int numColumns) {
        TreeMap map = new TreeMap();
        int nzmax = 0;
        Iterator jc = inputMap.entrySet().iterator();

        while(jc.hasNext()) {
            Map.Entry ir = (Map.Entry)jc.next();
            if(((Double)ir.getValue()).doubleValue() != 0.0D) {
                map.put(Pair.of((Integer)((Pair)ir.getKey()).second, (Integer)((Pair)ir.getKey()).first), (Double)ir.getValue());
                ++nzmax;
            }
        }

        int[] var14 = new int[nzmax];
        int[] var15 = new int[numColumns + 1];
        double[] pr = new double[nzmax];
        boolean rIdx = true;
        boolean cIdx = true;
        int k = 0;
        var15[0] = 0;
        int currentColumn = 0;

        for(Iterator var13 = map.entrySet().iterator(); var13.hasNext(); ++k) {
            Map.Entry entry = (Map.Entry)var13.next();
            int var16 = ((Integer)((Pair)entry.getKey()).second).intValue();
            int var17 = ((Integer)((Pair)entry.getKey()).first).intValue();
            pr[k] = ((Double)entry.getValue()).doubleValue();

            for(var14[k] = var16; currentColumn < var17; ++currentColumn) {
                var15[currentColumn + 1] = k;
            }
        }

        while(currentColumn < numColumns) {
            var15[currentColumn + 1] = k;
            ++currentColumn;
        }

        var15[numColumns] = k;
        return createSparseMatrixByCSCArrays(var14, var15, pr, numRows, numColumns, nzmax);
    }

    public static SparseMatrix createSparseMatrix(int[] rIndices, int[] cIndices, double[] values, int numRows, int numColumns, int nzmax) {
        boolean k = true;
        TreeMap map = new TreeMap();

        int var16;
        for(var16 = 0; var16 < values.length; ++var16) {
            if(values[var16] != 0.0D) {
                map.put(Pair.of(Integer.valueOf(cIndices[var16]), Integer.valueOf(rIndices[var16])), Double.valueOf(values[var16]));
            }
        }

        int[] ir = new int[nzmax];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nzmax];
        boolean rIdx = true;
        boolean cIdx = true;
        var16 = 0;
        jc[0] = 0;
        int currentColumn = 0;

        for(Iterator var15 = map.entrySet().iterator(); var15.hasNext(); ++var16) {
            Map.Entry entry = (Map.Entry)var15.next();
            int var17 = ((Integer)((Pair)entry.getKey()).second).intValue();
            int var18 = ((Integer)((Pair)entry.getKey()).first).intValue();
            pr[var16] = ((Double)entry.getValue()).doubleValue();

            for(ir[var16] = var17; currentColumn < var18; ++currentColumn) {
                jc[currentColumn + 1] = var16;
            }
        }

        while(currentColumn < numColumns) {
            jc[currentColumn + 1] = var16;
            ++currentColumn;
        }

        jc[numColumns] = var16;
        return createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nzmax);
    }

    public SparseMatrix(int[] rIndices, int[] cIndices, double[] values, int numRows, int numColumns) {
        this(rIndices, cIndices, values, numRows, numColumns, values.length);
    }

    public static SparseMatrix createSparseMatrixByCSCArrays(int[] ir, int[] jc, double[] pr, int M, int N, int nzmax) {
        SparseMatrix res = new SparseMatrix();
        res.ir = ir;
        res.jc = jc;
        res.pr = pr;
        res.M = M;
        res.N = N;
        res.nzmax = nzmax;
        int[] ic = new int[pr.length];
        int[] jr = new int[M + 1];
        int[] valCSRIndices = new int[pr.length];
        int[] rIndices = ir;
        int[] cIndices = new int[ic.length];
        int k = 0;
        int j = 0;

        while(k < ir.length && j < N) {
            if(jc[j] <= k && k < jc[j + 1]) {
                cIndices[k] = j;
                ++k;
            } else {
                ++j;
            }
        }

        TreeMap map = new TreeMap();

        for(k = 0; k < pr.length; ++k) {
            if(pr[k] != 0.0D) {
                map.put(Pair.of(Integer.valueOf(rIndices[k]), Integer.valueOf(cIndices[k])), Integer.valueOf(k));
            }
        }

        boolean rIdx = true;
        boolean cIdx = true;
        boolean vIdx = true;
        k = 0;
        jr[0] = 0;
        int currentRow = 0;

        for(Iterator var20 = map.entrySet().iterator(); var20.hasNext(); ++k) {
            Map.Entry entry = (Map.Entry)var20.next();
            int var21 = ((Integer)((Pair)entry.getKey()).first).intValue();
            int var22 = ((Integer)((Pair)entry.getKey()).second).intValue();
            int var23 = ((Integer)entry.getValue()).intValue();
            ic[k] = var22;

            for(valCSRIndices[k] = var23; currentRow < var21; ++currentRow) {
                jr[currentRow + 1] = k;
            }
        }

        while(currentRow < M) {
            jr[currentRow + 1] = k;
            ++currentRow;
        }

        jr[M] = k;
        res.ic = ic;
        res.jr = jr;
        res.valCSRIndices = valCSRIndices;
        res.nnz = map.size();
        return res;
    }

    public static SparseMatrix createSparseMatrixByCSRArrays(int[] ic, int[] jr, double[] pr, int M, int N, int nzmax) {
        int[] rIndices = new int[ic.length];
        int[] cIndices = ic;
        int k = 0;
        int i = 0;

        while(k < ic.length && i < M) {
            if(jr[i] <= k && k < jr[i + 1]) {
                rIndices[k] = i;
                ++k;
            } else {
                ++i;
            }
        }

        TreeMap map = new TreeMap();

        for(k = 0; k < pr.length; ++k) {
            if(pr[k] != 0.0D) {
                map.put(Pair.of(Integer.valueOf(cIndices[k]), Integer.valueOf(rIndices[k])), Double.valueOf(pr[k]));
            }
        }

        int numColumns = N;
        int[] ir = new int[nzmax];
        int[] jc = new int[N + 1];
        pr = new double[nzmax];
        boolean rIdx = true;
        boolean cIdx = true;
        k = 0;
        jc[0] = 0;
        int currentColumn = 0;

        for(Iterator var19 = map.entrySet().iterator(); var19.hasNext(); ++k) {
            Map.Entry entry = (Map.Entry)var19.next();
            int var20 = ((Integer)((Pair)entry.getKey()).second).intValue();
            int var21 = ((Integer)((Pair)entry.getKey()).first).intValue();
            pr[k] = ((Double)entry.getValue()).doubleValue();

            for(ir[k] = var20; currentColumn < var21; ++currentColumn) {
                jc[currentColumn + 1] = k;
            }
        }

        while(currentColumn < numColumns) {
            jc[currentColumn + 1] = k;
            ++currentColumn;
        }

        jc[numColumns] = k;
        return createSparseMatrixByCSCArrays(ir, jc, pr, M, numColumns, nzmax);
    }

    public int[] getIr() {
        return this.ir;
    }

    public int[] getJc() {
        return this.jc;
    }

    public int[] getIc() {
        return this.ic;
    }

    public int[] getJr() {
        return this.jr;
    }

    public double[] getPr() {
        return this.pr;
    }

    public int[] getValCSRIndices() {
        return this.valCSRIndices;
    }

    public int getRowDimension() {
        return this.M;
    }

    public int getColumnDimension() {
        return this.N;
    }

    public int getNZMax() {
        return this.nzmax;
    }

    public int getNNZ() {
        return this.nnz;
    }

    public double[][] getData() {
        double[][] data = new double[this.M][];

        for(int i = 0; i < this.M; ++i) {
            double[] rowData = ArrayOperator.allocateVector(this.N, 0.0D);

            for(int k = this.jr[i]; k < this.jr[i + 1]; ++k) {
                rowData[this.ic[k]] = this.pr[this.valCSRIndices[k]];
            }

            data[i] = rowData;
        }

        return data;
    }

    public Matrix mtimes(Matrix A) {
        Object res = null;
        int NA = A.getColumnDimension();
        boolean rr;
        int kl;
        int kr;
        int var27;
        if(A instanceof DenseMatrix) {
            double[][] ir = new double[this.M][];

            for(int jc = 0; jc < this.M; ++jc) {
                ir[jc] = new double[NA];
            }

            double[] var24 = (double[])null;
            double[][] pr = ((DenseMatrix)A).getData();
            rr = true;
            double cl = 0.0D;

            for(int i = 0; i < this.M; ++i) {
                var24 = ir[i];

                for(kl = 0; kl < NA; ++kl) {
                    cl = 0.0D;

                    for(kr = this.jr[i]; kr < this.jr[i + 1]; ++kr) {
                        var27 = this.ic[kr];
                        cl += this.pr[this.valCSRIndices[kr]] * pr[var27][kl];
                    }

                    var24[kl] = cl;
                }
            }

            res = new DenseMatrix(ir);
        } else if(A instanceof SparseMatrix) {
            int[] var23 = (int[])null;
            int[] var25 = (int[])null;
            double[] var26 = (double[])null;
            var23 = ((SparseMatrix)A).getIr();
            var25 = ((SparseMatrix)A).getJc();
            var26 = ((SparseMatrix)A).getPr();
            rr = true;
            boolean var28 = true;
            double s = 0.0D;
            boolean var30 = false;
            boolean var31 = false;
            int nzmax = 0;
            TreeMap map = new TreeMap();

            int numRows;
            int numColumns;
            for(numRows = 0; numRows < this.M; ++numRows) {
                for(numColumns = 0; numColumns < NA; ++numColumns) {
                    s = 0.0D;
                    kl = this.jr[numRows];
                    kr = var25[numColumns];

                    while(kl < this.jr[numRows + 1] && kr < var25[numColumns + 1]) {
                        int var29 = this.ic[kl];
                        var27 = var23[kr];
                        if(var29 < var27) {
                            ++kl;
                        } else if(var29 > var27) {
                            ++kr;
                        } else {
                            s += this.pr[this.valCSRIndices[kl]] * var26[kr];
                            ++kl;
                            ++kr;
                        }
                    }

                    if(s != 0.0D) {
                        ++nzmax;
                        map.put(Pair.of(Integer.valueOf(numColumns), Integer.valueOf(numRows)), Double.valueOf(s));
                    }
                }
            }

            numRows = this.M;
            numColumns = NA;
            var23 = new int[nzmax];
            var25 = new int[NA + 1];
            var26 = new double[nzmax];
            boolean rIdx = true;
            boolean cIdx = true;
            int k = 0;
            var25[0] = 0;
            int currentColumn = 0;

            for(Iterator var22 = map.entrySet().iterator(); var22.hasNext(); ++k) {
                Map.Entry entry = (Map.Entry)var22.next();
                int var32 = ((Integer)((Pair)entry.getKey()).second).intValue();
                int var33 = ((Integer)((Pair)entry.getKey()).first).intValue();
                var26[k] = ((Double)entry.getValue()).doubleValue();

                for(var23[k] = var32; currentColumn < var33; ++currentColumn) {
                    var25[currentColumn + 1] = k;
                }
            }

            while(currentColumn < numColumns) {
                var25[currentColumn + 1] = k;
                ++currentColumn;
            }

            var25[numColumns] = k;
            res = createSparseMatrixByCSCArrays(var23, var25, var26, numRows, numColumns, nzmax);
        }

        return (Matrix)res;
    }


    public double getEntry(int r, int c) {
        if(r < 0 || r > this.M - 1 || c < 0 || c > this.N - 1) {
            System.err.println("Wrong index.");
            System.exit(1);
        }

        double res = 0.0D;
        boolean idx = true;
        int var10000;
        int u;
        int l;
        int k;
        int idx1;
        if(r <= c) {
            u = this.jc[c + 1] - 1;
            l = this.jc[c];
            if(u < l) {
                return 0.0D;
            }

            var10000 = this.jc[c];

            while(l <= u) {
                k = (u + l) / 2;
                idx1 = this.ir[k];
                if(idx1 == r) {
                    return this.pr[k];
                }

                if(idx1 < r) {
                    l = k + 1;
                } else {
                    u = k - 1;
                }
            }
        } else {
            u = this.jr[r + 1] - 1;
            l = this.jr[r];
            if(u < l) {
                return 0.0D;
            }

            var10000 = this.jr[r];

            while(l <= u) {
                k = (u + l) / 2;
                idx1 = this.ic[k];
                if(idx1 == c) {
                    return this.pr[this.valCSRIndices[k]];
                }

                if(idx1 < c) {
                    l = k + 1;
                } else {
                    u = k - 1;
                }
            }
        }

        return res;
    }

    public void setEntry(int r, int c, double v) {
        if(r < 0 || r > this.M - 1 || c < 0 || c > this.N - 1) {
            System.err.println("Wrong index.");
            System.exit(1);
        }

        int u = this.jc[c + 1] - 1;
        int l = this.jc[c];
        if(u < l) {
            this.insertEntry(r, c, v, this.jc[c]);
        } else {
            boolean idx = true;
            int k = this.jc[c];
            byte flag = 0;

            while(l <= u) {
                k = (u + l) / 2;
                int var10 = this.ir[k];
                if(var10 == r) {
                    if(v == 0.0D) {
                        this.deleteEntry(r, c, k);
                    } else {
                        this.pr[k] = v;
                    }

                    return;
                }

                if(var10 < r) {
                    l = k + 1;
                    flag = 1;
                } else {
                    u = k - 1;
                    flag = 2;
                }
            }

            if(flag == 1) {
                ++k;
            }

            this.insertEntry(r, c, v, k);
        }
    }

    private void insertEntry(int r, int c, double v, int pos) {
        if(v != 0.0D) {
            int len_old = this.pr.length;
            int new_space = len_old < this.M * this.N - 10?10:this.M * this.N - len_old;
            int var14;
            if(this.nnz + 1 > len_old) {
                double[] u = new double[len_old + new_space];
                System.arraycopy(this.pr, 0, u, 0, pos);
                u[pos] = v;
                if(pos < len_old) {
                    System.arraycopy(this.pr, pos, u, pos + 1, len_old - pos);
                }

                this.pr = u;
            } else {
                for(var14 = this.nnz - 1; var14 >= pos; --var14) {
                    this.pr[var14 + 1] = this.pr[var14];
                }

                this.pr[pos] = v;
            }

            if(this.nnz + 1 > len_old) {
                int[] var15 = new int[len_old + new_space];
                System.arraycopy(this.ir, 0, var15, 0, pos);
                var15[pos] = r;
                if(pos < len_old) {
                    System.arraycopy(this.ir, pos, var15, pos + 1, len_old - pos);
                }

                this.ir = var15;
            } else {
                for(var14 = this.nnz - 1; var14 >= pos; --var14) {
                    this.ir[var14 + 1] = this.ir[var14];
                }

                this.ir[pos] = r;
            }

            for(var14 = c + 1; var14 < this.N + 1; ++var14) {
                ++this.jc[var14];
            }

            var14 = this.jr[r + 1] - 1;
            int l = this.jr[r];
            int k = this.jr[r];
            boolean idx = true;
            byte flag = 0;

            while(l <= var14) {
                k = (var14 + l) / 2;
                int var16 = this.ic[k];
                if(var16 != c) {
                    if(var16 < c) {
                        l = k + 1;
                        flag = 1;
                    } else {
                        var14 = k - 1;
                        flag = 2;
                    }
                }
            }

            if(flag == 1) {
                ++k;
            }

            int[] i;
            int var17;
            if(this.nnz + 1 > len_old) {
                i = new int[len_old + new_space];
                System.arraycopy(this.ic, 0, i, 0, k);
                i[k] = c;
                if(k < len_old) {
                    System.arraycopy(this.ic, k, i, k + 1, len_old - k);
                }

                this.ic = i;
            } else {
                for(var17 = this.nnz - 1; var17 >= k; --var17) {
                    this.ic[var17 + 1] = this.ic[var17];
                }

                this.ic[k] = c;
            }

            for(var17 = r + 1; var17 < this.M + 1; ++var17) {
                ++this.jr[var17];
            }

            for(var17 = 0; var17 < this.nnz; ++var17) {
                if(this.valCSRIndices[var17] >= pos) {
                    ++this.valCSRIndices[var17];
                }
            }

            if(this.nnz + 1 > len_old) {
                i = new int[len_old + new_space];
                System.arraycopy(this.valCSRIndices, 0, i, 0, k);
                i[k] = pos;
                if(k < len_old) {
                    System.arraycopy(this.valCSRIndices, k, i, k + 1, len_old - k);
                }

                this.valCSRIndices = i;
            } else {
                for(var17 = this.nnz - 1; var17 >= k; --var17) {
                    this.valCSRIndices[var17 + 1] = this.valCSRIndices[var17];
                }

                this.valCSRIndices[k] = pos;
            }

            ++this.nnz;
            if(this.nnz > len_old) {
                this.nzmax = len_old + new_space;
            }

        }
    }

    private void deleteEntry(int r, int c, int pos) {
        int u;
        for(u = pos; u < this.nnz - 1; ++u) {
            this.pr[u] = this.pr[u + 1];
            this.ir[u] = this.ir[u + 1];
        }

        for(u = c + 1; u < this.N + 1; ++u) {
            --this.jc[u];
        }

        u = this.jr[r + 1] - 1;
        int l = this.jr[r];
        int k = this.jr[r];
        boolean idx = true;

        while(l <= u) {
            k = (u + l) / 2;
            int var9 = this.ic[k];
            if(var9 == c) {
                break;
            }

            if(var9 < c) {
                l = k + 1;
            } else {
                u = k - 1;
            }
        }

        int i;
        for(i = 0; i < this.valCSRIndices.length; ++i) {
            if(this.valCSRIndices[i] > pos) {
                --this.valCSRIndices[i];
            }
        }

        for(i = k; i < this.nnz - 1; ++i) {
            this.ic[i] = this.ic[i + 1];
            this.valCSRIndices[i] = this.valCSRIndices[i + 1];
        }

        for(i = r + 1; i < this.M + 1; ++i) {
            --this.jr[i];
        }

        --this.nnz;
    }

    /** @deprecated */
    @Deprecated
    private Matrix transpose0() {
        double[] values = this.pr;
        int[] rIndices = this.ir;
        int[] cIndices = new int[this.ic.length];
        int k = 0;
        int j = 0;

        while(k < this.nnz && j < this.N) {
            if(this.jc[j] <= k && k < this.jc[j + 1]) {
                cIndices[k] = j;
                ++k;
            } else {
                ++j;
            }
        }

        return createSparseMatrix(cIndices, rIndices, values, this.N, this.M, this.nzmax);
    }

    public Matrix transpose() {
        SparseMatrix res = new SparseMatrix();
        res.M = this.N;
        res.N = this.M;
        res.nnz = this.nnz;
        res.nzmax = this.nzmax;
        res.ir = (int[])this.ic.clone();
        res.jc = (int[])this.jr.clone();
        res.ic = (int[])this.ir.clone();
        res.jr = (int[])this.jc.clone();
        double[] pr_new = new double[this.nzmax];
        boolean k = false;

        int var10;
        for(var10 = 0; var10 < this.nnz; ++var10) {
            pr_new[var10] = this.pr[this.valCSRIndices[var10]];
        }

        res.pr = pr_new;
        int[] valCSRIndices_new = new int[this.nzmax];
        int j = 0;
        boolean rIdx = true;
        boolean cIdx = true;
        var10 = 0;
        boolean k2 = false;
        boolean numBeforeThisEntry = false;

        while(var10 < this.nnz && j < this.N) {
            if(this.jc[j] <= var10 && var10 < this.jc[j + 1]) {
                int var11 = this.ir[var10];
                int var12 = j;
                int var14 = this.jr[var11];

                for(int var13 = this.jr[var11]; var13 < this.jr[var11 + 1] && this.ic[var13] != var12; ++var13) {
                    ++var14;
                }

                valCSRIndices_new[var10] = var14;
                ++var10;
            } else {
                ++j;
            }
        }

        res.valCSRIndices = valCSRIndices_new;
        return res;
    }

    public Matrix plus(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            Object res = null;
            if(A instanceof DenseMatrix) {
                res = A.plus(this);
            } else if(A instanceof SparseMatrix) {
                int[] ir = (int[])null;
                int[] jc = (int[])null;
                double[] pr = (double[])null;
                ir = ((SparseMatrix)A).getIr();
                jc = ((SparseMatrix)A).getJc();
                pr = ((SparseMatrix)A).getPr();
                boolean k1 = false;
                boolean k2 = false;
                boolean r1 = true;
                boolean r2 = true;
                int nzmax = 0;
                boolean i = true;
                double v = 0.0D;
                TreeMap map = new TreeMap();

                int numRows;
                for(numRows = 0; numRows < this.N; ++numRows) {
                    int var23 = this.jc[numRows];
                    int var24 = jc[numRows];
                    if(var23 != this.jc[numRows + 1] || var24 != jc[numRows + 1]) {
                        while(var23 < this.jc[numRows + 1] || var24 < jc[numRows + 1]) {
                            int var27;
                            if(var24 == jc[numRows + 1]) {
                                var27 = this.ir[var23];
                                v = this.pr[var23];
                                ++var23;
                            } else if(var23 == this.jc[numRows + 1]) {
                                var27 = ir[var24];
                                v = pr[var24];
                                ++var24;
                            } else {
                                int var25 = this.ir[var23];
                                int var26 = ir[var24];
                                if(var25 < var26) {
                                    var27 = var25;
                                    v = this.pr[var23];
                                    ++var23;
                                } else if(var25 == var26) {
                                    var27 = var25;
                                    v = this.pr[var23] + pr[var24];
                                    ++var23;
                                    ++var24;
                                } else {
                                    var27 = var26;
                                    v = pr[var24];
                                    ++var24;
                                }
                            }

                            if(v != 0.0D) {
                                map.put(Pair.of(Integer.valueOf(numRows), Integer.valueOf(var27)), Double.valueOf(v));
                                ++nzmax;
                            }
                        }
                    }
                }

                numRows = this.M;
                int numColumns = this.N;
                ir = new int[nzmax];
                jc = new int[numColumns + 1];
                pr = new double[nzmax];
                boolean rIdx = true;
                boolean cIdx = true;
                int k = 0;
                jc[0] = 0;
                int currentColumn = 0;

                for(Iterator var22 = map.entrySet().iterator(); var22.hasNext(); ++k) {
                    Map.Entry entry = (Map.Entry)var22.next();
                    int var28 = ((Integer)((Pair)entry.getKey()).second).intValue();
                    int var29 = ((Integer)((Pair)entry.getKey()).first).intValue();
                    pr[k] = ((Double)entry.getValue()).doubleValue();

                    for(ir[k] = var28; currentColumn < var29; ++currentColumn) {
                        jc[currentColumn + 1] = k;
                    }
                }

                while(currentColumn < numColumns) {
                    jc[currentColumn + 1] = k;
                    ++currentColumn;
                }

                jc[numColumns] = k;
                res = createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nzmax);
            }

            return (Matrix)res;
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix minus(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            Object res = null;
            int k1;
            int k2;
            if(A instanceof DenseMatrix) {
                res = A.copy();
                double[][] ir = ((DenseMatrix)res).getData();
                boolean jc = true;
                boolean pr = false;

                for(k1 = 0; k1 < this.N; ++k1) {
                    for(k2 = 0; k2 < this.M; ++k2) {
                        ir[k2][k1] = -ir[k2][k1];
                    }

                    for(int var26 = this.jc[k1]; var26 < this.jc[k1 + 1]; ++var26) {
                        int var24 = this.ir[var26];
                        ir[var24][k1] += this.pr[var26];
                    }
                }
            } else if(A instanceof SparseMatrix) {
                int[] var23 = (int[])null;
                int[] var25 = (int[])null;
                double[] var27 = (double[])null;
                var23 = ((SparseMatrix)A).getIr();
                var25 = ((SparseMatrix)A).getJc();
                var27 = ((SparseMatrix)A).getPr();
                boolean var28 = false;
                boolean var29 = false;
                boolean r1 = true;
                boolean r2 = true;
                int nzmax = 0;
                boolean i = true;
                double v = 0.0D;
                TreeMap map = new TreeMap();

                int numRows;
                for(numRows = 0; numRows < this.N; ++numRows) {
                    k1 = this.jc[numRows];
                    k2 = var25[numRows];
                    if(k1 != this.jc[numRows + 1] || k2 != var25[numRows + 1]) {
                        while(k1 < this.jc[numRows + 1] || k2 < var25[numRows + 1]) {
                            int var32;
                            if(k2 == var25[numRows + 1]) {
                                var32 = this.ir[k1];
                                v = this.pr[k1];
                                ++k1;
                            } else if(k1 == this.jc[numRows + 1]) {
                                var32 = var23[k2];
                                v = -var27[k2];
                                ++k2;
                            } else {
                                int var30 = this.ir[k1];
                                int var31 = var23[k2];
                                if(var30 < var31) {
                                    var32 = var30;
                                    v = this.pr[k1];
                                    ++k1;
                                } else if(var30 == var31) {
                                    var32 = var30;
                                    v = this.pr[k1] - var27[k2];
                                    ++k1;
                                    ++k2;
                                } else {
                                    var32 = var31;
                                    v = -var27[k2];
                                    ++k2;
                                }
                            }

                            if(v != 0.0D) {
                                map.put(Pair.of(Integer.valueOf(numRows), Integer.valueOf(var32)), Double.valueOf(v));
                                ++nzmax;
                            }
                        }
                    }
                }

                numRows = this.M;
                int numColumns = this.N;
                var23 = new int[nzmax];
                var25 = new int[numColumns + 1];
                var27 = new double[nzmax];
                boolean rIdx = true;
                boolean cIdx = true;
                int k = 0;
                var25[0] = 0;
                int currentColumn = 0;

                for(Iterator var22 = map.entrySet().iterator(); var22.hasNext(); ++k) {
                    Map.Entry entry = (Map.Entry)var22.next();
                    int var33 = ((Integer)((Pair)entry.getKey()).second).intValue();
                    int var34 = ((Integer)((Pair)entry.getKey()).first).intValue();
                    var27[k] = ((Double)entry.getValue()).doubleValue();

                    for(var23[k] = var33; currentColumn < var34; ++currentColumn) {
                        var25[currentColumn + 1] = k;
                    }
                }

                while(currentColumn < numColumns) {
                    var25[currentColumn + 1] = k;
                    ++currentColumn;
                }

                var25[numColumns] = k;
                res = createSparseMatrixByCSCArrays(var23, var25, var27, numRows, numColumns, nzmax);
            }

            return (Matrix)res;
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix times(Matrix A) {
        if(A.getRowDimension() == this.M && A.getColumnDimension() == this.N) {
            Object res = null;
            if(A instanceof DenseMatrix) {
                res = A.times(this);
            } else if(A instanceof SparseMatrix) {
                int[] ir = (int[])null;
                int[] jc = (int[])null;
                double[] pr = (double[])null;
                ir = ((SparseMatrix)A).getIr();
                jc = ((SparseMatrix)A).getJc();
                pr = ((SparseMatrix)A).getPr();
                boolean k1 = false;
                boolean k2 = false;
                boolean r1 = true;
                boolean r2 = true;
                int nzmax = 0;
                boolean i = true;
                double v = 0.0D;
                TreeMap map = new TreeMap();

                int numRows;
                for(numRows = 0; numRows < this.N; ++numRows) {
                    int var23 = this.jc[numRows];
                    int var24 = jc[numRows];
                    if(var23 != this.jc[numRows + 1] && var24 != jc[numRows + 1]) {
                        while(var23 < this.jc[numRows + 1] && var24 < jc[numRows + 1]) {
                            int var25 = this.ir[var23];
                            int var26 = ir[var24];
                            if(var25 < var26) {
                                ++var23;
                            } else if(var25 == var26) {
                                v = this.pr[var23] * pr[var24];
                                ++var23;
                                ++var24;
                                if(v != 0.0D) {
                                    map.put(Pair.of(Integer.valueOf(numRows), Integer.valueOf(var25)), Double.valueOf(v));
                                    ++nzmax;
                                }
                            } else {
                                ++var24;
                            }
                        }
                    }
                }

                numRows = this.M;
                int numColumns = this.N;
                ir = new int[nzmax];
                jc = new int[numColumns + 1];
                pr = new double[nzmax];
                boolean rIdx = true;
                boolean cIdx = true;
                int k = 0;
                jc[0] = 0;
                int currentColumn = 0;

                for(Iterator var22 = map.entrySet().iterator(); var22.hasNext(); ++k) {
                    Map.Entry entry = (Map.Entry)var22.next();
                    int var27 = ((Integer)((Pair)entry.getKey()).second).intValue();
                    int var28 = ((Integer)((Pair)entry.getKey()).first).intValue();
                    pr[k] = ((Double)entry.getValue()).doubleValue();

                    for(ir[k] = var27; currentColumn < var28; ++currentColumn) {
                        jc[currentColumn + 1] = k;
                    }
                }

                while(currentColumn < numColumns) {
                    jc[currentColumn + 1] = k;
                    ++currentColumn;
                }

                jc[numColumns] = k;
                res = createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nzmax);
            }

            return (Matrix)res;
        } else {
            System.err.println("Dimension doesn\'t match.");
            return null;
        }
    }

    public Matrix times(double v) {
        if(v == 0.0D) {
            return new SparseMatrix(this.M, this.N);
        } else {
            SparseMatrix res = (SparseMatrix)this.copy();

            for(int k = 0; k < this.nnz; ++k) {
                res.pr[k] = v * this.pr[k];
            }

            return res;
        }
    }

    public Matrix copy() {
        SparseMatrix res = new SparseMatrix();
        res.ir = (int[])this.ir.clone();
        res.jc = (int[])this.jc.clone();
        res.pr = (double[])this.pr.clone();
        res.ic = (int[])this.ic.clone();
        res.jr = (int[])this.jr.clone();
        res.valCSRIndices = (int[])this.valCSRIndices.clone();
        res.M = this.M;
        res.N = this.N;
        res.nzmax = this.nzmax;
        res.nnz = this.nnz;
        return res;
    }

    public Matrix clone() {
        return this.copy();
    }

    public Matrix plus(double v) {
        DenseMatrix res = new DenseMatrix(this.M, this.N, v);
        double[][] data = ((DenseMatrix)res).getData();

        for(int j = 0; j < this.N; ++j) {
            for(int k = this.jc[j]; k < this.jc[j + 1]; ++k) {
                data[this.ir[k]][j] += this.pr[k];
            }
        }

        return res;
    }

    public Matrix minus(double v) {
        return this.plus(-v);
    }

    public Vector operate(Vector b) {
        if(this.N != b.getDim()) {
            System.err.println("Dimension does not match.");
            System.exit(1);
        }

        Object res = null;
        double[] pr;
        int kl;
        int kr;
        if(b instanceof DenseVector) {
            double[] ir = new double[this.M];
            pr = ((DenseVector)b).getPr();
            double nnz = 0.0D;
            boolean c = false;

            for(kl = 0; kl < this.M; ++kl) {
                nnz = 0.0D;

                for(kr = this.jr[kl]; kr < this.jr[kl + 1]; ++kr) {
                    int var18 = this.ic[kr];
                    nnz += this.pr[this.valCSRIndices[kr]] * pr[var18];
                }

                ir[kl] = nnz;
            }

            res = new DenseVector(ir);
        } else if(b instanceof SparseVector) {
            int[] var16 = ((SparseVector)b).getIr();
            pr = ((SparseVector)b).getPr();
            int var17 = ((SparseVector)b).getNNZ();
            double s = 0.0D;
            boolean var19 = false;
            boolean var20 = false;
            boolean cl = false;
            boolean rr = false;
            TreeMap map = new TreeMap();

            int ind;
            for(ind = 0; ind < this.M; ++ind) {
                kl = this.jr[ind];
                kr = 0;
                s = 0.0D;

                while(kl < this.jr[ind + 1] && kr < var17) {
                    int var21 = this.ic[kl];
                    int var22 = var16[kr];
                    if(var21 < var22) {
                        ++kl;
                    } else if(var21 > var22) {
                        ++kr;
                    } else {
                        s += this.pr[this.valCSRIndices[kl]] * pr[kr];
                        ++kl;
                        ++kr;
                    }
                }

                if(s != 0.0D) {
                    map.put(Integer.valueOf(ind), Double.valueOf(s));
                }
            }

            var17 = map.size();
            var16 = new int[var17];
            pr = new double[var17];
            ind = 0;

            for(Iterator var15 = map.entrySet().iterator(); var15.hasNext(); ++ind) {
                Map.Entry entry = (Map.Entry)var15.next();
                var16[ind] = ((Integer)entry.getKey()).intValue();
                pr[ind] = ((Double)entry.getValue()).doubleValue();
            }

            res = new SparseVector(var16, pr, var17, this.M);
        }

        return (Vector)res;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        if(this.getNNZ() == 0) {
            sb.append("Empty sparse matrix." + System.lineSeparator());
            return sb.toString();
        } else {
            byte p = 4;
            int[] ic = this.getIc();
            int[] jr = this.getJr();
            double[] pr = this.getPr();
            int[] valCSRIndices = this.getValCSRIndices();
            int M = this.getRowDimension();
            String valueString = "";

            for(int r = 0; r < M; ++r) {
                sb.append("  ");
                boolean currentColumn = false;
                int lastColumn = -1;

                for(int k = jr[r]; k < jr[r + 1]; ++k) {
                    int var17;
                    for(var17 = ic[k]; lastColumn < var17 - 1; ++lastColumn) {
                        sb.append(Printer.sprintf(String.format("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{" "}));
                        sb.append("  ");
                    }

                    lastColumn = var17;
                    double v = pr[valCSRIndices[k]];
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
        this.nzmax = 0;
        this.nnz = 0;
        this.jc = new int[this.N + 1];

        int i;
        for(i = 0; i < this.N + 1; ++i) {
            this.jc[i] = 0;
        }

        this.jr = new int[this.M + 1];

        for(i = 0; i < this.M + 1; ++i) {
            this.jr[i] = 0;
        }

        this.ir = new int[0];
        this.pr = new double[0];
        this.ic = new int[0];
        this.valCSRIndices = new int[0];
    }

    public void clean() {
        TreeSet set = new TreeSet();

        for(int pair = 0; pair < this.nnz; ++pair) {
            if(this.pr[pair] == 0.0D) {
                set.add(Pair.of(Integer.valueOf(this.ir[pair]), Integer.valueOf(this.ic[pair])));
            }
        }

        Iterator var3 = set.iterator();

        while(var3.hasNext()) {
            Pair var4 = (Pair)var3.next();
            this.setEntry(((Integer)var4.first).intValue(), ((Integer)var4.second).intValue(), 0.0D);
        }

    }

    public void appendAnEmptyRow() {
        int[] jr = new int[this.M + 2];
        System.arraycopy(this.jr, 0, jr, 0, this.M + 1);
        jr[this.M + 1] = this.jr[this.M];
        ++this.M;
        this.jr = jr;
    }

    public void appendAnEmptyColumn() {
        int[] jc = new int[this.N + 2];
        System.arraycopy(this.jc, 0, jc, 0, this.N + 1);
        jc[this.N + 1] = this.jc[this.N];
        ++this.N;
        this.jc = jc;
    }

    public Matrix getSubMatrix(int startRow, int endRow, int startColumn, int endColumn) {
        int nnz = 0;
        int numRows = endRow - startRow + 1;
        int numColumns = endColumn - startColumn + 1;
        boolean rowIdx = true;

        int var19;
        for(int nzmax = startColumn; nzmax <= endColumn; ++nzmax) {
            for(int ir = this.jc[nzmax]; ir < this.jc[nzmax + 1]; ++ir) {
                var19 = this.ir[ir];
                if(var19 >= startRow) {
                    if(var19 > endRow) {
                        break;
                    }

                    ++nnz;
                }
            }
        }

        int[] var20 = new int[nnz];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        boolean cIdx = true;
        int k = 0;
        jc[0] = 0;
        int currentColumn = startColumn;

        for(int j = startColumn; j <= endColumn; ++j) {
            for(int t = this.jc[j]; t < this.jc[j + 1]; ++t) {
                var19 = this.ir[t];
                if(var19 >= startRow) {
                    if(var19 > endRow) {
                        break;
                    }

                    int var21 = var19 - startRow;
                    int var22 = j;
                    pr[k] = this.pr[t];

                    for(var20[k] = var21; currentColumn < var22; ++currentColumn) {
                        jc[currentColumn + 1 - startColumn] = k;
                    }

                    ++k;
                }
            }
        }

        while(currentColumn < numColumns) {
            jc[currentColumn + 1 - startColumn] = k;
            ++currentColumn;
        }

        jc[numColumns] = k;
        return createSparseMatrixByCSCArrays(var20, jc, pr, numRows, numColumns, nnz);
    }

    public Matrix getSubMatrix(int[] selectedRows, int[] selectedColumns) {
        int nRow = selectedRows.length;
        int nCol = selectedColumns.length;
        double v = 0.0D;
        SparseMatrix res = new SparseMatrix(nRow, nCol);

        for(int j = 0; j < nCol; ++j) {
            int c = selectedColumns[j];

            for(int i = 0; i < nRow; ++i) {
                int r = selectedRows[i];
                v = this.getEntry(r, c);
                if(v != 0.0D) {
                    res.setEntry(i, j, v);
                }
            }
        }

        return res;
    }

    public Matrix getColumnMatrix(int c) {
        int nnz = this.jc[c + 1] - this.jc[c];
        if(nnz == 0) {
            return new SparseMatrix(this.M, 1);
        } else {
            int[] ir = new int[nnz];
            int[] jc = new int[]{0, nnz};
            double[] pr = new double[nnz];
            int k = this.jc[c];

            for(int i = 0; k < this.jc[c + 1]; ++i) {
                ir[i] = this.ir[k];
                pr[i] = this.pr[k];
                ++k;
            }

            return createSparseMatrixByCSCArrays(ir, jc, pr, this.M, 1, nnz);
        }
    }

    public Vector getColumnVector(int c) {
        int dim = this.M;
        int nnz = this.jc[c + 1] - this.jc[c];
        if(nnz == 0) {
            return new SparseVector(dim);
        } else {
            int[] ir = new int[nnz];
            double[] pr = new double[nnz];
            int k = this.jc[c];

            for(int i = 0; k < this.jc[c + 1]; ++i) {
                ir[i] = this.ir[k];
                pr[i] = this.pr[k];
                ++k;
            }

            return new SparseVector(ir, pr, nnz, dim);
        }
    }

    public Matrix getRowMatrix(int r) {
        int nnz = this.jr[r + 1] - this.jr[r];
        if(nnz == 0) {
            return new SparseMatrix(1, this.N);
        } else {
            int[] ic = new int[nnz];
            int[] jr = new int[]{0, nnz};
            double[] pr = new double[nnz];
            int k = this.jr[r];

            for(int j = 0; k < this.jr[r + 1]; ++j) {
                ic[j] = this.ic[k];
                pr[j] = this.pr[this.valCSRIndices[k]];
                ++k;
            }

            return createSparseMatrixByCSRArrays(ic, jr, pr, 1, this.N, nnz);
        }
    }

    public Vector getRowVector(int r) {
        int dim = this.N;
        int nnz = this.jr[r + 1] - this.jr[r];
        if(nnz == 0) {
            return new SparseVector(dim);
        } else {
            int[] ir = new int[nnz];
            double[] pr = new double[nnz];
            int k = this.jr[r];

            for(int j = 0; k < this.jr[r + 1]; ++j) {
                ir[j] = this.ic[k];
                pr[j] = this.pr[this.valCSRIndices[k]];
                ++k;
            }

            return new SparseVector(ir, pr, nnz, dim);
        }
    }

    public Matrix getRows(int startRow, int endRow) {
        int nnz = 0;
        int numRows = endRow - startRow + 1;
        int numColumns = this.N;
        boolean rowIdx = true;

        int var17;
        for(int nzmax = 0; nzmax < numColumns; ++nzmax) {
            for(int ir = this.jc[nzmax]; ir < this.jc[nzmax + 1]; ++ir) {
                var17 = this.ir[ir];
                if(var17 >= startRow) {
                    if(var17 > endRow) {
                        break;
                    }

                    ++nnz;
                }
            }
        }

        int[] var18 = new int[nnz];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        boolean cIdx = true;
        int k = 0;
        jc[0] = 0;
        int currentColumn = 0;

        for(int j = 0; j < numColumns; ++j) {
            for(int t = this.jc[j]; t < this.jc[j + 1]; ++t) {
                var17 = this.ir[t];
                if(var17 >= startRow) {
                    if(var17 > endRow) {
                        break;
                    }

                    int var19 = var17 - startRow;
                    int var20 = j;
                    pr[k] = this.pr[t];

                    for(var18[k] = var19; currentColumn < var20; ++currentColumn) {
                        jc[currentColumn + 1] = k;
                    }

                    ++k;
                }
            }
        }

        while(currentColumn < numColumns) {
            jc[currentColumn + 1] = k;
            ++currentColumn;
        }

        jc[numColumns] = k;
        return createSparseMatrixByCSCArrays(var18, jc, pr, numRows, numColumns, nnz);
    }

    public Matrix getRows(int... selectedRows) {
        int nRow = selectedRows.length;
        int nCol = this.N;
        double v = 0.0D;
        SparseMatrix res = new SparseMatrix(nRow, nCol);

        for(int j = 0; j < nCol; ++j) {
            int c = j;

            for(int i = 0; i < nRow; ++i) {
                int r = selectedRows[i];
                v = this.getEntry(r, c);
                if(v != 0.0D) {
                    res.setEntry(i, j, v);
                }
            }
        }

        return res;
    }

    public Vector[] getRowVectors(int startRow, int endRow) {
        int numRows = endRow - startRow + 1;
        Vector[] res = new Vector[numRows];
        int r = startRow;

        for(int i = 0; r <= endRow; ++i) {
            res[i] = this.getRowVector(r);
            ++r;
        }

        return res;
    }

    public Vector[] getRowVectors(int... selectedRows) {
        int numRows = selectedRows.length;
        Vector[] res = new Vector[numRows];

        for(int i = 0; i < numRows; ++i) {
            res[i] = this.getRowVector(selectedRows[i]);
        }

        return res;
    }

    public Matrix getColumns(int startColumn, int endColumn) {
        int nnz = 0;
        int numRows = this.M;
        int numColumns = endColumn - startColumn + 1;
        boolean rowIdx = true;

        for(int nzmax = startColumn; nzmax <= endColumn; ++nzmax) {
            nnz += this.jc[nzmax + 1] - this.jc[nzmax];
        }

        int[] ir = new int[nnz];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        boolean cIdx = true;
        int k = 0;
        jc[0] = 0;
        int currentColumn = startColumn;

        for(int j = startColumn; j <= endColumn; ++j) {
            for(int t = this.jc[j]; t < this.jc[j + 1]; ++t) {
                int var17 = this.ir[t];
                int var18 = j;
                pr[k] = this.pr[t];

                for(ir[k] = var17; currentColumn < var18; ++currentColumn) {
                    jc[currentColumn + 1 - startColumn] = k;
                }

                ++k;
            }
        }

        while(currentColumn <= endColumn) {
            jc[currentColumn + 1 - startColumn] = k;
            ++currentColumn;
        }

        return createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nnz);
    }

    public Matrix getColumns(int... selectedColumns) {
        int nnz = 0;
        int numRows = this.M;
        int numColumns = selectedColumns.length;
        boolean rowIdx = true;
        boolean j = true;

        int var16;
        for(int nzmax = 0; nzmax < numColumns; ++nzmax) {
            var16 = selectedColumns[nzmax];
            nnz += this.jc[var16 + 1] - this.jc[var16];
        }

        int[] ir = new int[nnz];
        int[] jc = new int[numColumns + 1];
        double[] pr = new double[nnz];
        boolean rIdx = true;
        int k = 0;
        jc[0] = 0;

        for(int c = 0; c < numColumns; ++c) {
            var16 = selectedColumns[c];
            jc[c + 1] = jc[c] + this.jc[var16 + 1] - this.jc[var16];

            for(int t = this.jc[var16]; t < this.jc[var16 + 1]; ++t) {
                int var15 = this.ir[t];
                pr[k] = this.pr[t];
                ir[k] = var15;
                ++k;
            }
        }

        return createSparseMatrixByCSCArrays(ir, jc, pr, numRows, numColumns, nnz);
    }

    public Vector[] getColumnVectors(int startColumn, int endColumn) {
        int numColumns = endColumn - startColumn + 1;
        Vector[] res = new Vector[numColumns];
        int c = startColumn;

        for(int i = 0; c <= endColumn; ++i) {
            res[i] = this.getColumnVector(c);
            ++c;
        }

        return res;
    }

    public Vector[] getColumnVectors(int... selectedColumns) {
        int numColumns = selectedColumns.length;
        Vector[] res = new Vector[numColumns];

        for(int j = 0; j < numColumns; ++j) {
            res[j] = this.getColumnVector(selectedColumns[j]);
        }

        return res;
    }

    public void setRowMatrix(int r, Matrix A) {
        if(A.getRowDimension() != 1) {
            Printer.err("Input matrix should be a row matrix.");
//            Utility.exit(1);
        }

        for(int j = 0; j < this.N; ++j) {
            this.setEntry(r, j, A.getEntry(0, j));
        }

    }

    public void setRowVector(int r, Vector V) {
        for(int j = 0; j < this.N; ++j) {
            this.setEntry(r, j, V.get(j));
        }

    }

    public void setColumnMatrix(int c, Matrix A) {
        if(A.getColumnDimension() != 1) {
            Printer.err("Input matrix should be a column matrix.");
//            Utility.exit(1);
        }

        for(int i = 0; i < this.M; ++i) {
            this.setEntry(i, c, A.getEntry(i, 0));
        }

    }

    public void setColumnVector(int c, Vector V) {
        for(int i = 0; i < this.M; ++i) {
            this.setEntry(i, c, V.get(i));
        }

    }

    public double norm1(){
        double cur;
        double max = 0;
        Vector curvector;
        for (int i = 0; i < this.N; i++) {
            cur = 0;
            curvector = this.getColumnVector(i);
            for (int j = 0; j < this.M; j++) {
                cur+=Math.abs(curvector.get(j));
            }
            max=(cur>max)?cur:max;
        }

        return max;
    }

    public double normInf(){
        double cur;
        double max = 0;
        Vector curvector;
        for (int i = 0; i < this.M; i++) {
            cur = 0;
            curvector = this.getRowVector(i);
            for (int j = 0; j < this.N; j++) {
                cur+=Math.abs(curvector.get(j));
            }
            max=(cur>max)?cur:max;
        }
        return max;
    }

    public double cond1(){
        return this.norm1()* Matlab.inv(this).norm1();
    }

    public double condInf(){
        return this.normInf()* Matlab.inv(this).normInf();
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
            QR= QRDecomposition.decompose(AtA, false);
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

    public double normE(){
        double temp=0;
        for (int i = 0; i < this.getRowDimension(); i++) {
            for (int j = 0; j < this.getColumnDimension(); j++) {
                temp+=this.getEntry(i,j)*this.getEntry(i,j);
            }
        }
        return Math.sqrt(temp);
    }

    public double condE(){
        return this.normE()*Matlab.inv(this).normE();
    }

    public double cond2(){
        return this.norm2()*Matlab.inv(this).norm2();
    }

}
