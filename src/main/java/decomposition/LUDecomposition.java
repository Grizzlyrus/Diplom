package decomposition;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.ArrayOperator;
import utils.Matlab;
import utils.Printer;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

/**
 * Created by Admin on 08.03.2016.
 */
public class LUDecomposition {
    private Matrix L;
    private Matrix U;
    private Matrix P;
    private int numRowExchange;
    private int Lfilling=0;
    private int Ufilling=0;

    public Matrix getL() {
        return this.L;
    }

    public Matrix getU() {
        return this.U;
    }

    public Matrix getP() {
        return this.P;
    }

    public LUDecomposition(Matrix A) {
        Matrix[] LUP = this.run(A);
        this.L = LUP[0];
        this.U = LUP[1];
        this.P = LUP[2];
    }

    private Matrix[] run(Matrix A) {
        int n = A.getRowDimension();
        if(n != A.getColumnDimension()) {
            System.err.println("A should be a square matrix.");
            System.exit(1);
        }

        this.numRowExchange = 0;
        Matrix[] LUP = new Matrix[3];
        int i;
        double maxVal;
        int j;
        int A_i;
        double L_ki;
        int k;
        int l;
        if(A instanceof DenseMatrix) {
            double[][] LVs = (new DenseMatrix(n, n, 0.0D)).getData();
            double[][] AVs = ((DenseMatrix)A.copy()).getData();
            double[][] PVs = (new DenseMatrix(n, n, 0.0D)).getData();

            for(i = 0; i < n; ++i) {
                PVs[i][i] = 1.0D;
            }

            for(i = 0; i < n; ++i) {
                maxVal = AVs[i][i];
                j = i;

                for(A_i = i + 1; A_i < n; ++A_i) {
                    if(Math.abs(maxVal) < Math.abs(AVs[A_i][i])) {
                        j = A_i;
                        maxVal = AVs[A_i][i];
                    }
                }

//                TODO
                if(maxVal == 0.0D) {
                    System.err.println("Matrix A is singular.");
                    LUP[0] = null;
                    LUP[1] = null;
                    LUP[2] = null;
                    return LUP;
                }

                double[] var20;
                if(j != i) {
                    var20 = (double[])null;
                    var20 = AVs[i];
                    AVs[i] = AVs[j];
                    AVs[j] = var20;
                    var20 = LVs[i];
                    LVs[i] = LVs[j];
                    LVs[j] = var20;
                    var20 = PVs[i];
                    PVs[i] = PVs[j];
                    PVs[j] = var20;
                    ++this.numRowExchange;
                }

                LVs[i][i] = 1.0D;
                var20 = AVs[i];
                L_ki = 0.0D;

                for(k = i + 1; k < n; ++k) {
                    L_ki = AVs[k][i] / maxVal;
                    LVs[k][i] = L_ki;
                    double[] A_k = AVs[k];
                    A_k[i] = 0.0D;

                    for(l = i + 1; l < n; ++l) {
                        A_k[l] -= L_ki * var20[l];
                    }
                }
            }

            LUP[0] = new DenseMatrix(LVs);
            LUP[1] = new DenseMatrix(AVs);
            LUP[2] = new DenseMatrix(PVs);
        } else if(A instanceof SparseMatrix) {
            Vector[] var17 = Matlab.sparseMatrix2SparseRowVectors(new SparseMatrix(n, n));
            Vector[] var18 = Matlab.sparseMatrix2SparseRowVectors(A);
            Vector[] var19 = Matlab.sparseMatrix2SparseRowVectors(new SparseMatrix(n, n));

            for(i = 0; i < n; ++i) {
                var19[i].set(i, 1.0D);
            }

            for(i = 0; i < n; ++i) {
                maxVal = var18[i].get(i);
                j = i;

                for(A_i = i + 1; A_i < n; ++A_i) {
                    L_ki = var18[A_i].get(i);
                    if(Math.abs(maxVal) < Math.abs(L_ki)) {
                        j = A_i;
                        maxVal = L_ki;
                    }
                }

                if(maxVal == 0.0D) {
                    System.err.println("Matrix A is singular.");
                    LUP[0] = null;
                    LUP[1] = null;
                    LUP[2] = null;
                    return LUP;
                }

                Vector var21;
                if(j != i) {
                    var21 = null;
                    var21 = var18[i];
                    var18[i] = var18[j];
                    var18[j] = var21;
                    var21 = var17[i];
                    var17[i] = var17[j];
                    var17[j] = var21;
                    var21 = var19[i];
                    var19[i] = var19[j];
                    var19[j] = var21;
                    ++this.numRowExchange;
                }

                var17[i].set(i, 1.0D);
                var21 = var18[i];
                L_ki = 0.0D;

                for(k = i + 1; k < n; ++k) {
                    L_ki = var18[k].get(i) / maxVal;
                    var17[k].set(i, L_ki);
                    Vector var22 = var18[k];
                    var22.set(i, 0.0D);

                    for(l = i + 1; l < n; ++l) {
                        var22.set(l, var22.get(l) - L_ki * var21.get(l));
                    }
                }
            }

            LUP[0] = Matlab.sparseRowVectors2SparseMatrix(var17);
            LUP[1] = Matlab.sparseRowVectors2SparseMatrix(var18);
            LUP[2] = Matlab.sparseRowVectors2SparseMatrix(var19);
        }
        for (int q = 0; q < A.getRowDimension(); q++) {
            for (int s = 0; s < A.getColumnDimension(); s++) {
                if(LUP[0].getEntry(q,s)!=0.0D && A.getEntry(q,s)==0.0D){
                    Lfilling++;
                }
                if(LUP[1].getEntry(q,s)!=0.0D && A.getEntry(q,s)==0.0D){
                    Ufilling++;
                }
            }
        }
        return LUP;
    }

    public static Matrix[] decompose(Matrix A) {
        int n = A.getRowDimension();
        if(n != A.getColumnDimension()) {
            System.err.println("A should be a square matrix.");
            System.exit(1);
        }

        Matrix[] LUP = new Matrix[3];
        int i;
        double maxVal;
        int j;
        int A_i;
        double L_ki;
        int k;
        int l;
        if(A instanceof DenseMatrix) {
            double[][] LVs = (new DenseMatrix(n, n, 0.0D)).getData();
            double[][] AVs = ((DenseMatrix)A.copy()).getData();
            double[][] PVs = (new DenseMatrix(n, n, 0.0D)).getData();

            for(i = 0; i < n; ++i) {
                PVs[i][i] = 1.0D;
            }

            for(i = 0; i < n; ++i) {
                maxVal = AVs[i][i];
                j = i;

                for(A_i = i + 1; A_i < n; ++A_i) {
                    if(maxVal < AVs[A_i][i]) {
                        j = A_i;
                        maxVal = AVs[A_i][i];
                    }
                }

                if(maxVal == 0.0D) {
                    System.err.println("Matrix A is singular.");
                    LUP[0] = null;
                    LUP[1] = null;
                    LUP[2] = null;
                    return LUP;
                }

                double[] var19;
                if(j != i) {
                    var19 = (double[])null;
                    var19 = AVs[i];
                    AVs[i] = AVs[j];
                    AVs[j] = var19;
                    var19 = LVs[i];
                    LVs[i] = LVs[j];
                    LVs[j] = var19;
                    var19 = PVs[i];
                    PVs[i] = PVs[j];
                    PVs[j] = var19;
                }

                LVs[i][i] = 1.0D;
                var19 = AVs[i];
                L_ki = 0.0D;

                for(k = i + 1; k < n; ++k) {
                    L_ki = AVs[k][i] / maxVal;
                    LVs[k][i] = L_ki;
                    double[] A_k = AVs[k];
                    A_k[i] = 0.0D;

                    for(l = i + 1; l < n; ++l) {
                        A_k[l] -= L_ki * var19[l];
                    }
                }
            }

            LUP[0] = new DenseMatrix(LVs);
            LUP[1] = new DenseMatrix(AVs);
            LUP[2] = new DenseMatrix(PVs);
        } else if(A instanceof SparseMatrix) {
            Vector[] var16 = Matlab.sparseMatrix2SparseRowVectors(new SparseMatrix(n, n));
            Vector[] var17 = Matlab.sparseMatrix2SparseRowVectors(A);
            Vector[] var18 = Matlab.sparseMatrix2SparseRowVectors(new SparseMatrix(n, n));

            for(i = 0; i < n; ++i) {
                var18[i].set(i, 1.0D);
            }

            for(i = 0; i < n; ++i) {
                maxVal = var17[i].get(i);
                j = i;

                for(A_i = i + 1; A_i < n; ++A_i) {
                    L_ki = var17[A_i].get(i);
                    if(maxVal < L_ki) {
                        j = A_i;
                        maxVal = L_ki;
                    }
                }

                if(maxVal == 0.0D) {
                    System.err.println("Matrix A is singular.");
                    LUP[0] = null;
                    LUP[1] = null;
                    LUP[2] = null;
                    return LUP;
                }

                Vector var20;
                if(j != i) {
                    var20 = null;
                    var20 = var17[i];
                    var17[i] = var17[j];
                    var17[j] = var20;
                    var20 = var16[i];
                    var16[i] = var16[j];
                    var16[j] = var20;
                    var20 = var18[i];
                    var18[i] = var18[j];
                    var18[j] = var20;
                }

                var16[i].set(i, 1.0D);
                var20 = var17[i];
                L_ki = 0.0D;

                for(k = i + 1; k < n; ++k) {
                    L_ki = var17[k].get(i) / maxVal;
                    var16[k].set(i, L_ki);
                    Vector var21 = var17[k];
                    var21.set(i, 0.0D);

                    for(l = i + 1; l < n; ++l) {
                        var21.set(l, var21.get(l) - L_ki * var20.get(l));
                    }
                }
            }

            LUP[0] = Matlab.sparseRowVectors2SparseMatrix(var16);
            LUP[1] = Matlab.sparseRowVectors2SparseMatrix(var17);
            LUP[2] = Matlab.sparseRowVectors2SparseMatrix(var18);
        }

        return LUP;
    }

    private static void swap(double[] V1, double[] V2, int start, int end) {
        double temp = 0.0D;

        for(int i = start; i < end; ++i) {
            temp = V1[i];
            V1[i] = V2[i];
            V2[i] = temp;
        }

    }

    public Vector solve(double[] b) {
        return this.solve((Vector)(new DenseVector(b)));
    }

    public Vector solve(Vector b) {
        DenseVector res = null;
        double[] d;
        int n;
        double[] y;
        double v;
        int UVs;
        double[] x;
        int i;
        int ir;
        if(this.L instanceof DenseMatrix) {
            d = Matlab.full(this.P.operate(b)).getPr();
            n = this.L.getColumnDimension();
            double[][] LVs = Matlab.full(this.L).getData();
            double[] LRow_i = (double[])null;
            y = new double[n];
            v = 0.0D;

            for(UVs = 0; UVs < n; ++UVs) {
                v = d[UVs];
                LRow_i = LVs[UVs];

                for(int URow_i = 0; URow_i < UVs; ++URow_i) {
                    v -= LRow_i[URow_i] * y[URow_i];
                }

                y[UVs] = v;
            }

            double[][] var21 = Matlab.full(this.U).getData();
            double[] var22 = (double[])null;
            x = new double[n];
            v = 0.0D;

            for(i = n - 1; i > -1; --i) {
                var22 = var21[i];
                v = y[i];

                for(ir = n - 1; ir > i; --ir) {
                    v -= var22[ir] * x[ir];
                }

                x[i] = v;
            }

            res = new DenseVector(x);
        } else if(this.L instanceof SparseMatrix) {
            d = Matlab.full(this.P.operate(b)).getPr();
            n = this.L.getColumnDimension();
            Vector[] var19 = Matlab.sparseMatrix2SparseRowVectors(this.L);
            Vector var20 = null;
            y = new double[n];
            v = 0.0D;

            for(UVs = 0; UVs < n; ++UVs) {
                v = d[UVs];
                var20 = var19[UVs];
                int[] var24 = ((SparseVector)var20).getIr();
                x = ((SparseVector)var20).getPr();
                i = ((SparseVector)var20).getNNZ();
                boolean var26 = true;

                for(int pr = 0; pr < i; ++pr) {
                    ir = var24[pr];
                    if(ir >= UVs) {
                        break;
                    }

                    v -= x[pr] * y[ir];
                }

                y[UVs] = v;
            }

            Vector[] var23 = Matlab.sparseMatrix2SparseRowVectors(this.U);
            Vector var25 = null;
            x = new double[n];
            v = 0.0D;

            for(i = n - 1; i > -1; --i) {
                var25 = var23[i];
                v = y[i];
                int[] var27 = ((SparseVector)var25).getIr();
                double[] var28 = ((SparseVector)var25).getPr();
                int nnz = ((SparseVector)var25).getNNZ();
                boolean idx = true;
                int k = nnz - 1;

                while(true) {
                    int var29 = var27[k];
                    if(var29 <= i) {
                        x[i] = v / var25.get(i);
                        break;
                    }

                    v -= var28[k] * x[var29];
                    --k;
                }
            }

            res = new DenseVector(x);
        }

        return res;
    }

    public Matrix solve(double[][] B) {
        return this.solve((Matrix)(new DenseMatrix(B)));
    }

    public Matrix solve(Matrix B) {
        DenseMatrix res = null;
        double[][] D;
        double[] DRow_i;
        int n;
        double[][] Y;
        double[] YRow_i;
        double v;
        int UVs;
        double[] XRow_i;
        int i;
        int ir;
        int pr;
        double[][] var29;
        if(this.L instanceof DenseMatrix) {
            D = Matlab.full(this.P.mtimes(B)).getData();
            DRow_i = (double[])null;
            n = this.L.getColumnDimension();
            double[][] LVs = Matlab.full(this.L).getData();
            double[] LRow_i = (double[])null;
            Y = ArrayOperator.allocate2DArray(n, B.getColumnDimension(), 0.0D);
            YRow_i = (double[])null;
            v = 0.0D;

            for(UVs = 0; UVs < n; ++UVs) {
                LRow_i = LVs[UVs];
                DRow_i = D[UVs];
                YRow_i = Y[UVs];

                for(int URow_i = 0; URow_i < B.getColumnDimension(); ++URow_i) {
                    v = DRow_i[URow_i];

                    for(int X = 0; X < UVs; ++X) {
                        v -= LRow_i[X] * Y[X][URow_i];
                    }

                    YRow_i[URow_i] = v;
                }
            }

            double[][] var25 = Matlab.full(this.U).getData();
            double[] var26 = (double[])null;
            var29 = ArrayOperator.allocate2DArray(n, B.getColumnDimension(), 0.0D);
            XRow_i = (double[])null;

            for(i = n - 1; i > -1; --i) {
                var26 = var25[i];
                YRow_i = Y[i];
                XRow_i = var29[i];

                for(ir = 0; ir < B.getColumnDimension(); ++ir) {
                    v = YRow_i[ir];

                    for(pr = n - 1; pr > i; --pr) {
                        v -= var26[pr] * var29[pr][ir];
                    }

                    XRow_i[ir] = v / var26[i];
                }
            }

            res = new DenseMatrix(var29);
        } else if(this.L instanceof SparseMatrix) {
            D = Matlab.full(this.P.mtimes(B)).getData();
            DRow_i = (double[])null;
            n = this.L.getColumnDimension();
            Vector[] var23 = Matlab.sparseMatrix2SparseRowVectors(this.L);
            Vector var24 = null;
            Y = ArrayOperator.allocate2DArray(n, B.getColumnDimension(), 0.0D);
            YRow_i = (double[])null;
            v = 0.0D;

            for(UVs = 0; UVs < n; ++UVs) {
                var24 = var23[UVs];
                int[] var28 = ((SparseVector)var24).getIr();
                double[] var31 = ((SparseVector)var24).getPr();
                int var32 = ((SparseVector)var24).getNNZ();
                boolean var33 = true;
                DRow_i = D[UVs];
                YRow_i = Y[UVs];

                for(ir = 0; ir < B.getColumnDimension(); ++ir) {
                    v = DRow_i[ir];

                    for(pr = 0; pr < var32; ++pr) {
                        i = var28[pr];
                        if(i >= UVs) {
                            break;
                        }

                        v -= var31[pr] * Y[i][ir];
                    }

                    YRow_i[ir] = v;
                }
            }

            Vector[] var27 = Matlab.sparseMatrix2SparseRowVectors(this.U);
            Vector var30 = null;
            var29 = ArrayOperator.allocate2DArray(n, B.getColumnDimension(), 0.0D);
            XRow_i = (double[])null;

            for(i = n - 1; i > -1; --i) {
                var30 = var27[i];
                int[] var34 = ((SparseVector)var30).getIr();
                double[] var35 = ((SparseVector)var30).getPr();
                int nnz = ((SparseVector)var30).getNNZ();
                boolean idx = true;
                YRow_i = Y[i];
                XRow_i = var29[i];

                for(int j = 0; j < B.getColumnDimension(); ++j) {
                    v = YRow_i[j];
                    int k = nnz - 1;

                    while(true) {
                        int var36 = var34[k];
                        if(var36 <= i) {
                            XRow_i[j] = v / var30.get(i);
                            break;
                        }

                        v -= var35[k] * var29[var36][j];
                        --k;
                    }
                }
            }

            res = new DenseMatrix(var29);
        }

        return res;
    }

    public Matrix inverse() {
        if(this.U == null) {
            return null;
        } else {
            int n = this.L.getColumnDimension();
            double[][] AInverseTransposeData = new double[n][];
            double[][] eye = new double[n][];

            int i;
            for(i = 0; i < n; ++i) {
                eye[i] = ArrayOperator.allocateVector(n, 0.0D);
                eye[i][i] = 1.0D;
            }

            for(i = 0; i < n; ++i) {
                AInverseTransposeData[i] = Matlab.full(this.solve(eye[i])).getPr();
            }

            return (new DenseMatrix(AInverseTransposeData)).transpose();
        }
    }

    public double det() {
        if(this.U == null) {
            return 0.0D;
        } else {
            double s = 1.0D;

            for(int k = 0; k < this.U.getColumnDimension(); ++k) {
                s *= this.U.getEntry(k, k);
                if(s == 0.0D) {
                    break;
                }
            }

            return this.numRowExchange % 2 == 0?s:-s;
        }
    }

    public int getLfilling(){
        return Ufilling;
    }

    public int getUfilling(){
        return Lfilling;
    }
}