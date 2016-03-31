package decomposition;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.ArrayOperator;
import utils.Matlab;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

/**
 * Created by Admin on 13.03.2016.
 */
public class CholesskyDecomposition {
    private Matrix L;

    public Matrix getL() {
        return this.L;
    }

    public Matrix getLt() {
        return this.L.transpose();
    }

    public CholesskyDecomposition(Matrix A){
        this.L = run(A);
    }

    private Matrix run(Matrix A){
        int n = A.getRowDimension();
        if(n != A.getColumnDimension()) {
            System.err.println("A should be a square matrix.");
            System.exit(1);
        }
        Matrix L;

        if(A instanceof DenseMatrix) {
            L = new DenseMatrix(n,n,0.0D);

            double temp;
            double[][] AVs = ((DenseMatrix)A.copy()).getData();
            for(int j=0;j<n;j++){
                temp=0;
                for (int k = 0; k <= j-1; k++) {
                    temp+=L.getEntry(j,k)*L.getEntry(j,k);
                }
                L.setEntry(j,j,Math.sqrt(AVs[j][j]-temp));

                for (int i = j+1; i < n; i++) {
                    temp=0;
                    for (int k = 0; k <= j - 1; k++) {
                        temp+=L.getEntry(i,k)*L.getEntry(j,k);
                    }
                    L.setEntry(i,j,(AVs[i][j]-temp)/L.getEntry(j,j));
                }
            }
        } else /*if(A instanceof SparseMatrix)*/ {
            L = new SparseMatrix(n,n);

            double temp;
            for(int j=0;j<n;j++){
                temp=0;
                for (int k = 0; k <= j-1; k++) {
                    temp+=L.getEntry(j,k)*L.getEntry(j,k);
                }
                L.setEntry(j,j,Math.sqrt(A.getEntry(j,j)-temp));

                for (int i = j+1; i < n; i++) {
                    temp=0;
                    for (int k = 0; k <= j - 1; k++) {
                        temp+=L.getEntry(i,k)*L.getEntry(j,k);
                    }
                    L.setEntry(i,j,(A.getEntry(i,j)-temp)/L.getEntry(j,j));
                }
            }
        }
        return L;
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
            d = Matlab.full(b).getPr();
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
//TODO
                y[UVs] = v/LRow_i[UVs];
            }

            double[][] var21 = Matlab.full(this.L.transpose()).getData();
            double[] var22 = (double[])null;
            x = new double[n];
            v = 0.0D;

            for(i = n - 1; i > -1; --i) {
                var22 = var21[i];
                v = y[i];

                for(ir = n - 1; ir > i; --ir) {
                    v -= var22[ir] * x[ir];
                }

                x[i] = v / var22[i];
            }

            res = new DenseVector(x);
        } else if(this.L instanceof SparseMatrix) {
            d = Matlab.full(b).getPr();
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

                y[UVs] = v/var20.get(UVs);
            }

            Vector[] var23 = Matlab.sparseMatrix2SparseRowVectors(this.L.transpose());
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
            D = Matlab.full(B).getData();
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

            double[][] var25 = Matlab.full(this.L.transpose()).getData();
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
            D = Matlab.full(B).getData();
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

            Vector[] var27 = Matlab.sparseMatrix2SparseRowVectors(this.L.transpose());
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

    public double det() {
        if (this.L == null) {
            return 0.0D;
        } else {
            double s = 1.0D;

            for (int k = 0; k < this.L.getColumnDimension(); ++k) {
                s *= this.L.getEntry(k, k);
                if (s == 0.0D) {
                    break;
                }
            }

            return s * s;
        }
    }
}
