package decomposition;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.ArrayOperator;
import utils.Matlab;
import utils.Pair;
import utils.Printer;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.util.TreeMap;

/**
 * Created by Admin on 08.03.2016.
 */
public class QRDecomposition {
    private Matrix Q;
    private Matrix R;
    private Matrix P;


    public Matrix getQ() {
        return this.Q;
    }

    public Matrix getR() {
        return this.R;
    }

    public Matrix getP() {
        return this.P;
    }

    public QRDecomposition(Matrix A) {
        Matrix[] QRP = this.run(A);
        this.Q = QRP[0];
        this.R = QRP[1];
        this.P = QRP[2];
    }

    public QRDecomposition(double[][] A) {
        Matrix[] QRP = this.run(new DenseMatrix(A));
        this.Q = QRP[0];
        this.R = QRP[1];
        this.P = QRP[2];
    }

    private Matrix[] run(Matrix A) {
        return decompose(A,true);
    }

    public Vector solve(double[] b) {
        return this.solve((Vector)(new DenseVector(b)));
    }

    public Vector solve(Vector b) {
        double[] d = Matlab.full(this.Q.transpose().operate(b)).getPr();
        int rank = 0;
        int m = this.R.getRowDimension();
        int n = this.R.getColumnDimension();
        int len = Math.min(m, n);

        for(int y = 0; y < len; ++y) {
            if(this.R.getEntry(y, y) == 0.0D) {
                rank = y;
                break;
            }

            ++rank;
        }

        double[] var18 = ArrayOperator.allocate1DArray(n, 0.0D);
        double v;
        int i;
        if(this.R instanceof DenseMatrix) {
            double[][] x = ((DenseMatrix)this.R).getData();
            double[] RRow_i = (double[])null;
            v = 0.0D;

            for(i = rank - 1; i > -1; --i) {
                RRow_i = x[i];
                v = d[i];

                for(int ir = n - 1; ir > i; --ir) {
                    v -= RRow_i[ir] * var18[ir];
                }

                var18[i] = v / RRow_i[i];
            }
        } else if(this.R instanceof SparseMatrix) {
            Vector[] var19 = Matlab.sparseMatrix2SparseRowVectors(this.R);
            Vector var21 = null;
            v = 0.0D;

            for(i = rank - 1; i > -1; --i) {
                var21 = var19[i];
                v = d[i];
                int[] var22 = ((SparseVector)var21).getIr();
                double[] pr = ((SparseVector)var21).getPr();
                int nnz = ((SparseVector)var21).getNNZ();
                boolean idx = true;
                int k = nnz - 1;

                while(true) {
                    int var23 = var22[k];
                    if(var23 <= i) {
                        var18[i] = v / var21.get(i);
                        break;
                    }

                    v -= pr[k] * var18[var23];
                    --k;
                }
            }
        }

        Vector var20 = this.P.operate(new DenseVector(var18));
        return var20;
    }

    public Matrix solve(double[][] B) {
        return this.solve((Matrix)(new DenseMatrix(B)));
    }

    public Matrix solve(Matrix B) {
        double[][] D = Matlab.full(this.Q.transpose().mtimes(B)).getData();
        double[] DRow_i = (double[])null;
        int rank = 0;
        int m = this.R.getRowDimension();
        int n = this.R.getColumnDimension();
        int len = Math.min(m, n);

        for(int Y = 0; Y < len; ++Y) {
            if(this.R.getEntry(Y, Y) == 0.0D) {
                rank = Y;
                break;
            }

            ++rank;
        }

        double[][] var21 = ArrayOperator.allocate2DArray(n, B.getColumnDimension(), 0.0D);
        double[] YRow_i = (double[])null;
        double v;
        int i;
        int j;
        if(this.R instanceof DenseMatrix) {
            double[][] X = ((DenseMatrix)this.R).getData();
            double[] RRow_i = (double[])null;
            v = 0.0D;

            for(i = rank - 1; i > -1; --i) {
                RRow_i = X[i];
                DRow_i = D[i];
                YRow_i = var21[i];

                for(j = 0; j < B.getColumnDimension(); ++j) {
                    v = DRow_i[j];

                    for(int ir = n - 1; ir > i; --ir) {
                        v -= RRow_i[ir] * var21[ir][j];
                    }

                    YRow_i[j] = v / RRow_i[i];
                }
            }
        } else if(this.R instanceof SparseMatrix) {
            Vector[] var22 = Matlab.sparseMatrix2SparseRowVectors(this.R);
            Vector var24 = null;
            v = 0.0D;

            for(i = rank - 1; i > -1; --i) {
                var24 = var22[i];
                DRow_i = D[i];
                YRow_i = var21[i];

                for(j = 0; j < B.getColumnDimension(); ++j) {
                    v = DRow_i[j];
                    int[] var25 = ((SparseVector)var24).getIr();
                    double[] pr = ((SparseVector)var24).getPr();
                    int nnz = ((SparseVector)var24).getNNZ();
                    boolean idx = true;
                    int k = nnz - 1;

                    while(true) {
                        int var26 = var25[k];
                        if(var26 <= i) {
                            YRow_i[j] = v / var24.get(i);
                            break;
                        }

                        v -= pr[k] * var21[var26][j];
                        --k;
                    }
                }
            }
        }

        Matrix var23 = this.P.mtimes(new DenseMatrix(var21));
        return var23;
    }

    public static Matrix[] decompose(Matrix A, boolean pivoting) {
        A = A.copy();
        int m = A.getRowDimension();
        int n = A.getColumnDimension();
        Matrix[] QRP = new Matrix[3];
        double[] d = ArrayOperator.allocateVector(n, 0.0D);
        Vector[] PVs = Matlab.sparseMatrix2SparseColumnVectors(new SparseMatrix(n, n));

        for(int AVs = 0; AVs < n; ++AVs) {
            PVs[AVs].set(AVs, 1.0D);
        }

        double[] c;
        int j;
        int i;
        double maxVal;
        int s;
        int idx;
        int r;
        double var37;
        Vector var38;
        double var39;
        if(A instanceof DenseMatrix) {
            double[][] var34 = ((DenseMatrix)A).getData();
            c = ArrayOperator.allocateVector(n, 0.0D);

            for(j = 0; j < n && j < m; ++j) {
                for(i = j; i < n; ++i) {
                    maxVal = 0.0D;

                    for(s = j; s < m; ++s) {
                        maxVal += Math.pow(var34[s][i], 2.0D);
                    }

                    c[i] = maxVal;
                }

                i = j;
                maxVal = c[j];

                for(s = j + 1; s < n; ++s) {
                    if(maxVal < c[s]) {
                        i = s;
                        maxVal = c[s];
                    }
                }

                if(maxVal == 0.0D) {
                    System.out.println("Rank(A) < n.");
                    QRP[0] = computeQ(A);
                    QRP[1] = computeR(A, d);
                    QRP[2] = Matlab.sparseRowVectors2SparseMatrix(PVs);
                    return QRP;
                }

                if(i != j && pivoting) {
                    var37 = 0.0D;

                    for(int A_j = 0; A_j < m; ++A_j) {
                        var37 = var34[A_j][i];
                        var34[A_j][i] = var34[A_j][j];
                        var34[A_j][j] = var37;
                    }

                    var37 = c[i];
                    c[i] = c[j];
                    c[j] = var37;
                    var38 = PVs[i];
                    PVs[i] = PVs[j];
                    PVs[j] = var38;
                }

                var37 = Math.sqrt(c[j]);
                d[j] = var34[j][j] > 0.0D?-var37:var37;
                var39 = Math.sqrt(var37 * (var37 + Math.abs(var34[j][j])));
                var34[j][j] -= d[j];

                for(idx = j; idx < m; ++idx) {
                    var34[idx][j] /= var39;
                }

                for(idx = j + 1; idx < n; ++idx) {
                    var37 = 0.0D;

                    for(r = j; r < m; ++r) {
                        var37 += var34[r][j] * var34[r][idx];
                    }

                    for(r = j; r < m; ++r) {
                        var34[r][idx] -= var37 * var34[r][j];
                    }
                }
            }
        } else if(A instanceof SparseMatrix) {
            Vector[] var35 = Matlab.sparseMatrix2SparseColumnVectors(A);
            c = ArrayOperator.allocateVector(n, 0.0D);

            for(j = 0; j < n && j < m; ++j) {
                for(i = j; i < n; ++i) {
                    SparseVector var36 = (SparseVector)var35[i];
                    double[] pr = var36.getPr();
                    int[] var40 = var36.getIr();
                    int nnz = var36.getNNZ();
                    var39 = 0.0D;
                    boolean var42 = true;

                    for(r = 0; r < nnz; ++r) {
                        idx = var40[r];
                        if(idx >= j) {
                            var39 += Math.pow(pr[r], 2.0D);
                        }
                    }

                    c[i] = var39;
                }

                i = j;
                maxVal = c[j];

                for(s = j + 1; s < n; ++s) {
                    if(maxVal < c[s]) {
                        i = s;
                        maxVal = c[s];
                    }
                }

                if(maxVal == 0.0D) {
                    System.out.println("Rank(A) < n.");
                    QRP[0] = computeQ(A);
                    QRP[1] = computeR(A, d);
                    QRP[2] = Matlab.sparseRowVectors2SparseMatrix(PVs);
                    return QRP;
                }

                if(i != j && pivoting) {
                    var37 = 0.0D;
                    var38 = null;
                    var38 = var35[i];
                    var35[i] = var35[j];
                    var35[j] = var38;
                    var37 = c[i];
                    c[i] = c[j];
                    c[j] = var37;
                    var38 = PVs[i];
                    PVs[i] = PVs[j];
                    PVs[j] = var38;
                }

                var37 = Math.sqrt(c[j]);
                SparseVector var41 = (SparseVector)var35[j];
                double Ajj = var41.get(j);
                d[j] = Ajj > 0.0D?-var37:var37;
                double var43 = Math.sqrt(var37 * (var37 + Math.abs(Ajj)));
                var41.set(j, Ajj - d[j]);
                int[] ir = var41.getIr();
                double[] pr1 = var41.getPr();
                int nnz1 = var41.getNNZ();
                boolean idx1 = false;

                int k;
                for(k = 0; k < nnz1; ++k) {
                    int var44 = ir[k];
                    if(var44 >= j) {
                        pr1[k] /= var43;
                    }
                }

                for(k = j + 1; k < n; ++k) {
                    SparseVector A_k = (SparseVector)var35[k];
                    var37 = 0.0D;
                    int[] ir2 = A_k.getIr();
                    double[] pr2 = A_k.getPr();
                    int nnz2 = A_k.getNNZ();
                    int t;
                    if(nnz1 != 0 && nnz2 != 0) {
                        t = 0;
                        int k2 = 0;
                        boolean r1 = false;
                        boolean r2 = false;
                        double v = 0.0D;

                        while(t < nnz1 && k2 < nnz2) {
                            int var45 = ir[t];
                            int var46 = ir2[k2];
                            if(var45 < var46) {
                                ++t;
                            } else if(var45 == var46) {
                                v = pr1[t] * pr2[k2];
                                ++t;
                                ++k2;
                                if(var45 >= j) {
                                    var37 += v;
                                }
                            } else {
                                ++k2;
                            }
                        }
                    }

                    for(t = j; t < m; ++t) {
                        A_k.set(t, A_k.get(t) - var37 * var41.get(t));
                    }
                }
            }

            A = Matlab.sparseColumnVectors2SparseMatrix(var35);
        }

        QRP[0] = computeQ(A);
        QRP[1] = computeR(A, d);
        QRP[2] = Matlab.sparseColumnVectors2SparseMatrix(PVs);
        return QRP;
    }

    private static Matrix computeQ(Matrix A) {
        int m = A.getRowDimension();
        int n = A.getColumnDimension();
        double[][] Q = (new DenseMatrix(m, m, 0.0D)).getData();
        double s = 0.0D;
        double[] y = (double[])null;

        for(int i = 0; i < m; ++i) {
            y = Q[i];
            y[i] = 1.0D;

            for(int j = 0; j < n; ++j) {
                s = 0.0D;

                int k;
                for(k = j; k < m; ++k) {
                    s += A.getEntry(k, j) * y[k];
                }

                for(k = j; k < m; ++k) {
                    y[k] -= A.getEntry(k, j) * s;
                }
            }
        }

        return new DenseMatrix(Q);
    }

    private static Matrix computeR(Matrix A, double[] d) {
        int m = A.getRowDimension();
        int n = A.getColumnDimension();
        Object R = null;
        int i;
        if(A instanceof DenseMatrix) {
            double[][] map = ((DenseMatrix)A).getData();

            for(i = 0; i < m; ++i) {
                double[] j = map[i];
                if(i < n) {
                    j[i] = d[i];
                }

                int len = Math.min(i, n);

                for(int j1 = 0; j1 < len; ++j1) {
                    j[j1] = 0.0D;
                }
            }

            R = A;
        } else if(A instanceof SparseMatrix) {
            TreeMap var10 = new TreeMap();

            for(i = 0; i < m; ++i) {
                if(i < n) {
                    var10.put(Pair.of(Integer.valueOf(i), Integer.valueOf(i)), Double.valueOf(d[i]));
                }

                for(int var11 = i + 1; var11 < n; ++var11) {
                    var10.put(Pair.of(Integer.valueOf(i), Integer.valueOf(var11)), Double.valueOf(A.getEntry(i, var11)));
                }
            }

            R = SparseMatrix.createSparseMatrix(var10, m, n);
        }

        return (Matrix)R;
    }
}