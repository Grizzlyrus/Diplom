package vector;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.ArrayOperator;
import utils.Pair;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * Created by Admin on 08.03.2016.
 */
public class DenseVector implements Vector, Serializable {
    private static final long serialVersionUID = -6411390717530519480L;
    private double[] pr;

    public DenseVector() {
    }

    public DenseVector(int dim, double v) {
        this.pr = ArrayOperator.allocateVector(dim, v);
    }

    public DenseVector(int dim) {
        this.pr = ArrayOperator.allocateVector(dim, 0.0D);
    }

    public DenseVector(double[] pr) {
        this.pr = pr;
    }

    public static DenseVector buildDenseVector(double[] pr) {
        DenseVector res = new DenseVector();
        res.pr = pr;
        return res;
    }

    public double[] getPr() {
        return this.pr;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer(100);
        sb.append('[');

        for (int k = 0; k < this.pr.length; ++k) {
            sb.append(String.format("%.4f", new Object[]{Double.valueOf(this.pr[k])}));
            if (k < this.pr.length - 1) {
                sb.append(", ");
            }
        }

        sb.append(']');
        return sb.toString();
    }


    public int getDim() {
        return this.pr.length;
    }

    public Vector copy() {
        return new DenseVector((double[]) this.pr.clone());
    }

    public Vector clone() {
        return this.copy();
    }

    public Vector times(Vector V) {
        if (V.getDim() != this.pr.length) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if (V instanceof DenseVector) {
            return new DenseVector(ArrayOperator.times(this.pr, ((DenseVector) V).getPr()));
        } else if (!(V instanceof SparseVector)) {
            return null;
        } else {
            ArrayList list = new ArrayList();
            int[] ir = ((SparseVector) V).getIr();
            double[] pr = ((SparseVector) V).getPr();
            boolean idx = true;
            double v = 0.0D;

            int nnz;
            for (nnz = 0; nnz < ((SparseVector) V).getNNZ(); ++nnz) {
                int var15 = ir[nnz];
                v = this.pr[var15] * pr[nnz];
                if (v != 0.0D) {
                    list.add(Pair.of(Integer.valueOf(var15), Double.valueOf(v)));
                }
            }

            nnz = list.size();
            int dim = this.getDim();
            int[] ir_res = new int[nnz];
            double[] pr_res = new double[nnz];
            int k = 0;

            for (Iterator var14 = list.iterator(); var14.hasNext(); ++k) {
                Pair pair = (Pair) var14.next();
                ir_res[k] = ((Integer) pair.first).intValue();
                pr_res[k] = ((Double) pair.second).doubleValue();
            }

            return new SparseVector(ir_res, pr_res, nnz, dim);
        }
    }

    public Vector plus(Vector V) {
        if (V.getDim() != this.pr.length) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if (V instanceof DenseVector) {
            return new DenseVector(ArrayOperator.plus(this.pr, ((DenseVector) V).getPr()));
        } else if (!(V instanceof SparseVector)) {
            return null;
        } else {
            double[] res = (double[]) this.pr.clone();
            int[] ir = ((SparseVector) V).getIr();
            double[] pr = ((SparseVector) V).getPr();
            boolean idx = true;

            for (int k = 0; k < ((SparseVector) V).getNNZ(); ++k) {
                int var7 = ir[k];
                res[ir[k]] = this.pr[var7] + pr[k];
            }

            return new DenseVector(res);
        }
    }

    public Vector minus(Vector V) {
        if (V.getDim() != this.pr.length) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if (V instanceof DenseVector) {
            return new DenseVector(ArrayOperator.minus(this.pr, ((DenseVector) V).getPr()));
        } else if (!(V instanceof SparseVector)) {
            return null;
        } else {
            double[] res = (double[]) this.pr.clone();
            int[] ir = ((SparseVector) V).getIr();
            double[] pr = ((SparseVector) V).getPr();
            boolean idx = true;

            for (int k = 0; k < ((SparseVector) V).getNNZ(); ++k) {
                int var7 = ir[k];
                res[ir[k]] = this.pr[var7] - pr[k];
            }

            return new DenseVector(res);
        }
    }

    public double get(int i) {
        return this.pr[i];
    }

    public void set(int i, double v) {
        this.pr[i] = v;
    }

    public Vector operate(Matrix A) {
        int dim = this.getDim();
        int M = A.getRowDimension();
        int N = A.getColumnDimension();
        if (M != dim) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        double[] res = ArrayOperator.allocate1DArray(N, 0.0D);
        int k;
        if (A instanceof DenseMatrix) {
            double[][] ir = ((DenseMatrix) A).getData();
            double[] jc = (double[]) null;
            double pr = 0.0D;

            for (k = 0; k < M; ++k) {
                jc = ir[k];
                pr = this.pr[k];

                for (int j1 = 0; j1 < N; ++j1) {
                    res[j1] += pr * jc[j1];
                }
            }
        } else if (A instanceof SparseMatrix) {
            int[] var12 = ((SparseMatrix) A).getIr();
            int[] var13 = ((SparseMatrix) A).getJc();
            double[] var14 = ((SparseMatrix) A).getPr();

            for (int j = 0; j < N; ++j) {
                for (k = var13[j]; k < var13[j + 1]; ++k) {
                    res[j] += this.pr[var12[k]] * var14[k];
                }
            }
        }

        return new DenseVector(res);
    }

    public void clear() {
        ArrayOperator.clearVector(this.pr);
    }

    public Vector times(double v) {
        if (v == 0.0D) {
            return new DenseVector(this.getDim(), 0.0D);
        } else {
            double[] resData = (double[]) this.pr.clone();

            for (int i = 0; i < this.pr.length; ++i) {
                resData[i] *= v;
            }

            return new DenseVector(resData);
        }
    }

    public double norm1(){
        double sum=0;
        for (int i = 0; i < this.getDim(); i++) {
            sum+=Math.abs(this.get(i));
        }
        return sum;
    }

    public double norm2(){
        double sum=0;
        for (int i = 0; i < this.getDim(); i++) {
            sum+=Math.abs(this.get(i)*this.get(i));
        }
        return Math.sqrt(sum);
    }

    public double normInf(){
        double temp = 0;
        for (int i = 0; i < this.getDim(); i++) {
            temp = Math.abs(this.get(i))>temp ? Math.abs(this.get(i)) : temp;
        }
        return temp;
    }

}