package vector;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.ArrayOperator;
import utils.Pair;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by Admin on 08.03.2016.
 */
public class SparseVector implements Vector, Serializable {
    private static final long serialVersionUID = 4760099385084043335L;
    private int[] ir;
    private double[] pr;
    private int nnz;
    private int dim;

    public SparseVector(int dim) {
        this.ir = new int[0];
        this.pr = new double[0];
        this.nnz = 0;
        this.dim = dim;
    }

    public SparseVector(int[] ir, double[] pr, int nnz, int dim) {
        this.ir = ir;
        this.pr = pr;
        this.nnz = nnz;
        this.dim = dim;
    }

    public void assignSparseVector(SparseVector V) {
        this.ir = (int[])V.ir.clone();
        this.pr = (double[])V.pr.clone();
        this.nnz = V.nnz;
        this.dim = V.dim;
    }

    public int[] getIr() {
        return this.ir;
    }

    public double[] getPr() {
        return this.pr;
    }

    public int getNNZ() {
        return this.nnz;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer(100);
        sb.append('[');

        for(int k = 0; k < this.nnz; ++k) {
            sb.append(String.format("%d: %.4f", new Object[]{Integer.valueOf(this.ir[k]), Double.valueOf(this.pr[k])}));
            if(k < this.nnz - 1) {
                sb.append(", ");
            }
        }

        sb.append(']');
        return sb.toString();
    }

    public int getDim() {
        return this.dim;
    }

    public void setDim(int dim) {
        if(dim > this.dim) {
            this.dim = dim;
        }

    }

    public Vector copy() {
        return new SparseVector((int[])this.ir.clone(), (double[])this.pr.clone(), this.nnz, this.dim);
    }

    public Vector clone() {
        return this.copy();
    }

    public Vector times(Vector V) {
        if(V.getDim() != this.dim) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if(V instanceof DenseVector) {
            return V.times(this);
        } else if(!(V instanceof SparseVector)) {
            return null;
        } else {
            ArrayList list = new ArrayList();
            int[] ir = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            int nnz2 = ((SparseVector)V).getNNZ();
            int nnz;
            int dim;
            if(this.nnz != 0 && nnz2 != 0) {
                nnz = 0;
                dim = 0;
                boolean ir_res = false;
                boolean pr_res = false;
                double k = 0.0D;
                boolean i = true;

                while(nnz < this.nnz && dim < nnz2) {
                    int var13 = this.ir[nnz];
                    int var15 = ir[dim];
                    if(var13 < var15) {
                        ++nnz;
                    } else if(var13 == var15) {
                        k = this.pr[nnz] * pr[dim];
                        ++nnz;
                        ++dim;
                        if(k != 0.0D) {
                            list.add(Pair.of(Integer.valueOf(var13), Double.valueOf(k)));
                        }
                    } else {
                        ++dim;
                    }
                }
            }

            nnz = list.size();
            dim = this.getDim();
            int[] var14 = new int[nnz];
            double[] var16 = new double[nnz];
            int var17 = 0;

            for(Iterator var18 = list.iterator(); var18.hasNext(); ++var17) {
                Pair pair = (Pair)var18.next();
                var14[var17] = ((Integer)pair.first).intValue();
                var16[var17] = ((Double)pair.second).doubleValue();
            }

            return new SparseVector(var14, var16, nnz, dim);
        }
    }

    public Vector plus(Vector V) {
        if(V.getDim() != this.dim) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if(V instanceof DenseVector) {
            return V.plus(this);
        } else if(!(V instanceof SparseVector)) {
            return null;
        } else {
            ArrayList list = new ArrayList();
            int[] ir = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            int nnz2 = ((SparseVector)V).getNNZ();
            int nnz;
            int dim;
            if(this.nnz != 0 || nnz2 != 0) {
                nnz = 0;
                dim = 0;
                boolean ir_res = false;
                boolean pr_res = false;
                double k = 0.0D;
                boolean i = true;

                while(nnz < this.nnz || dim < nnz2) {
                    int var18;
                    if(dim == nnz2) {
                        var18 = this.ir[nnz];
                        k = this.pr[nnz];
                        ++nnz;
                    } else if(nnz == this.nnz) {
                        var18 = ir[dim];
                        k = pr[dim];
                        ++dim;
                    } else {
                        int var13 = this.ir[nnz];
                        int var15 = ir[dim];
                        if(var13 < var15) {
                            var18 = var13;
                            k = this.pr[nnz];
                            ++nnz;
                        } else if(var13 == var15) {
                            var18 = var13;
                            k = this.pr[nnz] + pr[dim];
                            ++nnz;
                            ++dim;
                        } else {
                            var18 = var15;
                            k = pr[dim];
                            ++dim;
                        }
                    }

                    if(k != 0.0D) {
                        list.add(Pair.of(Integer.valueOf(var18), Double.valueOf(k)));
                    }
                }
            }

            nnz = list.size();
            dim = this.getDim();
            int[] var14 = new int[nnz];
            double[] var16 = new double[nnz];
            int var17 = 0;

            for(Iterator var19 = list.iterator(); var19.hasNext(); ++var17) {
                Pair pair = (Pair)var19.next();
                var14[var17] = ((Integer)pair.first).intValue();
                var16[var17] = ((Double)pair.second).doubleValue();
            }

            return new SparseVector(var14, var16, nnz, dim);
        }
    }

    public Vector minus(Vector V) {
        if(V.getDim() != this.dim) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        if(V instanceof DenseVector) {
            double[] var13 = ((DenseVector)V).getPr();
            double[] var14 = new double[this.dim];

            int var15;
            for(var15 = 0; var15 < this.dim; ++var15) {
                var14[var15] = -var13[var15];
            }

            for(var15 = 0; var15 < this.nnz; ++var15) {
                var14[this.ir[var15]] += this.pr[var15];
            }

            return new DenseVector(var14);
        } else if(!(V instanceof SparseVector)) {
            return null;
        } else {
            ArrayList list = new ArrayList();
            int[] ir = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            int nnz2 = ((SparseVector)V).getNNZ();
            int nnz;
            int dim;
            if(this.nnz != 0 || nnz2 != 0) {
                nnz = 0;
                dim = 0;
                boolean ir_res = false;
                boolean pr_res = false;
                double k = 0.0D;
                boolean i = true;

                while(nnz < this.nnz || dim < nnz2) {
                    int var21;
                    if(dim == nnz2) {
                        var21 = this.ir[nnz];
                        k = this.pr[nnz];
                        ++nnz;
                    } else if(nnz == this.nnz) {
                        var21 = ir[dim];
                        k = -pr[dim];
                        ++dim;
                    } else {
                        int var16 = this.ir[nnz];
                        int var18 = ir[dim];
                        if(var16 < var18) {
                            var21 = var16;
                            k = this.pr[nnz];
                            ++nnz;
                        } else if(var16 == var18) {
                            var21 = var16;
                            k = this.pr[nnz] - pr[dim];
                            ++nnz;
                            ++dim;
                        } else {
                            var21 = var18;
                            k = -pr[dim];
                            ++dim;
                        }
                    }

                    if(k != 0.0D) {
                        list.add(Pair.of(Integer.valueOf(var21), Double.valueOf(k)));
                    }
                }
            }

            nnz = list.size();
            dim = this.getDim();
            int[] var17 = new int[nnz];
            double[] var19 = new double[nnz];
            int var20 = 0;

            for(Iterator var22 = list.iterator(); var22.hasNext(); ++var20) {
                Pair pair = (Pair)var22.next();
                var17[var20] = ((Integer)pair.first).intValue();
                var19[var20] = ((Double)pair.second).doubleValue();
            }

            return new SparseVector(var17, var19, nnz, dim);
        }
    }

    public double get(int i) {
        if(i < 0 || i >= this.dim) {
            System.err.println("Wrong index.");
            System.exit(1);
        }

        if(this.nnz == 0) {
            return 0.0D;
        } else {
            int u = this.nnz - 1;
            int l = 0;
            boolean idx = true;
            boolean k = false;

            while(l <= u) {
                int k1 = (u + l) / 2;
                int idx1 = this.ir[k1];
                if(idx1 == i) {
                    return this.pr[k1];
                }

                if(idx1 < i) {
                    l = k1 + 1;
                } else {
                    u = k1 - 1;
                }
            }

            return 0.0D;
        }
    }

    public void set(int i, double v) {
        if(i < 0 || i >= this.dim) {
            System.err.println("Wrong index.");
            System.exit(1);
        }

        if(this.nnz == 0) {
            this.insertEntry(i, v, 0);
        } else {
            int u = this.nnz - 1;
            int l = 0;
            boolean idx = true;
            int k = 0;
            byte flag = 0;

            while(l <= u) {
                k = (u + l) / 2;
                int var9 = this.ir[k];
                if(var9 == i) {
                    if(v == 0.0D) {
                        this.deleteEntry(k);
                    } else {
                        this.pr[k] = v;
                    }

                    return;
                }

                if(var9 < i) {
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

            this.insertEntry(i, v, k);
        }
    }

    private void insertEntry(int r, double v, int pos) {
        if(v != 0.0D) {
            int len_old = this.pr.length;
            int new_space = len_old < this.dim - 10?10:this.dim - len_old;
            int var8;
            if(this.nnz + 1 > len_old) {
                double[] i = new double[len_old + new_space];
                System.arraycopy(this.pr, 0, i, 0, pos);
                i[pos] = v;
                if(pos < len_old) {
                    System.arraycopy(this.pr, pos, i, pos + 1, len_old - pos);
                }

                this.pr = i;
            } else {
                for(var8 = this.nnz - 1; var8 >= pos; --var8) {
                    this.pr[var8 + 1] = this.pr[var8];
                }

                this.pr[pos] = v;
            }

            if(this.nnz + 1 > len_old) {
                int[] var9 = new int[len_old + new_space];
                System.arraycopy(this.ir, 0, var9, 0, pos);
                var9[pos] = r;
                if(pos < len_old) {
                    System.arraycopy(this.ir, pos, var9, pos + 1, len_old - pos);
                }

                this.ir = var9;
            } else {
                for(var8 = this.nnz - 1; var8 >= pos; --var8) {
                    this.ir[var8 + 1] = this.ir[var8];
                }

                this.ir[pos] = r;
            }

            ++this.nnz;
        }
    }

    public void clean() {
        for(int k = this.nnz - 1; k >= 0; --k) {
            if(this.pr[k] == 0.0D) {
                this.deleteEntry(k);
            }
        }

    }

    private void deleteEntry(int pos) {
        for(int i = pos; i < this.nnz - 1; ++i) {
            this.pr[i] = this.pr[i + 1];
            this.ir[i] = this.ir[i + 1];
        }

        --this.nnz;
    }

    public Vector operate(Matrix A) {
        int M = A.getRowDimension();
        int N = A.getColumnDimension();
        if(M != this.dim) {
            System.err.println("Dimension doesn\'t match.");
            System.exit(1);
        }

        double[] pr;
        double s;
        int var20;
        int var21;
        int var22;
        if(A instanceof DenseMatrix) {
            double[] var18 = ArrayOperator.allocate1DArray(N, 0.0D);
            double[][] var19 = ((DenseMatrix)A).getData();
            pr = (double[])null;
            s = 0.0D;

            for(var20 = 0; var20 < this.nnz; ++var20) {
                var21 = this.ir[var20];
                pr = var19[var21];
                s = this.pr[var20];

                for(var22 = 0; var22 < N; ++var22) {
                    var18[var22] += s * pr[var22];
                }
            }

            return new DenseVector(var18);
        } else if(!(A instanceof SparseMatrix)) {
            return null;
        } else {
            int[] ir = ((SparseMatrix)A).getIr();
            int[] jc = ((SparseMatrix)A).getJc();
            pr = ((SparseMatrix)A).getPr();
            s = 0.0D;
            boolean k1 = false;
            boolean k2 = false;
            boolean c = false;
            boolean r = false;
            TreeMap map = new TreeMap();

            int nnz;
            for(nnz = 0; nnz < N; ++nnz) {
                var20 = 0;
                var21 = jc[nnz];
                s = 0.0D;

                while(var21 < jc[nnz + 1] && var20 < this.nnz) {
                    var22 = this.ir[var20];
                    int var23 = ir[var21];
                    if(var23 < var22) {
                        ++var21;
                    } else if(var23 > var22) {
                        ++var20;
                    } else {
                        s += this.pr[var20] * pr[var21];
                        ++var20;
                        ++var21;
                    }
                }

                if(s != 0.0D) {
                    map.put(Integer.valueOf(nnz), Double.valueOf(s));
                }
            }

            nnz = map.size();
            ir = new int[nnz];
            pr = new double[nnz];
            int ind = 0;

            for(Iterator var17 = map.entrySet().iterator(); var17.hasNext(); ++ind) {
                Map.Entry entry = (Map.Entry)var17.next();
                ir[ind] = ((Integer)entry.getKey()).intValue();
                pr[ind] = ((Double)entry.getValue()).doubleValue();
            }

            return new SparseVector(ir, pr, nnz, N);
        }
    }

    public void clear() {
        this.ir = new int[0];
        this.pr = new double[0];
        this.nnz = 0;
    }

    public Vector times(double v) {
        if(v == 0.0D) {
            return new SparseVector(this.dim);
        } else {
            SparseVector res = (SparseVector)this.copy();
            double[] pr = res.pr;

            for(int k = 0; k < this.nnz; ++k) {
                pr[k] *= v;
            }

            return res;
        }
    }

    public double norm1(){
        double sum=0;
        for (int i=0;i<this.dim;i++){
            sum+=Math.abs(pr[i]);
        }
        return sum;
    }

    public double norm2(){
        double sum=0;
        for (int i=0;i<this.dim;i++){
            sum+=Math.abs(pr[i]*pr[i]);
        }
        return Math.sqrt(sum);
    }

    public double normInf(){
        double temp=0;
        for (int i = 0; i < this.getDim(); i++) {
            temp = Math.abs(pr[i])>temp ? Math.abs(pr[i]) : temp;
        }
        return temp;
    }
}
