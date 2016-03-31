package utils;

import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

/**
 * Created by Admin on 08.03.2016.
 */
public class Printer {
    public Printer() {
    }

    public static void printSparseMatrix(Matrix A, int p) {
        if(!(A instanceof SparseMatrix)) {
            System.err.println("SparseMatrix input is expected.");
        } else if(((SparseMatrix)A).getNNZ() == 0) {
            System.out.println("Empty sparse matrix.");
            System.out.println();
        } else {
            int nRow = A.getRowDimension();
            int nCol = A.getColumnDimension();
            String leftFormat = String.format("  %%%ds, ", new Object[]{Integer.valueOf(String.valueOf(nRow).length() + 1)});
            String rightFormat = String.format("%%-%ds", new Object[]{Integer.valueOf(String.valueOf(nCol).length() + 2)});
            String format = leftFormat + rightFormat + sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)});
            SparseMatrix S = (SparseMatrix)A;
            int[] ir = S.getIr();
            int[] jc = S.getJc();
            double[] pr = S.getPr();
            int N = S.getColumnDimension();
            String valueString = "";
            boolean i = true;

            for(int j = 0; j < N; ++j) {
                for(int k = jc[j]; k < jc[j + 1]; ++k) {
                    System.out.print("  ");
                    int var21 = ir[k];
                    double v = pr[k];
                    int rv = (int)Math.round(v);
                    if(v != (double)rv) {
                        valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                    } else {
                        valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                    }

                    String leftString = String.format("(%d", new Object[]{Integer.valueOf(var21 + 1)});
                    String rightString = String.format("%d)", new Object[]{Integer.valueOf(j + 1)});
                    System.out.println(String.format(format, new Object[]{leftString, rightString, valueString}));
                }
            }

            System.out.println();
        }
    }

    public static void printSparseMatrix(Matrix A) {
        printSparseMatrix(A, 4);
    }

    public static void printDenseMatrix(Matrix A, int p) {
        if(!(A instanceof DenseMatrix)) {
            System.err.println("DenseMatrix input is expected.");
        } else if(((DenseMatrix)A).getData() == null) {
            System.out.println("Empty matrix.");
        } else {
            for(int i = 0; i < A.getRowDimension(); ++i) {
                System.out.print("  ");

                for(int j = 0; j < A.getColumnDimension(); ++j) {
                    String valueString = "";
                    double v = A.getEntry(i, j);
                    int rv = (int)Math.round(v);
                    if(v != (double)rv) {
                        valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                    } else {
                        valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                    }

                    System.out.print(sprintf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
                    System.out.print("  ");
                }

                System.out.println();
            }

            System.out.println();
        }
    }

    public static void printDenseMatrix(Matrix A) {
        printDenseMatrix(A, 4);
    }

    public static void printMatrix(Matrix A, int p) {
        if(A == null) {
            System.out.println("Empty matrix.");
        } else {
            int rv;
            if(A instanceof SparseMatrix) {
                if(((SparseMatrix)A).getNNZ() == 0) {
                    System.out.println("Empty sparse matrix.");
                } else {
                    SparseMatrix var16 = (SparseMatrix)A;
                    int[] var17 = var16.getIc();
                    int[] var18 = var16.getJr();
                    double[] var19 = var16.getPr();
                    int[] valCSRIndices = var16.getValCSRIndices();
                    rv = var16.getRowDimension();
                    String valueString1 = "";

                    for(int r = 0; r < rv; ++r) {
                        System.out.print("  ");
                        boolean currentColumn = false;
                        int lastColumn = -1;

                        for(int k = var18[r]; k < var18[r + 1]; ++k) {
                            int var20;
                            for(var20 = var17[k]; lastColumn < var20 - 1; ++lastColumn) {
                                System.out.printf(String.format("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{" "});
                                System.out.print("  ");
                            }

                            lastColumn = var20;
                            double v1 = var19[valCSRIndices[k]];
                            int rv1 = (int)Math.round(v1);
                            if(v1 != (double)rv1) {
                                valueString1 = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v1)});
                            } else {
                                valueString1 = sprintf("%d", new Object[]{Integer.valueOf(rv1)});
                            }

                            System.out.printf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString1});
                            System.out.print("  ");
                        }

                        System.out.println();
                    }

                    System.out.println();
                }
            } else {
                if(A instanceof DenseMatrix) {
                    if(((DenseMatrix)A).getData() == null) {
                        System.out.println("Empty matrix.");
                        return;
                    }

                    for(int i = 0; i < A.getRowDimension(); ++i) {
                        System.out.print("  ");

                        for(int j = 0; j < A.getColumnDimension(); ++j) {
                            String valueString = "";
                            double v = A.getEntry(i, j);
                            rv = (int)Math.round(v);
                            if(v != (double)rv) {
                                valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                            } else {
                                valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                            }

                            System.out.print(sprintf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
                            System.out.print("  ");
                        }

                        System.out.println();
                    }

                    System.out.println();
                }

            }
        }
    }

    public static void printMatrix(Matrix A) {
        printMatrix((Matrix)A, 4);
    }

    public static void printMatrix(double[] V, int p) {
        for(int i = 0; i < V.length; ++i) {
            System.out.print("  ");
            String valueString = "";
            double v = V[i];
            int rv = (int)Math.round(v);
            if(v != (double)rv) {
                valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
            } else {
                valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
            }

            System.out.print(sprintf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
            System.out.print("  ");
            System.out.println();
        }

        System.out.println();
    }

    public static void printMatrix(double[] V) {
        printMatrix((double[])V, 4);
    }

    public static void printVector(double[] V, int p) {
        for(int i = 0; i < V.length; ++i) {
            System.out.print("  ");
            String valueString = "";
            double v = V[i];
            int rv = (int)Math.round(v);
            if(v != (double)rv) {
                valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
            } else {
                valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
            }

            System.out.print(sprintf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
        }

        System.out.println();
        System.out.println();
    }

    public static void printVector(double[] V) {
        printVector((double[])V, 4);
    }

    public static void printVector(Vector V, int p) {
        if(V instanceof DenseVector) {
            printDenseVector(V, p);
        } else {
            printSparseVector(V, p);
        }

    }

    public static void printVector(Vector V) {
        printVector((Vector)V, 15);
    }

    public static void printDenseVector(Vector V, int p) {
        if(V instanceof DenseVector) {
            int dim = V.getDim();
            double[] pr = ((DenseVector)V).getPr();

            for(int k = 0; k < dim; ++k) {
                System.out.print("  ");
                double v = pr[k];
                int rv = (int)Math.round(v);
                String valueString;
                if(v != (double)rv) {
                    valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                } else {
                    valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                }

                System.out.print(sprintf(sprintf("%%%ds", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
                System.out.println();
            }

            System.out.println();
        } else {
            System.err.println("The input vector should be a DenseVector instance");
            System.exit(1);
        }

    }

    public static void printDenseVector(Vector V) {
        printDenseVector(V, 4);
    }

    public static void printSparseVector(Vector V, int p) {
        if(V instanceof SparseVector) {
            int[] ir = ((SparseVector)V).getIr();
            double[] pr = ((SparseVector)V).getPr();
            int nnz = ((SparseVector)V).getNNZ();

            for(int k = 0; k < nnz; ++k) {
                System.out.print("  ");
                int idx = ir[k];
                double v = pr[k];
                int rv = (int)Math.round(v);
                String valueString;
                if(v != (double)rv) {
                    valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                } else {
                    valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                }

                System.out.print(sprintf(sprintf("(%d, 1)%%%ds", new Object[]{Integer.valueOf(idx + 1), Integer.valueOf(8 + p - 4)}), new Object[]{valueString}));
                System.out.println();
            }

            System.out.println();
        } else {
            System.err.println("The input vector should be a SparseVector instance");
            System.exit(1);
        }

    }

    public static void printSparseVector(Vector V) {
        printSparseVector(V, 4);
    }

    public static void display(Vector V, int p) {
        printVector(V, p);
    }

    public static void display(Vector V) {
        display((Vector)V, 4);
    }

    public static void display(double[] V, int p) {
        printVector((Vector)(new DenseVector(V)), p);
    }

    public static void display(double[] V) {
        display((double[])V, 4);
    }

    public static void display(Matrix A, int p) {
        if(A instanceof DenseMatrix) {
            printDenseMatrix(A, p);
        } else if(A instanceof SparseMatrix) {
            printSparseMatrix(A, p);
        }

    }

    public static void display(Matrix A) {
        display((Matrix)A, 4);
    }

    public static void display(double[][] A, int p) {
        printMatrix((Matrix)(new DenseMatrix(A)), p);
    }

    public static void display(double[][] A) {
        display((double[][])A, 4);
    }

    public static void disp(Vector V, int p) {
        display(V, p);
    }

    public static void disp(Vector V) {
        display((Vector)V, 4);
    }

    public static void disp(double[] V, int p) {
        display((Vector)(new DenseVector(V)), p);
    }

    public static void disp(double[] V) {
        display((Vector)(new DenseVector(V)), 4);
    }

    public static void disp(Matrix A, int p) {
        display(A, p);
    }

    public static void disp(Matrix A) {
        display((Matrix)A, 4);
    }

    public static void disp(double[][] A, int p) {
        display((Matrix)(new DenseMatrix(A)), p);
    }

    public static void disp(double[][] A) {
        display((Matrix)(new DenseMatrix(A)), 4);
    }

    public static void disp(double v) {
        System.out.print("  ");
        System.out.println(v);
    }

    public static void disp(int[] V) {
        display(V);
    }

    public static void disp(int[][] M) {
        display(M);
    }

    public static void display(int[] V) {
        if(V == null) {
            System.out.println("Empty vector!");
        } else {
            for(int i = 0; i < V.length; ++i) {
                System.out.print("  ");
                String valueString = "";
                double v = (double)V[i];
                int rv = (int)Math.round(v);
                if(v != (double)rv) {
                    valueString = String.format("%.4f", new Object[]{Double.valueOf(v)});
                } else {
                    valueString = String.format("%d", new Object[]{Integer.valueOf(rv)});
                }

                System.out.print(String.format("%7s", new Object[]{valueString}));
                System.out.print("  ");
            }

            System.out.println();
        }
    }

    public static void display(int[][] M) {
        if(M == null) {
            System.out.println("Empty matrix!");
        } else {
            for(int i = 0; i < M.length; ++i) {
                System.out.print("  ");

                for(int j = 0; j < M[0].length; ++j) {
                    String valueString = "";
                    double v = (double)M[i][j];
                    int rv = (int)Math.round(v);
                    if(v != (double)rv) {
                        valueString = String.format("%.4f", new Object[]{Double.valueOf(v)});
                    } else {
                        valueString = String.format("%d", new Object[]{Integer.valueOf(rv)});
                    }

                    System.out.print(String.format("%7s", new Object[]{valueString}));
                    System.out.print("  ");
                }

                System.out.println();
            }

            System.out.println();
        }
    }

    public static void display(String str) {
        fprintf("%s%n", new Object[]{str});
    }

    public static void disp(String str) {
        fprintf("%s%n", new Object[]{str});
    }

    public static String sprintf(String format, Object... os) {
        return String.format(format, os);
    }

    public static void fprintf(String format, Object... os) {
        System.out.format(format, os);
    }

    public static void printf(String format, Object... os) {
        System.out.format(format, os);
    }

    public static void print(String content) {
        System.out.print(content);
    }

    public static void print(char c) {
        System.out.print(c);
    }

    public static void print(char[] s) {
        System.out.print(s);
    }

    public static void print(int[] A) {
        int n = A.length;

        for(int i = 0; i < n - 1; ++i) {
            System.out.print(A[i]);
            System.out.print(' ');
        }

        System.out.print(A[n - 1]);
    }

    public static void print(Object obj) {
        System.out.print(obj);
    }

    public static void println(String content) {
        System.out.println(content);
    }

    public static void println(char[] s) {
        System.out.println(s);
    }

    public static void println(int[] A) {
        int n = A.length;

        for(int i = 0; i < n - 1; ++i) {
            System.out.print(A[i]);
            System.out.print(' ');
        }

        System.out.println(A[n - 1]);
    }

    public static void println(Object obj) {
        System.out.println(obj);
    }

    public static void println() {
        System.out.println();
    }

    public static void err(String input) {
        System.err.println(input);
    }

    public static void errf(String format, Object... os) {
        System.err.format(format, os);
    }

    public static void showSpec(DenseVector V, String[] spec) {
        showSpec(V, spec, 4);
    }

    public static void showSpec(DenseVector V, String[] spec, int p) {
        if(V instanceof DenseVector) {
            int dim = V.getDim();
            double[] pr = V.getPr();

            for(int k = 0; k < dim; ++k) {
                print("  ");
                double v = pr[k];
                int rv = (int)Math.round(v);
                String valueString;
                if(v != (double)rv) {
                    valueString = sprintf(sprintf("%%.%df", new Object[]{Integer.valueOf(p)}), new Object[]{Double.valueOf(v)});
                } else {
                    valueString = sprintf("%d", new Object[]{Integer.valueOf(rv)});
                }

                println(sprintf(sprintf("%%%ds  %%s", new Object[]{Integer.valueOf(8 + p - 4)}), new Object[]{valueString, spec[k]}));
            }

            println();
        } else {
            System.err.println("The input vector should be a DenseVector instance");
//            Utility.exit(1);
        }

    }
}
