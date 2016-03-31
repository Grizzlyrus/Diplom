package utils;

import matrix.DenseMatrix;

/**
 * Created by Admin on 08.03.2016.
 */
public class ArrayOperator {
    public ArrayOperator() {
    }

    public static int[] sort(double[] V, String order) {
        int len = V.length;
        int[] indices = colon(0, 1, len - 1);
        byte start = 0;
        int end = len - 1;
        quickSort(V, (int[])indices, start, end, order);
        return indices;
    }

    public static int[] sort(double[] V) {
        return sort(V, "ascend");
    }

    public static void quickSort(double[] values, int[] indices, int start, int end, String order) {
        int i = start;
        int j = end;
        double temp = values[start];
        int tempV = indices[start];

        do {
            if(order.equals("ascend")) {
                while(values[j] >= temp && j > i) {
                    --j;
                }
            } else if(order.equals("descend")) {
                while(values[j] <= temp && j > i) {
                    --j;
                }
            }

            if(j > i) {
                values[i] = values[j];
                indices[i] = indices[j];
                ++i;
            }

            if(order.equals("ascend")) {
                while(values[i] <= temp && i < j) {
                    ++i;
                }
            } else if(order.equals("descend")) {
                while(values[i] >= temp && i < j) {
                    ++i;
                }
            }

            if(j > i) {
                values[j] = values[i];
                indices[j] = indices[i];
                --j;
            }
        } while(i != j);

        values[i] = temp;
        indices[i] = tempV;
        ++i;
        --j;
        if(start < j) {
            quickSort(values, indices, start, j, order);
        }

        if(i < end) {
            quickSort(values, indices, i, end, order);
        }

    }

    public static void quickSort(double[] values, double[] indices, int start, int end, String order) {
        int i = start;
        int j = end;
        double temp = values[start];
        double tempV = indices[start];

        do {
            if(order.equals("ascend")) {
                while(values[j] >= temp && j > i) {
                    --j;
                }
            } else if(order.equals("descend")) {
                while(values[j] <= temp && j > i) {
                    --j;
                }
            }

            if(j > i) {
                values[i] = values[j];
                indices[i] = indices[j];
                ++i;
            }

            if(order.equals("ascend")) {
                while(values[i] <= temp && i < j) {
                    ++i;
                }
            } else if(order.equals("descend")) {
                while(values[i] >= temp && i < j) {
                    ++i;
                }
            }

            if(j > i) {
                values[j] = values[i];
                indices[j] = indices[i];
                --j;
            }
        } while(i != j);

        values[i] = temp;
        indices[i] = tempV;
        ++i;
        --j;
        if(start < j) {
            quickSort(values, indices, start, j, order);
        }

        if(i < end) {
            quickSort(values, indices, i, end, order);
        }

    }

    public static int fix(double x) {
        return x > 0.0D?(int)Math.floor(x):(int)Math.ceil(x);
    }

    public static int[] colon(int begin, int d, int end) {
        int m = fix((double)((end - begin) / d));
        if(m < 0) {
            System.err.println("Difference error!");
            System.exit(1);
        }

        int[] res = new int[m + 1];

        for(int i = 0; i <= m; ++i) {
            res[i] = begin + i * d;
        }

        return res;
    }

    public static int[] colon(int begin, int end) {
        return colon(begin, 1, end);
    }

    public static double[] colon(double begin, double d, double end) {
        int m = fix((end - begin) / d);
        if(m < 0) {
            System.err.println("Difference error!");
            System.exit(1);
        }

        double[] res = new double[m + 1];

        for(int i = 0; i <= m; ++i) {
            res[i] = begin + (double)i * d;
        }

        return res;
    }

    public static double[] colon(double begin, double end) {
        return colon(begin, 1.0D, end);
    }

    public static int argmax(double[] V) {
        int maxIdx = 0;
        double maxVal = V[0];

        for(int i = 1; i < V.length; ++i) {
            if(maxVal < V[i]) {
                maxVal = V[i];
                maxIdx = i;
            }
        }

        return maxIdx;
    }

    public static int argmax(double[] V, int begin, int end) {
        int maxIdx = begin;
        double maxVal = V[begin];

        for(int i = begin + 1; i < end; ++i) {
            if(maxVal < V[i]) {
                maxVal = V[i];
                maxIdx = i;
            }
        }

        return maxIdx;
    }

    public static int argmin(double[] V) {
        int maxIdx = 0;
        double maxVal = V[0];

        for(int i = 1; i < V.length; ++i) {
            if(maxVal > V[i]) {
                maxVal = V[i];
                maxIdx = i;
            }
        }

        return maxIdx;
    }

    public static int argmin(double[] V, int begin, int end) {
        int maxIdx = begin;
        double maxVal = V[begin];

        for(int i = begin + 1; i < end; ++i) {
            if(maxVal > V[i]) {
                maxVal = V[i];
                maxIdx = i;
            }
        }

        return maxIdx;
    }

    public static double min(double[] V) {
        double res = V[0];

        for(int i = 1; i < V.length; ++i) {
            if(res > V[i]) {
                res = V[i];
            }
        }

        return res;
    }

    public static double max(double[] V) {
        double res = V[0];

        for(int i = 1; i < V.length; ++i) {
            if(res < V[i]) {
                res = V[i];
            }
        }

        return res;
    }

    public static double min(double[] V, int begin, int end) {
        double res = V[0];

        for(int i = begin + 1; i < end; ++i) {
            if(res > V[i]) {
                res = V[i];
            }
        }

        return res;
    }

    public static double max(double[] V, int begin, int end) {
        double res = V[begin];

        for(int i = begin + 1; i < end; ++i) {
            if(res < V[i]) {
                res = V[i];
            }
        }

        return res;
    }

    public static void assignVector(double[] V, double v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] = v;
        }

    }

    public static void assign(double[] V, double v) {
        assignVector(V, v);
    }

    public static void assignIntegerVector(int[] V, int v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] = v;
        }

    }

    public static void assign(int[] V, int v) {
        assign(V, v);
    }

    public static void clearVector(double[] V) {
        assignVector(V, 0.0D);
    }

    public static void clear(double[] V) {
        clearVector(V);
    }

    public static void clearMatrix(double[][] M) {
        for(int i = 0; i < M.length; ++i) {
            assignVector(M[i], 0.0D);
        }

    }

    public static void clear(double[][] M) {
        clearMatrix(M);
    }

    public static double[] allocate1DArray(int n) {
        return allocateVector(n, 0.0D);
    }

    public static double[] allocate1DArray(int n, double v) {
        return allocateVector(n, v);
    }

    public static double[] allocateVector(int n) {
        return allocateVector(n, 0.0D);
    }

    public static double[] allocateVector(int n, double v) {
        double[] res = new double[n];
        assignVector(res, v);
        return res;
    }

    public static double[][] allocate2DArray(int m, int n, double v) {
        double[][] res = new double[m][];

        for(int i = 0; i < m; ++i) {
            res[i] = new double[n];

            for(int j = 0; j < n; ++j) {
                res[i][j] = v;
            }
        }

        return res;
    }

    public static double[][] allocate2DArray(int m, int n) {
        return allocate2DArray(m, n, 0.0D);
    }

    public static int[] allocateIntegerVector(int n) {
        return allocateIntegerVector(n, 0);
    }

    public static int[] allocateIntegerVector(int n, int v) {
        int[] res = new int[n];
        assignIntegerVector(res, v);
        return res;
    }

    public static double[][] allocateMatrix(int nRows, int nCols) {
        double[][] res = new double[nRows][];

        for(int i = 0; i < nRows; ++i) {
            res[i] = allocateVector(nCols);
        }

        return res;
    }

    public static void divideAssign(double[] V, double v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] /= v;
        }

    }

    public static void divideAssign(double[] V1, double[] V2) {
        for(int i = 0; i < V1.length; ++i) {
            V1[i] /= V2[i];
        }

    }

    public static void divideAssign(double[][] res, double v) {
        for(int i = 0; i < res.length; ++i) {
            double[] row = res[i];

            for(int j = 0; j < row.length; ++j) {
                row[j] /= v;
            }
        }

    }

    public static void timesAssign(double[] V, double v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] *= v;
        }

    }

    public static void timesAssign(double[] V1, double[] V2) {
        for(int i = 0; i < V1.length; ++i) {
            V1[i] *= V2[i];
        }

    }

    public static void timesAssign(double[][] res, double v) {
        for(int i = 0; i < res.length; ++i) {
            double[] row = res[i];

            for(int j = 0; j < row.length; ++j) {
                row[j] *= v;
            }
        }

    }

    public static double sum(double[] V) {
        double res = 0.0D;

        for(int i = 0; i < V.length; ++i) {
            res += V[i];
        }

        return res;
    }

    public static double mean(double[] V) {
        return sum(V) / (double)V.length;
    }

    public static double std(double[] V, int flag) {
        int n = V.length;
        if(n == 1) {
            return 0.0D;
        } else {
            double mean = mean(V);
            double res = 0.0D;
            double[] var11 = V;
            int var10 = V.length;

            for(int var9 = 0; var9 < var10; ++var9) {
                double v = var11[var9];
                double diff = v - mean;
                res += diff * diff;
            }

            if(flag == 0) {
                res /= (double)(n - 1);
            } else if(flag == 1) {
                res /= (double)n;
            }

            res = Math.sqrt(res);
            return res;
        }
    }

    public static void sum2one(double[] V) {
        divideAssign(V, sum(V));
    }

    public static void plusAssign(double[] V, double v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] += v;
        }

    }

    public static void plusAssign(double[] res, double a, double[] V) {
        for(int i = 0; i < res.length; ++i) {
            res[i] += a * V[i];
        }

    }

    public static void plusAssign(double[] V1, double[] V2) {
        for(int i = 0; i < V1.length; ++i) {
            V1[i] += V2[i];
        }

    }

    public static void plusAssign(double[][] res, double v) {
        for(int i = 0; i < res.length; ++i) {
            double[] row = res[i];

            for(int j = 0; j < row.length; ++j) {
                row[j] += v;
            }
        }

    }

    public static void minusAssign(int[] V, int v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] -= v;
        }

    }

    public static void minusAssign(double[] V, double v) {
        for(int i = 0; i < V.length; ++i) {
            V[i] -= v;
        }

    }

    public static void minusAssign(double[] res, double a, double[] V) {
        for(int i = 0; i < res.length; ++i) {
            res[i] -= a * V[i];
        }

    }

    public static void minusAssign(double[] V1, double[] V2) {
        for(int i = 0; i < V1.length; ++i) {
            V1[i] -= V2[i];
        }

    }

    public static void minusAssign(double[][] res, double v) {
        for(int i = 0; i < res.length; ++i) {
            double[] row = res[i];

            for(int j = 0; j < row.length; ++j) {
                row[j] -= v;
            }
        }

    }

    public static void assignVector(double[] V1, double[] V2) {
        System.arraycopy(V2, 0, V1, 0, V1.length);
    }

    public static void assign(double[][] res, double[][] A) {
        for(int i = 0; i < res.length; ++i) {
            assignVector(res[i], A[i]);
        }

    }

    public static void assign(double[][] res, double v) {
        for(int i = 0; i < res.length; ++i) {
            assignVector(res[i], v);
        }

    }

    public static double[] operate(double[][] A, double[] V) {
        double[] res = new double[A.length];
        double s = 0.0D;

        for(int i = 0; i < res.length; ++i) {
            s = 0.0D;
            double[] A_i = A[i];

            for(int j = 0; j < V.length; ++j) {
                s += A_i[j] * V[j];
            }

            res[i] = s;
        }

        return res;
    }

    public static void operate(double[] V1, double[][] A, double[] V2) {
        double s = 0.0D;

        for(int i = 0; i < V1.length; ++i) {
            double[] ARow = A[i];
            s = 0.0D;

            for(int j = 0; j < V2.length; ++j) {
                s += ARow[j] * V2[j];
            }

            V1[i] = s;
        }

    }

    public static double[] operate(double[] V, double[][] A) {
        double[] res = new double[A[0].length];
        double s = 0.0D;

        for(int j = 0; j < res.length; ++j) {
            s = 0.0D;

            for(int i = 0; i < V.length; ++i) {
                s += V[i] * A[i][j];
            }

            res[j] = s;
        }

        return res;
    }

    public static void operate(double[] V1, double[] V2, double[][] A) {
        double s = 0.0D;

        for(int j = 0; j < V1.length; ++j) {
            s = 0.0D;

            for(int i = 0; i < V2.length; ++i) {
                s += V2[i] * A[i][j];
            }

            V1[j] = s;
        }

    }

    public static double innerProduct(double[] V1, double[] V2) {
        if(V1 != null && V2 != null) {
            double res = 0.0D;

            for(int i = 0; i < V1.length; ++i) {
                res += V1[i] * V2[i];
            }

            return res;
        } else {
            return 0.0D;
        }
    }

    public static double innerProduct(double[] V1, double[] V2, int from, int to) {
        if(V1 != null && V2 != null) {
            double res = 0.0D;

            for(int i = from; i < to; ++i) {
                res += V1[i] * V2[i];
            }

            return res;
        } else {
            return 0.0D;
        }
    }

    public static double[] times(double[] V1, double[] V2) {
        double[] res = new double[V1.length];

        for(int i = 0; i < V1.length; ++i) {
            res[i] = V1[i] * V2[i];
        }

        return res;
    }

    public static double[] plus(double[] V1, double[] V2) {
        double[] res = new double[V1.length];

        for(int i = 0; i < V1.length; ++i) {
            res[i] = V1[i] + V2[i];
        }

        return res;
    }

    public static double[] minus(double[] V1, double[] V2) {
        double[] res = new double[V1.length];

        for(int i = 0; i < V1.length; ++i) {
            res[i] = V1[i] - V2[i];
        }

        return res;
    }
}
