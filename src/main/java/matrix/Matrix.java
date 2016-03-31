package matrix;

import vector.Vector;

/**
 * Created by Admin on 08.03.2016.
 */
public interface Matrix {
    int getRowDimension();

    int getColumnDimension();

    double getEntry(int var1, int var2);

    void setEntry(int var1, int var2, double var3);

    Matrix getSubMatrix(int var1, int var2, int var3, int var4);

    Matrix getSubMatrix(int[] var1, int[] var2);

    Matrix getRows(int var1, int var2);

    Matrix getRows(int... var1);

    Vector[] getRowVectors(int var1, int var2);

    Vector[] getRowVectors(int... var1);

    Matrix getRowMatrix(int var1);

    void setRowMatrix(int var1, Matrix var2);

    Vector getRowVector(int var1);

    void setRowVector(int var1, Vector var2);

    Matrix getColumns(int var1, int var2);

    Matrix getColumns(int... var1);

    Vector[] getColumnVectors(int var1, int var2);

    Vector[] getColumnVectors(int... var1);

    Matrix getColumnMatrix(int var1);

    void setColumnMatrix(int var1, Matrix var2);

    Vector getColumnVector(int var1);

    void setColumnVector(int var1, Vector var2);

    Matrix mtimes(Matrix var1);

    Matrix times(Matrix var1);

    Matrix times(double var1);

    Matrix plus(Matrix var1);

    Matrix plus(double var1);

    Matrix minus(Matrix var1);

    Matrix minus(double var1);

    Matrix transpose();

    Matrix copy();

    Vector operate(Vector var1);

    void clear();

    double[][] getData();

    double norm1();

    double norm2();

    double normE();

    double normInf();

    double cond1();

    double cond2();

    double condE();

    double condInf();

}
