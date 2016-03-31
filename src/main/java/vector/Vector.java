package vector;

import matrix.Matrix;

/**
 * Created by Admin on 08.03.2016.
 */
public interface Vector {
    int getDim();

    Vector copy();

    Vector times(Vector var1);

    Vector times(double var1);

    Vector plus(Vector var1);

    Vector minus(Vector var1);

    double get(int var1);

    void set(int var1, double var2);

    Vector operate(Matrix var1);

    void clear();

    double norm1();

    double norm2();

    double normInf();
}

