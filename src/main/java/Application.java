import decomposition.CholesskyDecomposition;
import decomposition.LUDecomposition;
import decomposition.LUDecompositionMark;
import decomposition.QRDecomposition;
import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.*;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.io.IOException;

/**
 * Created by Admin on 09.03.2016.
 */
public class Application {
    public static void main(String args[]){
        MTXReader reader = new MTXReader();
        MTXWriter writer = new MTXWriter();
        try {


            String filename = "well1033";
            reader.read("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\matrixes\\"+filename+".mtx");
            SparseMatrix matrix = reader.getMatrix();


//            SparseMatrix A = (SparseMatrix)matrix.getSubMatrix(0,matrix.getRowDimension(),0,matrix.getColumnDimension()-250);
//            SparseMatrix A = (SparseMatrix)matrix.copy();
            DenseMatrix A = Matlab.full(matrix);
//            writer.write("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\matrixes\\"+filename+"-1.mtx",A);
            Portrait port = new Portrait(A,filename);


            double[] Xet = new double[A.getColumnDimension()];
            for (int i = 0; i < A.getColumnDimension(); i++) {
                Xet[i]=i+1;
            }

            DenseVector XetVector = new DenseVector(Xet);
            Vector b = A.operate(XetVector);

            long start = System.currentTimeMillis();
            QRDecomposition qrd = new QRDecomposition(A);
            Vector qrX = qrd.solve(b);
            Printer.fprintf("QR-разложение: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

//            SparseMatrix AtA = (SparseMatrix)A.transpose().mtimes(A);
            DenseMatrix AtA = (DenseMatrix)A.transpose().mtimes(A);
            DenseVector Atb = (DenseVector)A.transpose().operate(b);

            start = System.currentTimeMillis();
            CholesskyDecomposition chd = new CholesskyDecomposition(AtA);
//            CholesskyDecomposition chd = new CholesskyDecomposition(A);
            Vector chX = chd.solve(Atb);
//            Vector chX = chd.solve(b);
            Printer.fprintf("Разложение Холецкого: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecomposition lud = new LUDecomposition(AtA);
//            LUDecomposition lud = new LUDecomposition(A);
            Vector luX = lud.solve(Atb);
//            Vector luX = lud.solve(b);
            Printer.fprintf("LUP-разложение: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionMark ludM = new LUDecompositionMark(AtA);
//            LUDecompositionMark ludM = new LUDecompositionMark(A);
            Vector luXM = ludM.solve(Atb);
//            Vector luXM = ludM.solve(b);
            Printer.fprintf("LUP-разложение: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});


            System.out.println("------------------------------");
            System.out.println("Размерность матрицы: " + A.getRowDimension() + "x" + A.getColumnDimension());
//            System.out.println("Количество ненулевых элементов матрицы A: " + A.getNNZ());
            System.out.println("Число обусловленности матрицы A: cond2 = " + A.cond2());
//            System.out.println("Количество ненулевых элементов матрицы AtA: " + AtA.getNNZ());
            System.out.println("Число обусловленности матрицы AtA: cond2 = " + AtA.cond2());
            System.out.println("------------------------------");
            System.out.println("QR-разложение");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(qrX)).norm2()/XetVector.norm2());
            System.out.println("Количество ненулевых элементов в матрице Q: " + Matlab.sparse(qrd.getQ()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице R: " + Matlab.sparse(qrd.getR()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("Разложение Холецкого (матрица AtA)");
            System.out.println("Относительная погрешность: " + Atb.minus(AtA.operate(chX)).norm2()/XetVector.norm2());
//            System.out.println("Относительная погрешность: " + b.minus(A.operate(chX)).norm2()/XetVector.norm2());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(chd.getL()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (матрица AtA)");
            System.out.println("Относительная погрешность: " + Atb.minus(AtA.operate(luX)).norm2()/XetVector.norm2());
//            System.out.println("Относительная погрешность: " + b.minus(A.operate(luX)).norm2()/XetVector.norm2());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(lud.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(lud.getU()).getNNZ());
//            Printer.printVector(luX);
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (стратегия марковица) (матрица AtA)");
            System.out.println("Относительная погрешность: " + Atb.minus(AtA.operate(luXM)).norm2()/XetVector.norm2());
//            System.out.println("Относительная погрешность: " + b.minus(A.operate(luXM)).norm2()/XetVector.norm2());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludM.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludM.getU()).getNNZ());
//            Printer.printVector(luXM);
        }catch (IOException e){
            System.out.println("Error");
        }
    }

}
