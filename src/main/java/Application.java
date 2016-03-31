import decomposition.CholesskyDecomposition;
import decomposition.LUDecomposition;
import decomposition.QRDecomposition;
import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.MTXReader;
import utils.Matlab;
import utils.Portrait;
import utils.Printer;
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
        try {


            String filename = "bcsstk08";
            reader.read("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\matrixes\\"+filename+".mtx");
            SparseMatrix matrix = reader.getMatrix();

            SparseMatrix A = (SparseMatrix)matrix.getSubMatrix(0,matrix.getRowDimension(),0,matrix.getColumnDimension()-150);
//            Portrait port = new Portrait(A,filename);


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

            SparseMatrix AtA = (SparseMatrix)A.transpose().mtimes(A);
            SparseVector Atb = (SparseVector)A.transpose().operate(b);

            start = System.currentTimeMillis();
            CholesskyDecomposition chd = new CholesskyDecomposition(AtA);
            Vector chX = chd.solve(Atb);
            Printer.fprintf("Разложение Холецкого: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecomposition lud = new LUDecomposition(AtA);
            Vector luX = lud.solve(Atb);
            Printer.fprintf("LUP-разложение: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            System.out.println("------------------------------");
            System.out.println("Размерность матрицы: " + A.getRowDimension() + "x" + A.getColumnDimension());
            System.out.println("Количество ненулевых элементов матрицы A: " + A.getNNZ() + "; Общее количество элементов: " + A.getColumnDimension()*A.getRowDimension());
            System.out.println("Число обусловленности матрицы A: cond2 = " + A.cond2());
            System.out.println("Количество ненулевых элементов матрицы AtA: " + AtA.getNNZ() + "; Общее количество элементов: " + AtA.getColumnDimension()*AtA.getRowDimension());
            System.out.println("Число обусловленности матрицы AtA: cond2 = " + AtA.cond2());
            System.out.println("------------------------------");
            System.out.println("QR-разложение");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(qrX)).norm2()/XetVector.norm2());
            System.out.println("------------------------------");
            System.out.println("Разложение Холецкого (матрица AtA)");
            System.out.println("Относительная погрешность: " + Atb.minus(AtA.operate(chX)).norm2()/XetVector.norm2());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (матрица AtA)");
            System.out.println("Относительная погрешность: " + Atb.minus(AtA.operate(luX)).norm2()/XetVector.norm2());

        }catch (IOException e){
            System.out.println("Error");
        }
    }

}
