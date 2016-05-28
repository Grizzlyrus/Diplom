import decomposition.*;
import matrix.DenseMatrix;
import matrix.Matrix;
import matrix.SparseMatrix;
import utils.*;
import vector.DenseVector;
import vector.SparseVector;
import vector.Vector;

import java.io.IOException;
import java.math.BigDecimal;

/**
 * Created by Admin on 09.03.2016.
 */
public class Application {
    public static void main(String args[]){
        MTXReader reader = new MTXReader();
        MTXWriter writer = new MTXWriter();
        try {

//            TODO чтение и зпись матриц

            String filename = "ash958-1";
            reader.read("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\matrixes\\"+filename+".mtx");
            SparseMatrix matrix = reader.getMatrix();


//            SparseMatrix A = (SparseMatrix)matrix.getSubMatrix(0,matrix.getRowDimension(),0,matrix.getColumnDimension()-250);
            SparseMatrix A = (SparseMatrix)matrix.copy();
//            DenseMatrix A = Matlab.full(matrix);
//            writer.write("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\matrixes\\"+filename+"-1.mtx",A);
            Portrait port = new Portrait(A,filename);

//            TODO нормальная и расширенная системы

            SparseMatrix Arash = new SparseMatrix(A.getRowDimension()+A.getColumnDimension(), A.getRowDimension()+A.getColumnDimension());
            for(int i = 0; i<A.getRowDimension();i++){
                Arash.setEntry(i,i,1.0D);
            }

            for(int i = 0; i<A.getRowDimension();i++){
                for(int j= A.getRowDimension();j<A.getRowDimension()+A.getColumnDimension();j++){
                    Arash.setEntry(i,j,A.getEntry(i,j-A.getRowDimension()));
                }
            }
            for(int i = A.getRowDimension(); i<A.getRowDimension()+A.getColumnDimension();i++){
                for(int j= 0;j<A.getRowDimension();j++){
                    Arash.setEntry(i,j,A.getEntry(j,i-A.getRowDimension()));
                }
            }
            Portrait Arashport = new Portrait(Arash, "Arash"+filename);

            double[] Xet = new double[A.getColumnDimension()];
            for (int i = 0; i < A.getColumnDimension(); i++) {
                Xet[i]=i+1;
            }

            DenseVector XetVector = new DenseVector(Xet);
            Vector b = A.operate(XetVector);
            DenseVector brash = new DenseVector(A.getRowDimension()+A.getColumnDimension());
            for (int i = 0; i < b.getDim(); i++) {
                brash.set(i,b.get(i));
            }

            SparseMatrix AtA = (SparseMatrix)A.transpose().mtimes(A);
            Portrait Aport = new Portrait(AtA,"AtA"+filename);
            DenseVector Atb = (DenseVector)A.transpose().operate(b);

//            TODO разложения

            long start = System.currentTimeMillis();
            QRDecomposition qrd = new QRDecomposition(A);
            Vector qrX = qrd.solve(b);
            Printer.fprintf("QR-разложение: %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            QRDecomposition qrdata = new QRDecomposition(AtA);
            Vector qrXata = qrdata.solve(Atb);
            Printer.fprintf("QR-разложение (м-ца AtA): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            QRDecomposition qrdrash = new QRDecomposition(Arash);
            Vector qrXrash = qrdrash.solve(brash);
            Printer.fprintf("QR-разложение (расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            CholesskyDecomposition chd = new CholesskyDecomposition(AtA);
            Vector chX = chd.solve(Atb);
            Printer.fprintf("Разложение Холецкого (м-ца AtA): %.10f seconds.\n", new Object[]{Float.valueOf((float) (System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            CholesskyDecomposition chdrash = new CholesskyDecomposition(Arash);
            Vector chXrash = chdrash.solve(brash);
            Printer.fprintf("Разложение Холецкого (расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float) (System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionEx ludex = new LUDecompositionEx(AtA);
            Vector luXex = ludex.solve(Atb);
            Printer.fprintf("LUP-разложение (м-ца AtA): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionEx ludexrash = new LUDecompositionEx(Arash);
            Vector luXexrash = ludexrash.solve(brash);
            Printer.fprintf("LUP-разложение (расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});


            start = System.currentTimeMillis();
            LUDecomposition lud = new LUDecomposition(AtA);
            Vector luX = lud.solve(Atb);
            Printer.fprintf("LUP-разложение (выбор ведущего элемента по столбцу, м-ца AtA): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecomposition ludrash = new LUDecomposition(Arash);
            Vector luXrash = ludrash.solve(brash);
            Printer.fprintf("LUP-разложение (выбор ведущего элемента по столбцу, расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionMark ludM = new LUDecompositionMark(AtA);
            Vector luXM = ludM.solve(Atb);
            Printer.fprintf("LUP-разложение (стратегия Марковица): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionMark ludMrash = new LUDecompositionMark(Arash);
            Vector luXMrash = ludMrash.solve(brash);
            Printer.fprintf("LUP-разложение (стратегия Марковица, расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionTest ludt = new LUDecompositionTest(AtA);
            Vector luXT = ludt.solve(Atb);
            Printer.fprintf("LUP-разложение (выбор ведущего элемента по активной подматрице): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            start = System.currentTimeMillis();
            LUDecompositionTest ludtrash = new LUDecompositionTest(Arash);
            Vector luXTrash = ludtrash.solve(brash);
            Printer.fprintf("LUP-разложение (выбор ведущего элемента по активной подматрице, расширенная система): %.10f seconds.\n", new Object[]{Float.valueOf((float)(System.currentTimeMillis() - start) / 1000.0F)});

            System.out.println("------------------------------");
            System.out.println("Размерность матрицы: " + A.getRowDimension() + "x" + A.getColumnDimension());
            System.out.println("Количество ненулевых элементов матрицы A: " + A.getNNZ());
            System.out.println("Число обусловленности матрицы A: cond2 = " + A.cond2());
            System.out.println("Количество ненулевых элементов матрицы AtA: " + AtA.getNNZ());
//            System.out.println("Число обусловленности матрицы AtA: cond2 = " + AtA.cond2());
            System.out.println("------------------------------");
            System.out.println("QR-разложение");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(qrX)).normInf());
            System.out.println("Количество ненулевых элементов в матрице Q: " + Matlab.sparse(qrd.getQ()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице R: " + Matlab.sparse(qrd.getR()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("QR-разложение (матрица AtA)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(qrXata)).normInf());
            System.out.println("Количество ненулевых элементов в матрице Q: " + Matlab.sparse(qrdata.getQ()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице R: " + Matlab.sparse(qrdata.getR()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("QR-разложение (расширенная система)");
            DenseVector sol = new DenseVector(A.getRowDimension());
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i, qrXrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице Q: " + Matlab.sparse(qrdrash.getQ()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице R: " + Matlab.sparse(qrdrash.getR()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("Разложение Холецкого (матрица AtA)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(chX)).normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(chd.getL()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("Разложение Холецкого (расширенная система)");
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i,chXrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(chdrash.getL()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (матрица AtA) (выбор ведущего элемента по столбцу)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(luX)).normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(lud.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(lud.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (расширенная система) (выбор ведущего элемента по столбцу)");
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i, luXrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludrash.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludrash.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (стратегия марковица) (матрица AtA)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(luXM)).normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludM.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludM.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (стратегия марковица) (расширенная система)");
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i,luXMrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludMrash.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludMrash.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (выбо ведущего эл-та по активной подматрице) (матрица AtA)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(luXT)).normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludt.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludt.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение (выбо ведущего эл-та по активной подматрице) (расширенная система)");
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i,luXTrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludtrash.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludtrash.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение  (матрица AtA)");
            System.out.println("Относительная погрешность: " + b.minus(A.operate(luXex)).normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludex.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludex.getU()).getNNZ());
            System.out.println("------------------------------");
            System.out.println("LUP-разложение  (расширенная система)");
            for (int i = 0; i < sol.getDim(); i++) {
                sol.set(i,luXexrash.get(i));
            }
            System.out.println("Относительная погрешность: " + sol.normInf());
            System.out.println("Количество ненулевых элементов в матрице L: " + Matlab.sparse(ludexrash.getL()).getNNZ());
            System.out.println("Количество ненулевых элементов в матрице U: " + Matlab.sparse(ludexrash.getU()).getNNZ());
        }catch (IOException e){
            System.out.println("Error");
        }
    }

}
