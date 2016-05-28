package utils;

import matrix.Matrix;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * Created by ������ on 01.12.2015.
 */
public class Portrait {

    private static BufferedImage img;
    private int m, n;
    Matrix A;


    public Portrait(Matrix A, String filename) throws IOException {
        this.A = A;
        this.m = A.getRowDimension();
        this.n = A.getColumnDimension();

        img = new BufferedImage(n * 5, m * 5, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) img.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, n * 5, m * 5);
        g.setColor(Color.BLUE);


        for(int i = 0;i<m;i++) {
            for(int j = 0;j<n;j++){
                if(A.getEntry(i, j)!=0){
                    g.fillOval(5 * j + 5, 5 * i + 5, 5, 5);
                }
            }
        }

        ImageIO.write(img, "jpg", new File("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\portraits\\"+filename+".jpg"));
//        JFrame window = new JFrame("Matrix's portrait");
//        window.getContentPane().add(new JScrollPane(new JLabel(new ImageIcon(img))));
//        window.setSize(1024, 768);
//        window.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        window.setVisible(true);
    }

    public void compare(Matrix B, String filename) throws IOException{
        int number = 0;
        if(B.getRowDimension()!=A.getRowDimension() || B.getColumnDimension()!=A.getColumnDimension()) throw new IllegalArgumentException("Matrices must have same dimensions");
        img = new BufferedImage(n * 5, m * 5, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) img.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, n * 5, m * 5);
        g.setColor(Color.BLACK);


        for(int i = 0;i<m;i++) {
            for(int j = 0;j<n;j++){
                if(B.getEntry(i, j)!=0.0D && A.getEntry(i,j)==0.0D){
                    g.setColor(Color.GREEN);
                    g.fillOval(5 * j + 5, 5 * i + 5, 5, 5);
                    number++;
                }else if(B.getEntry(i, j)!=0.0D){
                    g.setColor(Color.BLACK);
                    g.fillOval(5 * j + 5, 5 * i + 5, 5, 5);
                }
            }
        }
        System.out.println(number);
        ImageIO.write(img, "jpg", new File("C:\\Users\\Admin\\IdeaProjects\\Diplom\\src\\main\\resources\\portraits\\"+filename+".jpg"));
    }
}
