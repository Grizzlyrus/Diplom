package utils;

import matrix.Matrix;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * Created by Кирилл on 01.12.2015.
 */
public class Portrait {

    private static BufferedImage img;
    private int m, n;


    public Portrait(Matrix A, String filename) throws IOException {

        this.m = A.getRowDimension();
        this.n = A.getColumnDimension();

        img = new BufferedImage(m * 10, n * 10, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) img.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, m * 10, n * 10);
        g.setColor(Color.BLACK);


        for(int i = 0;i<m;i++) {
            for(int j = 0;j<n;j++){
                if(A.getEntry(i, j)!=0){
                    g.fillOval(10 * i + 5, 10 * j + 5, 5, 5);
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
}
