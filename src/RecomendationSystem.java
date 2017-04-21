<<<<<<< HEAD
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Jama.Matrix;
import Jama.SingularValueDecomposition;


public class RecomendationSystem {

  /**
   * @param R
   * @return Array of average
   */
  private static ArrayList<Double> getColumnAverage(Matrix R) {
    // R.get(row,column)
    ArrayList<Double> columnAverage = new ArrayList<Double>(R.getColumnDimension());
    for (int j = 0; j < R.getColumnDimension(); j++) {
      Double average = 0.0;
      for (int i = 0; i < R.getRowDimension(); i++) {
        if (i + 1 == R.getRowDimension()) {
          average += R.get(i, j);
          columnAverage.add(average / R.getRowDimension());
        }
        average += R.get(i, j);
        // System.out.println("Column: " + j + " Row: " + i + " = " + r.get(i, j));
      }
    }
    return columnAverage;
  }

  private static ArrayList<Double> getRowAverage(Matrix R) {
    // R.get(row,column)
    ArrayList<Double> rowAverage = new ArrayList<Double>(R.getRowDimension());
    for (int i = 0; i < R.getRowDimension(); i++) {
      Double average = 0.0;
      for (int j = 0; j < R.getColumnDimension(); j++) {
        if (j + 1 == R.getColumnDimension()) {
          average += R.get(i, j);
          rowAverage.add(average / R.getColumnDimension());
        }
        average += R.get(i, j);
        // System.out.println("Column: " + j + " Row: " + i + " = " + r.get(i, j));
      }
    }
    return rowAverage;
  }

  public static void main(String[] args) throws IOException {
    /*
     * A will be a Matrix that represent the User and the Rating for each product
     * Columns will be the Products
     * Rows will be the Users
     */

    // Read File
    String filePath = new File("").getAbsolutePath();
    FileReader f = new FileReader(filePath + "/input/in.txt");
    BufferedReader b = new BufferedReader(f);
    // Create Matrix of File read
    Matrix A = Matrix.read(b);
    System.out.println("A = ");
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);


    ArrayList<Double> columnAverage = new ArrayList<Double>(A.getColumnDimension());

    // Get average per column
    columnAverage = getColumnAverage(A);

    // Assign the average if column is 0
    for (int j = 0; j < A.getColumnDimension(); j++) {
      for (int i = 0; i < A.getRowDimension(); i++) {
        if (A.get(i, j) == 0) {
          A.set(i, j, columnAverage.get(j));
        }
      }
    }

    System.out.println("Column Average R = ");
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);

    // Calculate Row Average
    ArrayList<Double> rowAverage = new ArrayList<Double>(A.getRowDimension());
    // Get average per column
    rowAverage = getRowAverage(A);

    // To each Rij assign (Rij-rowAverage)
    for (int i = 0; i < A.getRowDimension(); i++) {
      for (int j = 0; j < A.getColumnDimension(); j++) {
        A.set(i, j, A.get(i, j) - rowAverage.get(i));
      }
    }

    System.out.println("Row Average R = ");
    A.print(1, 2);


    // compute the singular value decomposition
    System.out.println("A = U S V^T");
    System.out.println();
    // Get Singular Value of A
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);
    SingularValueDecomposition s = A.svd();
    System.out.println("U = ");
    // Get U Matrix
    Matrix U = s.getU();
    // U.print(1, 2);
    System.out.println("Sigma = ");
    // Get Sigma Matrix
    Matrix S = s.getS();
    Matrix identidad = Matrix.identity(S.getRowDimension(), S.getColumnDimension());
    S.print(1, 2);

    System.out.println("Sigma k = ");
    Matrix Sk = newton(S, identidad);
    Sk.print(1, 2);

    System.out.println("V = ");
    // Get V Matrix
    Matrix V = s.getV();
    V.print(1, 2);
    System.out.println("rank = " + s.rank());
    System.out.println("condition number = " + s.cond());
    System.out.println("2-norm = " + s.norm2());
    // Calculate singular values
    System.out.println("singular values = ");
    Matrix svalues = new Matrix(s.getSingularValues(), 1);
    svalues.print(1, 2);


    // Method
    double sumaDiagonalCuadrado = 0;
    for (int ij = 0; ij < S.rank(); ij++) {
      sumaDiagonalCuadrado += S.get(ij, ij) * S.get(ij, ij);
    }

    double disminuir = sumaDiagonalCuadrado;
    int k;

    for (k = S.rank() - 1; k >= 0; k--) {
      System.out.println(disminuir / sumaDiagonalCuadrado);
      if (disminuir / sumaDiagonalCuadrado <= 0.9) {

        break;
      } else {
        disminuir -= S.get(k, k) * S.get(k, k);
      }
    }

    System.out.println(k);
    for (int i = k + 1; i < S.getColumnDimension(); i++) {
      S.set(i, i, 0.0);
    }

    S.print(1, 2);


    for (int j = k + 1; j < U.getColumnDimension(); j++) {
      for (int i = 0; i < U.getRowDimension(); i++) {
        U.set(i, j, 0.0);
      }
    }
    U.print(1, 2);
    V.print(1, 2);
    for (int i = k + 1; i < V.getRowDimension(); i++) {
      for (int j = 0; j < V.getColumnDimension(); j++) {
        V.set(i, j, 0.0);
      }
    }
    V.print(1, 2);

    Matrix result = newton(S, Matrix.identity(S.getRowDimension(), S.getColumnDimension()));
    result.print(1, 5);
  }


  public static Matrix newton(Matrix a, Matrix x0) {
    int contador = 0;

    Matrix x1;
    double error = 0.0;
    do {
      x1 = x0.plus(a.times(x0.inverse())).times(0.5);
      error = Math.abs(x1.normInf() - x0.normInf());
      x0 = x1;
      System.out.println(contador);
      System.out.println("--->" + error);
      contador++;
    } while (error > 0.000001);
    return x0;
  }

=======
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class RecomendationSystem {

  /**
   * @param R
   * @return Array of average
   */
  private static ArrayList<Double> getColumnAverage(Matrix R) {
    // R.get(row,column)
    ArrayList<Double> columnAverage = new ArrayList<Double>(R.getColumnDimension());
    for (int j = 0; j < R.getColumnDimension(); j++) {
      Double average = 0.0;
      for (int i = 0; i < R.getRowDimension(); i++) {
        if (i + 1 == R.getRowDimension()) {
          average += R.get(i, j);
          columnAverage.add(average / R.getRowDimension());
        }
        average += R.get(i, j);
        // System.out.println("Column: " + j + " Row: " + i + " = " + r.get(i, j));
      }
    }
    return columnAverage;
  }

  private static ArrayList<Double> getRowAverage(Matrix R) {
    // R.get(row,column)
    ArrayList<Double> rowAverage = new ArrayList<Double>(R.getRowDimension());
    for (int i = 0; i < R.getRowDimension(); i++) {
      Double average = 0.0;
      for (int j = 0; j < R.getColumnDimension(); j++) {
        if (j + 1 == R.getColumnDimension()) {
          average += R.get(i, j);
          rowAverage.add(average / R.getColumnDimension());
        }
        average += R.get(i, j);
        // System.out.println("Column: " + j + " Row: " + i + " = " + r.get(i, j));
      }
    }
    return rowAverage;
  }

  public static void main(String[] args) throws IOException {
    /*
     * A will be a Matrix that represent the User and the Rating for each product
     * Columns will be the Products
     * Rows will be the Users
     */

    // Read File
    FileReader f = new FileReader(
      "C:/Users/GAMMAX/git/RecomedationSystem/Matrix de entrada/in.txt");
    BufferedReader b = new BufferedReader(f);
    // Create Matrix of File read
    Matrix A = Matrix.read(b);
    System.out.println("A = ");
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);


    ArrayList<Double> columnAverage = new ArrayList<Double>(A.getColumnDimension());

    // Get average per column
    columnAverage = getColumnAverage(A);

    // Assign the average if column is 0
    for (int j = 0; j < A.getColumnDimension(); j++) {
      for (int i = 0; i < A.getRowDimension(); i++) {
        if (A.get(i, j) == 0) {
          A.set(i, j, columnAverage.get(j));
        }
      }
    }

    System.out.println("Column Average R = ");
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);

    // Calculate Row Average
    ArrayList<Double> rowAverage = new ArrayList<Double>(A.getRowDimension());
    // Get average per column
    rowAverage = getRowAverage(A);

    // To each Rij assign (Rij-rowAverage)
    for (int i = 0; i < A.getRowDimension(); i++) {
      for (int j = 0; j < A.getColumnDimension(); j++) {
        A.set(i, j, A.get(i, j) - rowAverage.get(i));
      }
    }

    System.out.println("Row Average R = ");
    A.print(1, 2);


    // compute the singular value decomposition
    System.out.println("A = U S V^T");
    System.out.println();
    // Get Singular Value of A
    // Print matrix A, first parameter is the width(for better read), second is the number of digits afer 0
    A.print(1, 2);
    SingularValueDecomposition s = A.svd();
    System.out.println("U = ");
    // Get U Matrix
    Matrix U = s.getU();
    U.print(1, 2);
    System.out.println("Sigma = ");
    // Get Sigma Matrix
    Matrix S = s.getS();
    Matrix identidad = Matrix.identity(S.getRowDimension(), S.getColumnDimension());
    S.print(1, 2);

//    System.out.println("Sigma k = ");
//    Matrix Sk = newton(S, identidad);
//    Sk.print(1, 2);

    System.out.println("V = ");
    // Get V Matrix
    Matrix V = s.getV();
    V.print(1, 2);
    System.out.println("rank = " + s.rank());
    System.out.println("condition number = " + s.cond());
    System.out.println("2-norm = " + s.norm2());
    
    // Calculate singular values
    System.out.println("singular values = ");
    Matrix svalues = new Matrix(s.getSingularValues(), 1);
    svalues.print(1, 2);

    //get k to reduce USV
    int k = getK(S);
    System.out.println("k = " + k);
    
    // Reduce USV with k to obtain Uk Sk Vk
    Matrix Uk = U.getMatrix(0, U.getRowDimension()-1, 0, k);
    Matrix Sk = S.getMatrix(0, k, 0, k);
    Matrix Vk = V.getMatrix(0, k, 0, V.getColumnDimension()-1);
    
    System.out.println("Uk: ");
    Uk.print(1, 2);
    System.out.println("Sk: ");
    Sk.print(1, 2);
    System.out.println("Vk: ");
    Vk.print(1, 2);
    

    //Compute square-root of Sk
    System.out.println("Sk^(1/2): ");
    Matrix squareRootSk = newton(Sk, Matrix.identity(Sk.getRowDimension(), Sk.getColumnDimension()));
    squareRootSk.print(1, 2);
    
    //Compute resultant matrices
    System.out.println("UkSk^(1/2): ");
    Matrix UkSquareRootSk = Uk.times(squareRootSk);
    UkSquareRootSk.print(1, 2);
    System.out.println("Sk^(1/2)Vk: ");
    Matrix SquareRootSkVk = squareRootSk.times(Vk);
    SquareRootSkVk.print(1, 2);

    //Compute UkSk^(1/2).Sk^(1/2)Vk (dot product)
    //double dotProduct = UkSquareRootSk.
    
    
//    for (int i = k + 1; i < S.getColumnDimension(); i++) {
//      S.set(i, i, 0.0);
//    }
//
//    S.print(1, 2);
//
//
//    for (int j = k + 1; j < U.getColumnDimension(); j++) {
//      for (int i = 0; i < U.getRowDimension(); i++) {
//        U.set(i, j, 0.0);
//      }
//    }
//    U.print(1, 2);
// //   V.print(1, 2);
//    
//
//    
//    for (int i = k + 1; i < V.getRowDimension(); i++) {
//      for (int j = 0; j < V.getColumnDimension(); j++) {
//        V.set(i, j, 0.0);
//      }
//    }
//    V.print(1, 2);
//
//    //Calculate Sk^(1/2)
//    Matrix result = newton(S, Matrix.identity(S.getRowDimension(), S.getColumnDimension()));
//    result.print(1, 5);
  }


  public static Matrix newton(Matrix a, Matrix x0) {
    int contador = 0;

    Matrix x1;
    double error = 0.0;
    do {
      x1 = x0.plus(a.times(x0.inverse())).times(0.5);
      error = Math.abs(x1.normInf() - x0.normInf());
      x0 = x1;
//      System.out.println(contador);
//      System.out.println("--->" + error);
      contador++;
    } while (error > 0.000001);
    return x0;
  }
  
  //Calculate k to reduce USV matrices
  public static int getK(Matrix S){
	    double sumaDiagonalCuadrado = 0;
	    for (int ij = 0; ij < S.rank(); ij++) {
	      sumaDiagonalCuadrado += S.get(ij, ij) * S.get(ij, ij);
	    }

	    double disminuir = sumaDiagonalCuadrado;
	    int k;

	    for (k = S.rank() - 1; k >= 0; k--) {
//	      System.out.println(disminuir / sumaDiagonalCuadrado);
	      if (disminuir / sumaDiagonalCuadrado <= 0.9) {

	        break;
	      } else {
	        disminuir -= S.get(k, k) * S.get(k, k);
	      }
	    }
	    
	    return k;
  }

>>>>>>> branch 'master' of https://github.com/AndresFVR/RecomedationSystem.git
}