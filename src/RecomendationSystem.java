import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Scanner;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

//Cambiado para proyecto final. 

public class RecomendationSystem {

  private static Scanner scanner;


  /**
   * Compute a dot product between a row c of UkSquareRootSk and a column p of SquareRootSkVk
   * 
   * @param UkSquareRootSk
   * @param SquareRootSkVk
   * @param c costumer selected
   * @param p poduct selected
   * @return dotProduct
   */
  private static double cpDotProduct(Matrix UkSquareRootSk, Matrix SquareRootSkVk, int c, int p) {
    double dotProduct = 0;
    for (int j = 0; j < UkSquareRootSk.getColumnDimension(); j++) {
      dotProduct += UkSquareRootSk.get(c, j) * SquareRootSkVk.get(j, p);
    }

    return dotProduct;
  }

  public static double dotProductWithMatrix(Matrix vectorA, Matrix vectorB) {
    double dotProduct = 0;

    for (int i = 0; i < vectorA.getColumnDimension(); i++) {
      dotProduct += vectorA.get(0, i) * vectorB.get(0, i);
    }

    return dotProduct;

  }

  /**
   * @param R: Matrix to calculate the average per column
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
      }
    }
    return columnAverage;
  }

  // Calculate k to reduce USV matrices
  public static int getK(Matrix S) {

    double sumaDiagonalCuadrado = 0;
    for (int ij = 0; ij < S.rank(); ij++) {
      sumaDiagonalCuadrado += S.get(ij, ij) * S.get(ij, ij);
    }

    double disminuir = sumaDiagonalCuadrado;
    int k;

    for (k = S.rank() - 1; k >= 0; k--) {
      // System.out.println(disminuir / sumaDiagonalCuadrado);
      if (disminuir / sumaDiagonalCuadrado <= 0.9) {

        break;
      } else {
        disminuir -= S.get(k, k) * S.get(k, k);
      }
    }

    return k;
  }


  public static Matrix getMoviesRated(Matrix backUpAdjacencyMatrix, Matrix top5Customers) {
    Matrix moviesRated = new Matrix(2, backUpAdjacencyMatrix.getColumnDimension());

    // assign index to moviesRated
    for (int i = 0; i < moviesRated.getColumnDimension(); i++) {
      moviesRated.set(0, i, i); // assgin index
      moviesRated.set(1, i, 0); // initial values

    }
    // customer i
    for (int i = 0; i < top5Customers.getColumnDimension(); i++) {

      int customerIndex = (int) top5Customers.get(0, i);

      for (int j = 0; j < backUpAdjacencyMatrix.getColumnDimension(); j++) {

        if (backUpAdjacencyMatrix.get(customerIndex, j) == 1.0) {
          double anterior = moviesRated.get(1, j);
          double siguiente = anterior + 1;
          moviesRated.set(1, j, siguiente); // add 1 to movie
        }
      }

    }

    return moviesRated;

  }

  /**
   * @param R: Matrix to calculate the average per row
   * @return Array of average
   */
  private static ArrayList<Double> getRowAverage(Matrix R) {
    ArrayList<Double> rowAverage = new ArrayList<Double>(R.getRowDimension());
    for (int i = 0; i < R.getRowDimension(); i++) {
      Double average = 0.0;
      for (int j = 0; j < R.getColumnDimension(); j++) {
        if (j + 1 == R.getColumnDimension()) {
          average += R.get(i, j);
          rowAverage.add(average / R.getColumnDimension());
        }
        average += R.get(i, j);
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

    // Create Matrix of FileReader
    Matrix A = Matrix.read(b);
    Matrix Atopn = A.copy();

    System.out.println("Matrix of Custumers i - Products j = ");
    // Print matrix A, first parameter is the width(for better read), second is the number of digits after 0
    A.print(4, 2);
    // Calculate the recommendation to customer c and product p
    int c = 0;
    int p = 0;

    scanner = new Scanner(System.in);

    System.out.println("Digite al cliente que deseas realizar la prediccion y el top-5: ");
    c = scanner.nextInt();
    System.out.println("Digite el producto (para la recomendacion): ");
    p = scanner.nextInt();


    // Get the column average of A
    ArrayList<Double> columnAverage = new ArrayList<Double>(A.getColumnDimension());

    // Get average per column
    columnAverage = getColumnAverage(A);

    // Assign the average to the Matrix A if the column is 0
    for (int j = 0; j < A.getColumnDimension(); j++) {
      for (int i = 0; i < A.getRowDimension(); i++) {
        if (A.get(i, j) == 0) {
          A.set(i, j, columnAverage.get(j));
        }
      }
    }

    // System.out.println("Column Average R = ");
    // A.print(4, 2);

    // Calculate Row Average
    ArrayList<Double> rowAverage = new ArrayList<Double>(A.getRowDimension());
    // Get average per row
    rowAverage = getRowAverage(A);

    // To each Aij assign (Aij-rowAverage)
    for (int i = 0; i < A.getRowDimension(); i++) {
      for (int j = 0; j < A.getColumnDimension(); j++) {
        A.set(i, j, A.get(i, j) - rowAverage.get(i));
      }
    }

    // System.out.println("Row Average R = ");
    // A.print(4, 2);

    // Compute the singular value decomposition to get A-norm
    System.out.println("A = U S V^T");
    System.out.println();
    // Get Singular Value of A
    SingularValueDecomposition s = A.svd();
    // Get U Matrix
    System.out.println("U = ");
    Matrix U = s.getU();
    U.print(4, 2);
    // Get S Matrix
    System.out.println("S = ");
    Matrix S = s.getS();
    S.print(4, 2);
    // Get V Matrix
    System.out.println("V = ");
    Matrix V = s.getV();
    V.print(4, 2);
    System.out.println("Effective numerical matrix rank = " + s.rank());
    System.out.println("Two norm condition number = " + s.cond());
    System.out.println("Two norm = " + s.norm2());

    // Calculate singular values
    System.out.println("singular values = ");
    Matrix svalues = new Matrix(s.getSingularValues(), 1);
    svalues.print(4, 2);

    // get k to reduce USV
    int k = getK(S);
    System.out.println("k = " + k);

    // Reduce USV with k to obtain Uk Sk Vk
    Matrix Uk = U.getMatrix(0, U.getRowDimension() - 1, 0, k);
    Matrix Sk = S.getMatrix(0, k, 0, k);
    Matrix Vk = V.getMatrix(0, k, 0, V.getColumnDimension() - 1);

    System.out.println("Uk: ");
    Uk.print(4, 2);
    System.out.println("Sk: ");
    Sk.print(4, 2);
    System.out.println("Vk: ");
    Vk.print(4, 2);


    // Compute square-root of Sk
    System.out.println("Sk^(1/2): ");
    Matrix squareRootSk = newton(Sk, Matrix.identity(Sk.getRowDimension(), Sk.getColumnDimension()));
    squareRootSk.print(4, 2);

    // Compute resultant matrices
    System.out.println("UkSk^(1/2): ");
    Matrix UkSquareRootSk = Uk.times(squareRootSk);
    UkSquareRootSk.print(4, 2);
    System.out.println("Sk^(1/2)Vk: ");
    Matrix SquareRootSkVk = squareRootSk.times(Vk);
    SquareRootSkVk.print(4, 2);

    // Calculate the dotProduct
    double dotProduct = cpDotProduct(UkSquareRootSk, SquareRootSkVk, c, p);
    double customerRatingAverge = rowAverage.get(c);
    // Get CP recommendation
    double recommendation = customerRatingAverge + dotProduct;
    System.out.println("Dot Product: " + dotProduct);
    System.out.println("customerRatingAverge: " + customerRatingAverge);

    // Recommendation TOP-N Products for client c
    System.out.println("A Top N:");
    Atopn.print(4, 2);


    // Change positive values to 1, focusing on consumption, instead of rating
    for (int i = 0; i < Atopn.getRowDimension(); i++) {
      for (int j = 0; j < Atopn.getColumnDimension(); j++) {
        if (Atopn.get(i, j) > 0) {
          Atopn.set(i, j, 1);
        }
      }
    }

    System.out.println("Positive values ​​changed to 1");
    Atopn.print(4, 2);

    Matrix backUpAdjacencyMatrix = Atopn.copy();


    // Get average per column
    ArrayList<Double> columnAverageTopn = getColumnAverage(Atopn);
    // Assign the average if column is 0
    for (int j = 0; j < Atopn.getColumnDimension(); j++) {
      for (int i = 0; i < Atopn.getRowDimension(); i++) {
        if (Atopn.get(i, j) == 0) {
          Atopn.set(i, j, columnAverageTopn.get(j));
        }
      }
    }

    System.out.println("Column Average R Top N = ");
    Atopn.print(4, 2);

    // Calculate Row Average
    ArrayList<Double> rowAverageTopn = new ArrayList<Double>(Atopn.getRowDimension());
    // Get average per column
    rowAverageTopn = getRowAverage(Atopn);

    // To each Rij assign (Rij-rowAverage)
    for (int i = 0; i < Atopn.getRowDimension(); i++) {
      for (int j = 0; j < Atopn.getColumnDimension(); j++) {
        Atopn.set(i, j, Atopn.get(i, j) - rowAverageTopn.get(i));
      }
    }

    System.out.println("Row Average R Top N = ");
    Atopn.print(4, 2);


    System.out.println("Atopn = Utopn Stopn Vtopn^T");
    System.out.println();

    // Print matrix Atopn, first parameter is the width(for better read), second is the number of digits afer 0
    Atopn.print(4, 2);


    SingularValueDecomposition stopn = Atopn.svd();

    System.out.println("U = ");
    // Get U Matrix
    Matrix Utopn = stopn.getU();
    Utopn.print(4, 2);

    System.out.println("Sigma = ");
    // Get Sigma Matrix
    Matrix Stopn = stopn.getS();
    S.print(4, 2);


    // get ktopn to reduce UtopnStopnVtopn
    int ktopn = getK(Stopn);
    System.out.println("k = " + k);

    // Reduce Uktopn Sktopn with ktopn to obtain Uktopn Sktopn
    Matrix Uktopn = Utopn.getMatrix(0, Utopn.getRowDimension() - 1, 0, ktopn);
    Matrix Sktopn = Stopn.getMatrix(0, ktopn, 0, ktopn);

    System.out.println("Uktopn: ");
    Uktopn.print(4, 2);
    System.out.println("Sktopn: ");
    Sktopn.print(4, 2);


    // Compute square-root of Sktopn
    System.out.println("Sktopn^(1/2): ");
    Matrix squareRootSktopn = newton(Sktopn, Matrix.identity(Sktopn.getRowDimension(), Sktopn.getColumnDimension()));
    squareRootSktopn.print(4, 2);

    // Compute resultant matrices
    System.out.println("UkSktopn^(1/2): ");
    Matrix UkSquareRootSktopn = Uktopn.times(squareRootSktopn);
    UkSquareRootSktopn.print(4, 2);

    // Unsorted customers relations
    System.out.println("Unsorted customers relations: ");
    Matrix relationCwithAll = similarCustomers(UkSquareRootSktopn, c);
    relationCwithAll.print(4, 2);

    // Sorted customers relations
    System.out.println("Sorted customers relations: ");
    sortMatrix2xN(relationCwithAll);
    relationCwithAll.print(4, 2);

    // Top-5-Customers
    System.out.println("Top 5 Customers: ");
    Matrix top5Customers = relationCwithAll.getMatrix(0, 0, 1, 5);
    top5Customers.print(4, 2);


    // Adjacency List Customer vs Product
    System.out.println("Adjacency List: ");
    backUpAdjacencyMatrix.print(4, 2);

    // Top-5 products
    // 2xn Matrix
    // 0 1 2 3 4 ... j
    // r0 r1 r2 r3 r4 ... rj
    // ri is times a movie is rated
    System.out.println("Movies Rated: ");
    Matrix moviesRated = getMoviesRated(backUpAdjacencyMatrix, top5Customers);
    moviesRated.print(4, 2);

    // Sorted movies Rated
    System.out.println("Sort Movies Rated: ");
    Matrix moviesRatedSorted = sortMatrix2xN(moviesRated);
    moviesRatedSorted.print(4, 2);

    // Recomendation
    DecimalFormat df = new DecimalFormat("#.00");

    System.out
      .println("Recomendation rating for client " + c + " of product " + p + " is " + df.format((recommendation)));
    System.out.println();

    // Top-5 Movies
    System.out.println("Top-5 Movies recommended to client " + c + ": ");
    Matrix top5movies = moviesRatedSorted.getMatrix(0, 0, 0, 4);
    top5movies.print(2, 0);


  }

  /**
   * Calculate square root of a Matrix
   * 
   * @param a matrix to calculate square root
   * @param x0 identity matrix of a
   * @return Matrix
   */
  public static Matrix newton(Matrix a, Matrix x0) {
    Matrix x1;
    double error = 0.0;
    do {
      x1 = x0.plus(a.times(x0.inverse())).times(0.5);
      error = Math.abs(x1.normInf() - x0.normInf());
      x0 = x1;
    } while (error > 0.000001);
    return x0;
  }


  // in: Matrix of
  public static Matrix similarCustomers(Matrix M, int c) {
    int columns = M.getColumnDimension();
    Matrix similarClients = new Matrix(2, M.getRowDimension());

    Matrix vectorC = M.getMatrix(c, c, 0, columns - 1);
    for (int i = 0; i < M.getRowDimension(); i++) {
      if (i != c) {
        Matrix vectorCi = M.getMatrix(i, i, 0, columns - 1);
        double dotProductCwithCi = dotProductWithMatrix(vectorC, vectorCi);
        double normProduct = vectorC.norm2() * vectorCi.norm2();
        double simcwithci = dotProductCwithCi / normProduct;
        similarClients.set(0, i, i);
        similarClients.set(1, i, simcwithci);
      } else {
        similarClients.set(0, i, i);
        similarClients.set(1, i, 9); // a different number of [-1, 1]
      }
    }

    return similarClients;
  }


  public static Matrix sortMatrix2xN(Matrix A) {
    for (int j = 1; j < A.getColumnDimension(); j++) {
      double key = A.get(1, j);
      int i = j - 1;
      while (i >= 0 && A.get(1, i) < key) {
        A.set(1, i + 1, A.get(1, i));
        A.set(0, i + 1, A.get(0, i));
        i = i - 1;
      }
      A.set(1, i + 1, key);
      A.set(0, i + 1, j);

    }

    return A;
  }


}