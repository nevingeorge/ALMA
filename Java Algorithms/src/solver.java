// modified version of LUDecomposition.java found in the Apache Commons Math 3 package
// performs LU decomposition, solves matrix equations, and computes inverses in the field modulo 2

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonSquareMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.util.FastMath;

public class solver {

    // Entries of LU decomposition.
    private final double[][] lu;
    // Pivot permutation associated with LU decomposition.
    private final int[] pivot;
    // Parity of the permutation associated with the LU decomposition.
    private boolean even;
    // Singularity indicator.
    private boolean singular;

    /*
     * Calculates the LU-decomposition of the given matrix.
     * @param matrix The matrix to decompose.
     * @throws NonSquareMatrixException if matrix is not square
     */
    public solver(RealMatrix matrix) {
    	if (!matrix.isSquare()) {
            throw new NonSquareMatrixException(matrix.getRowDimension(), matrix.getColumnDimension());
    	}

        final int m = matrix.getColumnDimension();
        lu = matrix.getData();
        pivot = new int[m];

        // Initialize permutation array and parity
        for (int row=0; row<m; row++) {
            pivot[row] = row;
        }
        even = true;
        singular = false;

        // Loop over columns
        for (int col=0; col<m; col++) {
            // upper
            for (int row=0; row<col; row++) {
                final double[] luRow = lu[row];
                double sum = luRow[col];
                for (int i=0; i<row; i++) {
                    sum -= luRow[i] * lu[i][col];
                }
                sum = Mod2_MA.mod2(sum);
                luRow[col] = sum;
            }

            // lower
            // permutation row
            int max = col;
            double largest = Double.NEGATIVE_INFINITY;
            for (int row=col; row<m; row++) {
                final double[] luRow = lu[row];
                double sum = luRow[col];
                for (int i=0; i<col; i++) {
                    sum -= luRow[i] * lu[i][col];
                }
                sum = Mod2_MA.mod2(sum);
                luRow[col] = sum;

                // maintain best permutation choice
                if (sum > largest) {
                    largest = FastMath.abs(sum);
                    max = row;
                }
            }

            // Singularity check
            if (lu[max][col] == 0) {
                singular = true;
                return;
            }

            // Pivot if necessary
            if (max != col) {
                double tmp = 0;
                final double[] luMax = lu[max];
                final double[] luCol = lu[col];
                for (int i=0; i<m; i++) {
                    tmp = luMax[i];
                    luMax[i] = luCol[i];
                    luCol[i] = tmp;
                }
                int temp = pivot[max];
                pivot[max] = pivot[col];
                pivot[col] = temp;
                even = !even;
            }
        }
    }

    /*
     * Get a solver for finding the A &times; X = B solution in exact linear sense.
     * @return a solver
     */
    public DecompositionSolver getSolver() {
        return new Solver(lu, pivot, singular);
    }

    // Specialized solver
    private static class Solver implements DecompositionSolver {

    	// Entries of LU decomposition
        private final double[][] lu;

        // Pivot permutation associated with LU decomposition
        private final int[] pivot;

        // Singularity indicator
        private final boolean singular;

        /*
         * Build a solver from decomposed matrix.
         * @param lu entries of LU decomposition
         * @param pivot pivot permutation associated with LU decomposition
         * @param singular singularity indicator
         */
        private Solver(final double[][] lu, final int[] pivot, final boolean singular) {
            this.lu       = lu;
            this.pivot    = pivot;
            this.singular = singular;
        }

        public boolean isNonSingular() {
            return !singular;
        }

        public RealVector solve(RealVector b) {
            final int m = pivot.length;
            if (b.getDimension() != m) {
                throw new DimensionMismatchException(b.getDimension(), m);
            }
            if (singular) {
                throw new SingularMatrixException();
            }

            final double[] bp = new double[m];

            // Apply permutations to b
            for (int row=0; row<m; row++) {
                bp[row] = b.getEntry(pivot[row]);
            }

            // Solve LY = b
            for (int col=0; col<m; col++) {
                final double bpCol = bp[col];
                for (int i=col+1; i<m; i++) {
                    bp[i] = Mod2_MA.mod2(bp[i] - bpCol * lu[i][col]);
                }
            }

            // Solve UX = Y
            for (int col=m-1; col>=0; col--) {
                bp[col] /= lu[col][col];
                final double bpCol = bp[col];
                for (int i=0; i<col; i++) {
                	bp[i] = Mod2_MA.mod2(bp[i] - bpCol * lu[i][col]);
                }
            }

            return new ArrayRealVector(bp, false);
        }

        public RealMatrix solve(RealMatrix b) {
            final int m = pivot.length;
            if (b.getRowDimension() != m) {
                throw new DimensionMismatchException(b.getRowDimension(), m);
            }
            if (singular) {
                throw new SingularMatrixException();
            }

            final int nColB = b.getColumnDimension();

            // Apply permutations to b
            final double[][] bp = new double[m][nColB];
            for (int row=0; row<m; row++) {
                final double[] bpRow = bp[row];
                final int pRow = pivot[row];
                for (int col=0; col<nColB; col++) {
                    bpRow[col] = b.getEntry(pRow, col);
                }
            }

            // Solve LY = b
            for (int col=0; col<m; col++) {
                final double[] bpCol = bp[col];
                for (int i=col+1; i<m; i++) {
                    final double[] bpI = bp[i];
                    final double luICol = lu[i][col];
                    for (int j=0; j<nColB; j++) {
                        bpI[j] = Mod2_MA.mod2(bpI[j] - bpCol[j] * luICol);
                    }
                }
            }

            // Solve UX = Y
            for (int col=m-1; col>=0; col--) {
                final double[] bpCol = bp[col];
                for (int i=0; i<col; i++) {
                    final double[] bpI = bp[i];
                    final double luICol = lu[i][col];
                    for (int j=0; j<nColB; j++) {
                    	bpI[j] = Mod2_MA.mod2(bpI[j] - bpCol[j] * luICol);
                    }
                }
            }

            return new Array2DRowRealMatrix(bp, false);
        }

        /*
         * Get the inverse of the decomposed matrix.
         * @return the inverse matrix.
         * @throws SingularMatrixException if the decomposed matrix is singular.
         */
        public RealMatrix getInverse() {
            return solve(MatrixUtils.createRealIdentityMatrix(pivot.length));
        }
    }

}
