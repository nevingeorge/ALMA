import java.util.ArrayList;
import java.util.Arrays;

public class LinearAlgebra {
	
	public static int mod2(int n) {
		if(n%2==0)
			return 0;
		else
			return 1;
	}
	
	public static int determinant(int[][] arr) {
		int result = 0;
		if (arr.length == 1) {
			result = arr[0][0];
			return result;
		}
		if (arr.length == 2) {
			result = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
			return result;
		}
		for (int i = 0; i < arr[0].length; i++) {
			int temp[][] = new int[arr.length - 1][arr[0].length - 1];

			for (int j = 1; j < arr.length; j++) {
				for (int k = 0; k < arr[0].length; k++) {

					if (k < i) {
						temp[j - 1][k] = arr[j][k];
					} else if (k > i) {
						temp[j - 1][k - 1] = arr[j][k];
					}
				}
			}
			result += arr[0][i] * Math.pow(-1, (int) i) * determinant(temp);
		}
		return result;
	}
	
	public static int[][] inverse(int[][] a){
		int n = a.length;
		int x[][] = new int[n][n];
		int b[][] = new int[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i) 
            b[i][i] = 1;
 
        gaussian(a, index);
 
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                    	    -= a[index[j]][i]*b[index[i]][k];
 
        for (int i=0; i<n; ++i) 
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j) 
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k) 
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        
        for(int i=0;i<n;i++) {
        	for(int j=0;j<n;j++)
        		x[i][j] = mod2(x[i][j]);
        }
 
        return x;
	}
	
	public static void gaussian(int a[][], int index[]) 
    {
        int n = index.length;
        int c[] = new int[n];
 
        for (int i=0; i<n; ++i) 
            index[i] = i;
 
        for (int i=0; i<n; ++i) 
        {
        	int c1 = 0;
            for (int j=0; j<n; ++j) 
            {
            	int c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }
 
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
        	int pi1 = 0;
            for (int i=j; i<n; ++i) 
            {
            	int pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
 
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
            	int pj = a[index[i]][j]/a[index[j]][j];
 
                a[index[i]][j] = pj;
 
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
	
	public static int[] rowxnxnMatrixMult(int[] arr1, int[][] arr2) {
		int size = arr1.length;
		int[] out = new int[size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++)
				out[i] += arr1[j]*arr2[j][i];
			out[i] = mod2(out[i]);
		}
		return out;
	}
	
	public static int[][] nxnMatrixMult(int[][] arr1, int[][] arr2){
		int size = arr1.length;
		int[][] out = new int[size][size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++) {
				for(int k=0;k<size;k++)
					out[i][j] += arr1[i][k]*arr2[k][j];
				out[i][j] = mod2(out[i][j]);
			}
		}
		return out;
	}
	
	public static int[] nxnxnx1MatrixMult(int[][] arr1, int[] arr2) {
		int size = arr1.length;
		int[] out = new int[size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++)
				out[i] += arr1[i][j]*arr2[j];
			out[i] = mod2(out[i]);
		}
		return out;
	}
	
	public static int dotProduct(int[] v1, int[] v2) {
		int sum = 0;
		for(int i=0;i<v1.length;i++)
			sum += v1[i]*v2[i];
		return mod2(sum);
	}
	
	public static int[][] rref(int[][] matrix) {
		int[][] rref = new int[matrix.length][];
		for (int i = 0; i < matrix.length; i++)
			rref[i] = Arrays.copyOf(matrix[i], matrix[i].length);

		int r = 0;
		for (int c = 0; c < rref[0].length && r < rref.length; c++) {
			int j = r;
			for (int i = r + 1; i < rref.length; i++)
				if (rref[i][c] > rref[j][c])
					j = i;
			if (rref[j][c] == 0)
				continue;

			int[] temp = rref[j];
			rref[j] = rref[r];
			rref[r] = temp;

			for (int i = 0; i < rref.length; i++) {
				if (i != r) {
					int t = rref[i][c];
					for (j = 0; j < rref[0].length; j++)
						rref[i][j] = mod2(rref[i][j] - (t * rref[r][j]));
				}
			}
			r++;
		}

		return rref;
	}
	
	public static boolean linInd(int[] w, ArrayList<int[]> B) {
		// forms augmented matrix B|w
		int numRows = w.length;
		int numCols = B.size()+1;
		int[][] arr = new int[numRows][numCols];
		for(int j=0;j<numCols-1;j++) {
			for(int i=0;i<numRows;i++)
				arr[i][j] = B.get(j)[i];
		}
		for(int i=0;i<numRows;i++)
			arr[i][numCols-1] = w[i];
		
		int[][] rrefArr = rref(arr);
		
		// finds the index of the last 1 in the last column (if exists)
		int index = -1;
		for(int i=numRows-1;i>=0;i--) {
			if(rrefArr[i][numCols-1] == 1) {
				index = i;
				break;
			}
		}
		
		// last vector is the 0 vector, in span(B)
		if(index == -1)
			return false;
		
		// checks whether in span
		for(int j=0;j<numCols-1;j++) {
			if(rrefArr[index][j] == 1)
				return false;
		}
		
		// linearly independent
		return true;
	}
}
