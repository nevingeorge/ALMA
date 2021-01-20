/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program contains the membership query functions that are used in arbitrary.java. Each
 * function takes in as input a word in the specified language and returns as output either 0 or 1.
 */

public class MQ {
	
	// MQ #1
	public static int MQ1(String w) {
		String[] wArr = w.split(" ");
		if (w.length() == 0) {
			wArr = new String[0];
		}
		
		if (wArr.length%2 == 0) {
			return 1;
		}
		return 0;
	}
	
	// MQ #2
	public static int MQ2(String w) {
		String[] wArr = w.split(" ");
		if (w.length() == 0) {
			wArr = new String[0];
		}
		
		if (wArr.length%10 == 0) {
			return 1;
		}
		return 0;
	}
}
