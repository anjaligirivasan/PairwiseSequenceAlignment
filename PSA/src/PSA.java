import java.io.*;
import java.util.*;
import java.lang.*;

public class PSA 
{
	// final arrays containing aligned sequences
	private static int qEnd[], sEnd[];
	
	// variables needed for calculating distance score
	private static int indelCount = 0;
	private static int matchCount = 0;
	private static int totalWOCol = 0;
	
	/**
	 * first part of algorithm based off of Needleman-Wunsch algorithm for calculating
	 * progress pairwise alignment of two sequences
	 * - fills in table -
	 * (altered slightly since negative values in arrays cause issues)
	 * @param pa : 2D array holding table for the algorithm
	 * @param query : first sequence inputted
	 * @param subject : second sequence inputted
	 * @param paNW : indicator of match/mismatch
	 * @param gap : indicator of indel
	 */
	static void needlemanWunsch(int[][] pa, String query, String subject, int paNW, int gap)
	{
		// initializing variables for for loop
		int i = 0;
		int j = 0;

		// length of query
		int qLength = query.length();
		// length of subject
		int sLength = subject.length(); 

		// initializing the table to 2,4,6,... in the first row and column
		for (i = 0; i <= (sLength + qLength); i++)
		{
			pa[i][0] = i * gap;
			pa[0][i] = i * gap;
		}

		// calculating the minimum penalty (filling up the table)
		for (i = 1; i <= qLength; i++)
		{
			for (j = 1; j <= sLength; j++)
			{
				if (query.charAt(i - 1) == subject.charAt(j - 1))
				{
					pa[i][j] = pa[i - 1][j - 1];
				}
				else
				{
					// choosing between mismatch nw, north or west
					pa[i][j] = Math.min(Math.min(pa[i - 1][j - 1] + paNW , pa[i - 1][j] + gap) , pa[i][j - 1] + gap );
				}
			}
		}
	}
	
	/**
	 * Second part of Needleman-Wunsch algorithm that extracts the scores
	 * from the table, inserts indels accordingly and prints aligned sequences
	 * - extracts from table -
	 * @param query
	 * @param subject
	 * @param paNW
	 * @param gap
	 */
	static void findAlignment(String query, String subject, int paNW, int gap)
	{
		int i, j; // initializing variables

		// length of query and subject
		int qLength = query.length();
		int sLength = subject.length();

		// create table to import created one for needleman-wunsch
		int pwa[][] = new int[sLength + qLength + 1][sLength + qLength + 1];
		needlemanWunsch(pwa, query, subject, paNW, gap);

		// maximum possible length (no optimal matches)
		int l = sLength + qLength; 

		//setting variables for loops
		i = qLength; 
		j = sLength;

		int qPosition = l;
		int sPosition = l;

		// ending results for optimal alignment of strings
		qEnd = new int[l + 1];
		sEnd = new int[l + 1];
	
		// Needleman Wunsch Algorithm - extracting characters from tables
		// based on optimal score traversal
		while ( !(i == 0 || j == 0))
		{
			// Case 1: if minimum is they match keep characters
			if (query.charAt(i - 1) == subject.charAt(j - 1))
			{
				qEnd[qPosition--] = (int)query.charAt(i - 1);
				sEnd[sPosition--] = (int)subject.charAt(j - 1);
				i--; 
				j--;
			}
			// Case 2: if minimum is they mismatch keep characters
			else if (pwa[i - 1][j - 1] + paNW == pwa[i][j])
			{
				qEnd[qPosition--] = (int)query.charAt(i - 1);
				sEnd[sPosition--] = (int)subject.charAt(j - 1);
				i--; 
				j--;
			}
			// Case 3: if minimum is west gap insert indel into subject
			else if (pwa[i - 1][j] + gap == pwa[i][j])
			{
				qEnd[qPosition--] = (int)query.charAt(i - 1);
				sEnd[sPosition--] = (int)'_';
				i--;
			}
			// Case 4: if minimum is north gap insert indel into query
			else if (pwa[i][j - 1] + gap == pwa[i][j])
			{
				qEnd[qPosition--] = (int)'_';
				sEnd[sPosition--] = (int)subject.charAt(j - 1);
				j--;
			}
		}
		
		// inserts the optimal query array into the final query array
		// inserts indels for the remaining spaces since highest possible
		// length is = query length + subject length
		while (qPosition > 0)
		{
			if (i > 0) 
				qEnd[qPosition--] = (int)query.charAt(--i);
			else 
				qEnd[qPosition--] = (int)'-';
		}
		
		// inserts the optimal subject array into the final subject array
		// inserts indels for the remaining spaces since highest possible
		// length is = query length + subject length
		while (sPosition > 0)
		{
			if (j > 0) 
				sEnd[sPosition--] = (int)subject.charAt(--j);
			else 
				sEnd[sPosition--] = (int)'-';
		}

		// answer is q length + s length long
		// so we need to catalogue parts where both sequences 
		// have a gap at the same time so that we can avoid those
		int doubleGap = 1;
		for (i = l; i >= 1; i--)
		{
			if ((char)sEnd[i] == '-' && (char)qEnd[i] == '-')
			{
				doubleGap = i + 1;
				break;
			}
		}
		
		// Printing the aligned result
		// stopping before we print out the double indel aligned portion
		System.out.println("\nThe aligned genes are :");
		for (i = doubleGap; i <= l; i++)
		{
			System.out.print((char)qEnd[i]);
		}
		System.out.print("\n");
		for (i = doubleGap; i <= l; i++)
		{
			System.out.print((char)sEnd[i]);
		}
		
		// calling the methods needed for calculating distance score of alignment
		countIndels();
		countMatches();
		countTotalCol();
		
		// calculating distance score
		calculateDistance(indelCount, matchCount, totalWOCol);
		return;
	}
	
	/**
	 * counts the number of matches among the two final sequences
	 */
	static void countMatches()
	{
		// cataloging parts where both sequences 
		// have a gap at the same time so that we can avoid those
		int both = 1;
		for(int i = qEnd.length - 1; i >= 1; i--)
		{
			if ((char)qEnd[i] == '-' && (char)sEnd[i] == '-')
			{
				both = i + 1;
				break;
			}
		}
		
		// checks for matches before the indel matches begin
		for (int j = both; j <= qEnd.length - 1; j++)
		{
			if(qEnd[j] == sEnd[j])
				matchCount++;
		}
		System.out.println("\nTotal match columns in aligned sequence: " + matchCount);
		return;
	}
	
	/**
	 * counts the number of indels present in either sequences
	 */
	static void countIndels()
	{
		// cataloging parts where both sequences 
		// have a gap at the same time so that we can avoid those
		int both = 1;
		for(int i = qEnd.length - 1; i >= 1; i--)
		{
			if ((char)qEnd[i] == '-' && (char)sEnd[i] == '-')
			{
				both = i + 1;
				break;
			}
		}
		
		// checks for indels in either sequence before indel matches begin
		for (int j = both; j <= qEnd.length - 1; j++)
		{
			if((char)qEnd[j] == '_' || (char)sEnd[j] == '_')
				indelCount++;
		}
		System.out.println("\n\nTotal columns with indels in aligned sequence: " + indelCount);
		return;
	}
	
	/**
	 * counts the total amount of aligned columns
	 */
	static void countTotalCol()
	{
		// cataloging parts where both sequences 
		// have a gap at the same time so that we can avoid those
		// increasing total columns count otherwise
		int both = 1;
		for(int i = qEnd.length - 1; i >= 1; i--)
		{
			if ((char)qEnd[i] == '-' && (char)sEnd[i] == '-')
			{
				both = i + 1;
				break;
			}
			else
			{
				totalWOCol++;
			}
		}
		return;
	}

	/**
	 * calculates the distance score of the final aligned sequences
	 * @param indel : subtracted from total columns to find total columns without indels
	 * @param match : divided by total columns without indels to find distance score
	 * @param total : total columns to be subtracted by indel columns
	 */
	static void calculateDistance(int indel, int match, int total)
	{
		// distance formula = matches/ total alignment without indel columns
		double distScore = 0;
		totalWOCol = total - indel;
		distScore = (double)match/totalWOCol;
		System.out.println("\nTotal columns w/o indels: " + totalWOCol);
		System.out.println("\nTotal Distance Score: " + distScore);
		return;
	}
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args)
	{
		boolean finish = false;
		// input strings
		while(finish == false)
		{
			Scanner fasta = new Scanner(System.in);
			System.out.println("\nInsert query fasta sequence without defline: \n");
			String query = fasta.next();
			System.out.println("\nInsert subject fasta sequence without defline: \n");
			String subject = fasta.next();

			// mismatch penalty has to be 3 instead of 1 because there
			// are no negative numbers
			int misMatchPenalty = 3;
			// gap penalty can be the same as before
			int gapPenalty = 2;
	
			// calling the function to
			// calculate the result
			findAlignment(query, subject,
					misMatchPenalty, gapPenalty);
			System.out.println("\nWould you like to enter a new alignment? y/n ");
			String ans = fasta.next();
			if(ans.contentEquals("n"))
			{
				finish = true;
				fasta.close();
			}
			indelCount = 0;
			matchCount = 0;
			totalWOCol = 0;
		}
	}
}
