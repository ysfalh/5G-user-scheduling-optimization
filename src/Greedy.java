import java.io.FileNotFoundException;

public class Greedy {

	public static void main(String[] args) throws FileNotFoundException {
//		String path = "/home/sariah/Documents/eclipse/testfiles/test3.txt";
		String path = args[0];
		boolean display = false;
		if (args.length >= 2)
			display = Boolean.parseBoolean(args[1]);

		Solution s = new Solution(path);
		try {
			long startTime = System.currentTimeMillis();
			s.Greedy();
			long endTime = System.currentTimeMillis();
			long timegreedy = endTime - startTime;
			System.out.println("Execution de Greedy en " + timegreedy + " ms.");
			System.out.println();
			System.out.println("Budget épuisé = " + s.p_used+"/"+ s.P+" Débit max = " + s.rate_opt );
			System.out.println();
			if (display) s.Display();
		} catch (Exception e) {
			e.printStackTrace();

		}

	}

}
