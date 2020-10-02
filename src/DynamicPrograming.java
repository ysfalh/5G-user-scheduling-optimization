import java.io.FileNotFoundException;

public class DynamicPrograming extends Solution {
	public static void main(String[] args) throws FileNotFoundException {
		String path = args[0];
//		String path = "/home/sariah/Documents/eclipse/testfiles/test5.txt";
		boolean display = false;
		if (args.length >= 2)
			display = Boolean.parseBoolean(args[1]);

		Solution s = new Solution(path);
		try { 
			System.out.println("DPPower en cours ...");
			long startTime = System.currentTimeMillis();
			s.DPPower();
			long endTime = System.currentTimeMillis();
			
			long timeDPPower = endTime - startTime;
			System.out.println("Execution de DPPower en " + timeDPPower + " ms.");
			System.out.println();
			int U = s.DPr + 500;
			
			if (display) s.Display();
			
			s = new Solution(path);
			
			System.out.println("DPRate(U = " + U + ") en cours ...");
			long startTime2 = System.currentTimeMillis();
			s.DPRate(U);
			long endTime2 = System.currentTimeMillis();
			long timeDPRate = endTime2 - startTime2;
			
			System.out.println("Execution de DPRate en " + timeDPRate + " ms.");
			
			if (display) s.Display();
			
		} catch (Exception e) {
			e.printStackTrace();

		}
	}
}
