import java.io.FileNotFoundException;

public class Preprocessing extends Solution {
	public static void main(String[] args) throws FileNotFoundException {
		String path = args[0];
		boolean display = false;
		if (args.length >= 2)
			display = Boolean.parseBoolean(args[1]);

		Solution s = new Solution(path);
		try {
			int elim0 = s.nElim;
			System.out.println("Quickprocessing en cours ...");
			System.out.println();
			s.quickPreprocessing();
			boolean feasible = s.FEASIBLE;
			if (!feasible) {
				System.out.println("Votre problème n'a pas de solutions.");
				System.out.println();
				return;
			}
			if (display)
				s.Display();
			int elim1 = s.nElim;
			System.out.println("IP-Processing en cours ...");
			System.out.println();
			s.IPProcessing();
			if (display)
				s.Display();
			int elim2 = s.nElim;
			System.out.println("LPProcessing en cours ...");
			System.out.println();
			s.LPProcessing();
			int elim3 = s.nElim;

			if (display)
				s.Display();
			
			System.out.println( (elim1 - elim0) +" variables éliminées au cours de Quickprocessing.");
			System.out.println( (elim2-elim1) +" variables éliminées au cours de IP-Processing.");
			System.out.println( (elim3-elim2) +" variables éliminées au cours de LP-Processing.");
			System.out.println( elim3 + " variables éliminées en tout au cours de ce prétraitement.");

		} catch (Exception e) {
			e.printStackTrace();

		}
	}
}
