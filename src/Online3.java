import java.io.FileNotFoundException;

public class Online3 extends Solution {

	public static void main(String[] args) throws FileNotFoundException {

		int iter = 10000;
		if (args.length >= 1)
			iter = (int) Integer.parseInt(args[0]);

		int pmax = 50, rmax = 100;
		double r = 0;
		double p = 0;
		int nAlert = 0;
		for (int i = 0; i < iter; i++) {
			Input ipt = new Input(4, 2, 10, 100);

			Solution sol = new Solution(ipt);
			Solution s = new Solution(ipt);
			try {
				sol.online3(pmax, rmax);
				nAlert += sol.nAlert;
				if (sol.nAlert == 1) {
					continue;
				}
			} catch (Exception e) {
				continue;
			}
			double p1 = sol.p_used;
			double r1 = sol.rate_opt;
			try {
				s.input.setPower(sol.input.getPower());
				s.input.setRate(sol.input.getRate());
				s.Greedy();
				
			} catch (Exception e) {
				continue;
			}
			double p2 = s.p_used;
			if (p2 == 0)
				continue;
			double r2 = s.rate_opt;

			p += p1 / p2;
			r += r1 / r2;
		}
		System.out.println("Ratio moyen de puissance sur " + iter + " it�rations : " + p / iter);
		System.out.println("Ratio moyen de rate sur " + iter + " it�rations : " + r / iter);
		System.out.println(
				"Nombre de dépassements de contraintes en pourcentage : " + 100 * (double) nAlert / iter + "%");
	}

}
