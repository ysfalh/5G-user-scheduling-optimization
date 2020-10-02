import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Stack;

import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

public class Solution extends Input {
	Input input;
	double[][][] X;
	int N, K, M, P;
	int nElim = 0;
	boolean FEASIBLE = true;
	double p_used; // puissance utilis�e
	double rate_opt; // rate produit
	static final int ELIM = -1000;
	int nAlert = 0;
	int DPr = 0;

	static final int BigInt = 100000;
	static final int LowInt = -100000;

	Solution() {
		input = new Input();
		this.N = input.getCst()[0];
		this.M = input.getCst()[1];
		this.K = input.getCst()[2];
		this.P = input.getCst()[3];
		X = new double[N][K][M];
	}

	Solution(Input input) {
		this.input = input;
		this.N = input.getCst()[0];
		this.M = input.getCst()[1];
		this.K = input.getCst()[2];
		this.P = input.getCst()[3];
		X = new double[N][K][M];
	}

	public Solution(String path) throws FileNotFoundException {
		this.input = new Input(path);
		this.N = this.input.getCst()[0];
		this.M = this.input.getCst()[1];
		this.K = this.input.getCst()[2];
		this.P = this.input.getCst()[3];
		X = new double[N][K][M];
	}

	void Display() {

		System.out.println(this.N);
		System.out.println(this.M);
		System.out.println(this.K);
		System.out.println(this.P);

		for (int n = 0; n < this.N; n++) {
			for (int k = 0; k < this.K; k++) {
				String toprint = new String();
				for (int m = 0; m < this.M; m++) {
					if (this.X[n][k][m] == ELIM) {
						toprint += " _";
					} else {
						toprint += " " + String.valueOf(this.X[n][k][m]);
					}
				}
				System.out.println(toprint);
			}
		System.out.println();
		}
	}

	void Display(int n) {

		System.out.println(this.N);
		System.out.println(this.M);
		System.out.println(this.K);
		System.out.println(this.P);

		for (int k = 0; k < this.K; k++) {
			String toprint = new String();
			for (int m = 0; m < this.M; m++) {
				if (this.X[n][k][m] == ELIM) {
					toprint += " _ ";
				} else {
					toprint += " " + String.valueOf(this.X[n][k][m]);
					System.out.println("n = " + n + " k = " + k + " m = " + m + "p = " + this.input.getPower()[n][k][m]
							+ " r = " + this.input.getRate()[n][k][m]);
				}
			}
			System.out.println(toprint);
		}
		System.out.println();

	}

	void quickPreprocessing() {

		int minpower = 0;
		int[] minpowerlist = this.minPowerAux();

		for (int n = 0; n < this.N; n++) {

			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					if (this.input.getPower()[n][k][m] > this.P) {
						if (this.X[n][k][m] != ELIM) {
							this.nElim++;
						}
						this.X[n][k][m] = ELIM;
					}

				}
			}
			minpower += minpowerlist[n];

		}

		if (minpower > this.P)
			this.FEASIBLE = false;
	}

	int[] minPowerAux() {
		int[][][] power = this.input.getPower();
		int N = this.N;
		int K = this.K;
		int[] minPowerList = new int[N];
		for (int i = 0; i < N; i++) {
			int m = power[i][0][0];
			for (int k = 0; k < K; k++) {
				if (m > power[i][k][0]) {
					m = power[i][k][0];
				}
			}
			minPowerList[i] = m;
		}
		return minPowerList;
	}

	void IPProcessing() throws Exception {
		this.quickPreprocessing();
		if (!this.FEASIBLE) {
			throw new Exception("Il n'y a pas de solution pour cette instance.");
		}

		ArrayList<Map<Integer[], Integer>> mapList = new ArrayList<>();

		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = new HashMap<>();
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					Integer[] indices = new Integer[2];
					indices[0] = k;
					indices[1] = m;
					if (!isElim(indices, n))
						map.put(indices, this.input.getPower()[n][k][m]);
				}
			}
			map = MapUtil.sortByValue(map);
			mapList.add(map);
		}

		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = mapList.get(n);

			Integer[][] myArray = new Integer[map.keySet().size()][];
			Integer[][] keyList = map.keySet().toArray(myArray);
			int l = 0;
			int i = 1;
			while (l + i < keyList.length) {

				Integer[][] indices = new Integer[2][2];
				indices[0] = keyList[l];
				indices[1] = keyList[l + i];

				int k1 = indices[0][0];
				int m1 = indices[0][1];
				int k2 = indices[1][0];
				int m2 = indices[1][1];

				int r1 = this.input.getRate()[n][k1][m1];
				int r2 = this.input.getRate()[n][k2][m2];

				int p1 = this.input.getPower()[n][k1][m1];
				int p2 = this.input.getPower()[n][k2][m2];

				if (r2 <= r1) {
					if (this.X[n][k2][m2] != ELIM) {
						this.nElim++;
					}
					this.X[n][k2][m2] = ELIM;
					i++;
				} else {
					if (p1 == p2) {
						if (this.X[n][k1][m1] != ELIM) {
							this.nElim++;
						}
						this.X[n][k1][m1] = ELIM;
						l += i;
						i = 1;
					} else {
						l += i;
						i = 1;
					}
				}
			}

		}

	}

	boolean angle(int n, Integer[] Pivot, Integer[] P1, Integer[] P2) {
		Integer kPivot = Pivot[0];
		Integer mPivot = Pivot[1];

		Integer kP1 = P1[0];
		Integer mP1 = P1[1];

		Integer kP2 = P2[0];
		Integer mP2 = P2[1];

		int pPivot = this.input.getPower()[n][kPivot][mPivot];
		int rPivot = this.input.getRate()[n][kPivot][mPivot];

		int pP1 = this.input.getPower()[n][kP1][mP1];
		int rP1 = this.input.getRate()[n][kP1][mP1];

		int pP2 = this.input.getPower()[n][kP2][mP2];
		int rP2 = this.input.getRate()[n][kP2][mP2];

		int det = (pP1 - pPivot) * (rP2 - rPivot) - (pP2 - pPivot) * (rP1 - rPivot);
		return (Math.signum(det) > 0);
	}

	void LPProcessing() throws Exception {
		this.IPProcessing();

		ArrayList<Map<Integer[], Integer>> mapList = new ArrayList<>();

		for (int n = 0; n < this.N; n++) {

			Map<Integer[], Integer> map = new HashMap<>();
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					if (this.X[n][k][m] != ELIM) {
						Integer[] indices = new Integer[2];
						indices[0] = k;
						indices[1] = m;
						if (!isElim(indices, n))
							map.put(indices, this.input.getPower()[n][k][m]);
					}
				}
			}
			map = MapUtil.sortByValue(map);
			mapList.add(map);
		}

		for (int n = 0; n < this.N; n++) {

			Map<Integer[], Integer> map = mapList.get(n);

			Integer[][] myArray = new Integer[map.keySet().size()][];
			Integer[][] keyList = map.keySet().toArray(myArray);
			int L = keyList.length;

			Stack<Integer[]> enveloppe = new Stack<Integer[]>();

			enveloppe.push(keyList[0]);
			enveloppe.push(keyList[1]);

			for (int l = 2; l < L; l++) {
				boolean valide = false;
				enveloppe.push(keyList[l]);

				while (!valide && enveloppe.size() > 2) {

					Integer[] P3 = enveloppe.pop();
					Integer[] P2 = enveloppe.pop();
					Integer[] P1 = enveloppe.pop();

					if (angle(n, P1, P2, P3)) {
						enveloppe.push(P1);
						enveloppe.push(P3);

						if (this.X[n][P2[0]][P2[1]] != ELIM) {
							this.nElim++;
						}
						this.X[n][P2[0]][P2[1]] = ELIM;
					} else {
						enveloppe.push(P1);
						enveloppe.push(P2);
						enveloppe.push(P3);
						valide = true;
					}
				}
			}

		}

	}

	void Greedy() throws Exception {
		this.LPProcessing();

		ArrayList<Map<Integer[], Integer>> mapList = new ArrayList<>();

		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = new HashMap<>();
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					Integer[] indices = new Integer[2];
					indices[0] = k;
					indices[1] = m;
					if (!this.isElim(indices, n))
						map.put(indices, this.input.getPower()[n][k][m]);
				}
			}
			map = MapUtil.sortByValue(map);
			mapList.add(map);
		}
//		Initialise le budget utilis�
		double budget = 0;
		double rate = 0;

		Integer[] Candidats = new Integer[this.N];
		double[] Efficacite = new double[this.N];
		Integer[][][] bigKeyList = new Integer[this.N][][];
		Integer[][] keyList;
		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = mapList.get(n);

			Integer[][] myArray = new Integer[map.keySet().size()][];
			keyList = map.keySet().toArray(myArray);
			bigKeyList[n] = keyList;

//			Initialise la liste des candidats pour la solution optimale
//			Candidats[n] est l'indice dans liste keyList de l'indice du candidat
			Candidats[n] = 0;
			int k = bigKeyList[n][Candidats[n]][0];
			int m = bigKeyList[n][Candidats[n]][1];
			budget += this.input.getPower()[n][k][m];
		}

//		target d�signe l'indice du candidat promu � l'�tape en question
		int target = 0;
		while (budget <= this.P) {
			target = 0;
			for (int n = 0; n < this.N; n++) {
				keyList = bigKeyList[n];
				if (Candidats[n] == keyList.length - 1)
					Efficacite[n] = -1;
				else {
					int k1 = keyList[Candidats[n]][0];
					int m1 = keyList[Candidats[n]][1];
					int k2 = keyList[Candidats[n] + 1][0];
					int m2 = keyList[Candidats[n] + 1][1];

					int r1 = this.input.getRate()[n][k1][m1];
					int r2 = this.input.getRate()[n][k2][m2];

					int p1 = this.input.getPower()[n][k1][m1];
					int p2 = this.input.getPower()[n][k2][m2];

					Efficacite[n] = ((double) (r2 - r1)) / (p2 - p1);
				}
				if (Efficacite[n] != -1 && Efficacite[target] <= Efficacite[n])
					target = n;
			}
			if (Efficacite[target] == -1)
				break;
			int k = bigKeyList[target][Candidats[target]][0];
			int m = bigKeyList[target][Candidats[target]][1];
			Candidats[target]++;
			int kpromu = bigKeyList[target][Candidats[target]][0];
			int mpromu = bigKeyList[target][Candidats[target]][1];

			budget += this.input.getPower()[target][kpromu][mpromu] - this.input.getPower()[target][k][m];

		}

		for (int n = 0; n < this.N; n++) {
			int k = bigKeyList[n][Candidats[n]][0];
			int m = bigKeyList[n][Candidats[n]][1];
			if (n == target && budget > this.P) {
				int kpaspromu = bigKeyList[target][Candidats[target] - 1][0];
				int mpaspromu = bigKeyList[target][Candidats[target] - 1][1];
				double p_dispo = this.P - budget + this.input.getPower()[target][k][m];
				this.X[n][k][m] = (p_dispo - this.input.getPower()[target][kpaspromu][mpaspromu])
						/ ((double) (this.input.getPower()[target][k][m]
								- this.input.getPower()[target][kpaspromu][mpaspromu]));
				this.X[n][kpaspromu][mpaspromu] = 1 - this.X[n][k][m];

				budget = budget - this.input.getPower()[target][k][m]
						+ this.X[n][k][m] * this.input.getPower()[target][k][m]
						+ this.X[n][kpaspromu][mpaspromu] * this.input.getPower()[target][kpaspromu][mpaspromu];
				rate += this.input.getRate()[n][k][m] * this.X[n][k][m]
						+ this.input.getRate()[n][kpaspromu][mpaspromu] * this.X[n][kpaspromu][mpaspromu];
			} else {
				this.X[n][k][m] = 1;
				rate += this.input.getRate()[n][k][m];
			}
		}
		this.p_used = budget;
		this.rate_opt = rate;

	}

	void LPSolver1() throws Exception {
		ArrayList<Integer[]> variables = new ArrayList<>();
		this.LPProcessing();
//		DECISION VARIABLES ARRAYLIST
		for (int n = 0; n < this.N; n++) {
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					Integer[] indices = new Integer[2];
					indices[0] = k;
					indices[1] = m;
					if (!this.isElim(indices, n)) {
						variables.add(new Integer[] { n, k, m });
					}
				}
			}
		}

		double[] coeff = new double[variables.size()];
		for (int i = 0; i < coeff.length; i++) {
			int n = variables.get(i)[0];
			int k = variables.get(i)[1];
			int m = variables.get(i)[2];
			coeff[i] = this.input.getRate()[n][k][m];
		}

		LinearProgram lp = new LinearProgram(coeff);
		for (int i = 0; i < coeff.length; i++) {
			int n = variables.get(i)[0];
			int k = variables.get(i)[1];
			int m = variables.get(i)[2];
			lp.addConstraint(new LinearBiggerThanEqualsConstraint(contrainte1(i, variables.size()), 0.0,
					"c[" + n + "][" + k + "][" + m + "]"));
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(contrainte1(i, variables.size()), 1.0,
					"d[" + n + "][" + k + "][" + m + "]"));
		}

		double[][] bigC2 = new double[this.N][variables.size()];
		int j = 0;
		while (j < variables.size()) {
			int n = variables.get(j)[0];
			bigC2[n][j] = 1;
			j++;
		}

		double[] c3 = new double[variables.size()];
		for (int i = 0; i < coeff.length; i++) {
			int n = variables.get(i)[0];
			int k = variables.get(i)[1];
			int m = variables.get(i)[2];
			c3[i] = this.input.getPower()[n][k][m];
		}

		for (int n = 0; n < this.N; n++) {
			lp.addConstraint(new LinearBiggerThanEqualsConstraint(bigC2[n], 1.0, "c21"));
			lp.addConstraint(new LinearSmallerThanEqualsConstraint(bigC2[n], 1.0, "c22"));
		}

		lp.addConstraint(new LinearSmallerThanEqualsConstraint(c3, this.P, "c3"));
		lp.setMinProblem(false);

		LinearProgramSolver solver = SolverFactory.newDefault();
		double[] sol = solver.solve(lp);

		double s = 0;
		double budget = 0;
		for (int i = 0; i < sol.length; i++) {
			int n = variables.get(i)[0];
			int k = variables.get(i)[1];
			int m = variables.get(i)[2];
			this.X[n][k][m] = sol[i];
			s += sol[i] * this.input.getRate()[n][k][m];
			budget += sol[i] * this.input.getPower()[n][k][m];
		}

		this.p_used = budget;
		this.rate_opt = s;
		System.out.println("Budget �puis� : " + budget + "/" + this.P);
		System.out.println("Rate optimal : " + s);

	}

	void powernp(int[][][] tableRate, int[][][] tableRatestatic, Integer[][] keyList, int n, int p) {
		// retourne la case (n,p) du tableau pour l'approche dynamique, pour n >= 1

		int max = 0;
		int x = 0;

		int kchoosen = 0;
		int mchoosen = 0;
		int prevpchoosen = 0;

		for (int q = 0; q <= p; q++) {
			int r = tableRatestatic[n][p - q][0];
			int k = tableRatestatic[n][p - q][1];
			int m = tableRatestatic[n][p - q][2];

			x = tableRate[n - 1][q][0] + r;
			int prevp = q;

			if (max <= x) {
				max = x;
				kchoosen = k;
				mchoosen = m;
				prevpchoosen = prevp;

			}
		}
		tableRate[n][p][0] = max;
		tableRate[n][p][1] = kchoosen;
		tableRate[n][p][2] = mchoosen;
		tableRate[n][p][3] = prevpchoosen;

	}

	void DPPower() throws Exception {

		this.LPProcessing();

		ArrayList<Map<Integer[], Integer>> mapList = new ArrayList<>();

		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = new HashMap<>();
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					Integer[] indices = new Integer[2];
					indices[0] = k;
					indices[1] = m;
					if (!this.isElim(indices, n))
						map.put(indices, this.input.getPower()[n][k][m]);
				}
			}
			map = MapUtil.sortByValue(map);
			mapList.add(map);
		}

		int[][][] tableRate = new int[this.N][this.P + 1][4];
		int[][][] tableRatestatic = new int[this.N][this.P + 1][3];

		// initialisation de tableRate

		int n = 0;
		Map<Integer[], Integer> map = mapList.get(n);
		Integer[][] keyList;
		Integer[][] myArray = new Integer[map.keySet().size()][];
		keyList = map.keySet().toArray(myArray);

		int l = 0;
		int k = keyList[l][0];
		int m = keyList[l][1];
		int r = LowInt;
		int Pmax = this.input.getPower()[n][k][m];
		int p = 0;
		while (p <= this.P && Pmax <= this.P) {
			while (p < Pmax) {
				tableRate[n][p][0] = r;
				tableRate[n][p][1] = k;
				tableRate[n][p][2] = m;
				p++;
				if (p == this.P) {
					tableRate[n][p][0] = r;
					tableRate[n][p][1] = k;
					tableRate[n][p][2] = m;
					p++;
				}
			}

			k = keyList[l][0];
			m = keyList[l][1];
			r = this.input.getRate()[n][k][m];

			if (l == keyList.length - 1) {
				Pmax = this.P;

			} else {
				l++;
				Pmax = this.input.getPower()[n][keyList[l][0]][keyList[l][1]];
			}
		}
		// initialisation de tableRatestatic

		for (n = 0; n < this.N; n++) {
			map = mapList.get(n);
			myArray = new Integer[map.keySet().size()][];
			keyList = map.keySet().toArray(myArray);

			l = 0;
			k = keyList[l][0];
			m = keyList[l][1];

			r = LowInt;
			Pmax = this.input.getPower()[n][k][m];
			p = 0;
			while (p <= this.P && Pmax <= this.P) {
				while (p < Pmax) {
					tableRatestatic[n][p][0] = r;
					tableRatestatic[n][p][1] = k;
					tableRatestatic[n][p][2] = m;
					p++;
					if (p == this.P) {
						tableRatestatic[n][p][0] = r;
						tableRatestatic[n][p][1] = k;
						tableRatestatic[n][p][2] = m;
						p++;
					}
				}

				k = keyList[l][0];
				m = keyList[l][1];
				r = this.input.getRate()[n][k][m];
				if (l == keyList.length - 1) {
					Pmax = this.P;

				} else {
					l++;
					Pmax = this.input.getPower()[n][keyList[l][0]][keyList[l][1]];
				}
			}
		}
		// It�ration par approche dynamique

		for (n = 1; n < this.N; n++) {
			map = mapList.get(n);
			myArray = new Integer[map.keySet().size()][];
			keyList = map.keySet().toArray(myArray);

			for (p = 0; p <= this.P; p++) {
				this.powernp(tableRate, tableRatestatic, keyList, n, p);
			}
		}

		p = this.P;
		for (n = this.N - 1; n >= 0; n--) {
			k = tableRate[n][p][1];
			m = tableRate[n][p][2];
			p = tableRate[n][p][3];

			this.X[n][k][m] = 1;
		}
		int powerused = 0;

		for (n = 0; n < this.N; n++) {
			for (k = 0; k < this.K; k++) {
				for (m = 0; m < this.M; m++) {
					if (this.X[n][k][m] == 1)
						powerused += this.input.getPower()[n][k][m];
				}
			}
		}

		System.out.println("DPPower: Puissance utilisée = " + powerused + "/" + this.P
				+ " Débit max = " + tableRate[this.N - 1][this.P][0]);
		this.DPr = tableRate[this.N - 1][this.P][0];
	}

	void ratenp(int[][][] tablePower, int[][][] tablePowerstatic, Integer[][] keyList, int n, int r, int U) {
		// retourne la case (n,p) du tableau pour l'approche dynamique, pour n >= 1

		int min = tablePowerstatic[n][U][0];
		int x;
		int kchoosen = 0;
		int mchoosen = 0;
		int prevrchoosen = 0;

		for (int q = 1; q < r; q++) {
			int p = tablePowerstatic[n][r - q][0];
			int k = tablePowerstatic[n][r - q][1];
			int m = tablePowerstatic[n][r - q][2];

			x = tablePower[n - 1][q][0] + p;
			int prevr = q;

			if (x != 0) {
				if (min >= x) {
					min = x;
					kchoosen = k;
					mchoosen = m;
					prevrchoosen = prevr;
				}
			}

		}
		tablePower[n][r][0] = min;
		tablePower[n][r][1] = kchoosen;
		tablePower[n][r][2] = mchoosen;
		tablePower[n][r][3] = prevrchoosen;

	}

	void DPRate(int U) throws Exception {

		this.LPProcessing();

		ArrayList<Map<Integer[], Integer>> mapList = new ArrayList<>();

		for (int n = 0; n < this.N; n++) {
			Map<Integer[], Integer> map = new HashMap<>();
			for (int k = 0; k < this.K; k++) {
				for (int m = 0; m < this.M; m++) {
					Integer[] indices = new Integer[2];
					indices[0] = k;
					indices[1] = m;
					if (!this.isElim(indices, n))
						map.put(indices, this.input.getRate()[n][k][m]);
				}
			}
			map = MapUtil.sortByValue(map);
			mapList.add(map);
		}

		int[][][] tablePower = new int[this.N][U + 1][4];
		int[][][] tablePowerstatic = new int[this.N][U + 1][3];

		// initialisation de tableRate

		int n = 0;
		Map<Integer[], Integer> map = mapList.get(n);
		Integer[][] keyList;
		Integer[][] myArray = new Integer[map.keySet().size()][];
		keyList = map.keySet().toArray(myArray);

		int l = 0;
		int k = keyList[l][0];
		int m = keyList[l][1];
		int p = this.input.getPower()[n][k][m];
		int Rmax = this.input.getRate()[n][k][m];
		int r = 0;
		while (r <= U && Rmax <= U) {
			while (r <= Rmax) {
				tablePower[n][r][0] = p;
				tablePower[n][r][1] = k;
				tablePower[n][r][2] = m;
				r++;
				if (r == U) {
					tablePower[n][r][0] = p;
					tablePower[n][r][1] = k;
					tablePower[n][r][2] = m;
					r++;
				}
			}

			if (l == keyList.length - 1) {
				Rmax = U;
				p = BigInt;

			} else {
				l++;
				k = keyList[l][0];
				m = keyList[l][1];
				Rmax = this.input.getRate()[n][k][m];
				p = this.input.getPower()[n][k][m];
			}

		}

		for (n = 1; n < this.N; n++) {
			map = mapList.get(n);
			myArray = new Integer[map.keySet().size()][];
			keyList = map.keySet().toArray(myArray);

			l = 0;
			k = keyList[l][0];
			m = keyList[l][1];
			p = this.input.getPower()[n][k][m];
			Rmax = this.input.getRate()[n][k][m];
			r = 0;
			while (r <= U && Rmax <= U) {
				while (r <= Rmax) {
					tablePowerstatic[n][r][0] = p;
					tablePowerstatic[n][r][1] = k;
					tablePowerstatic[n][r][2] = m;
					r++;
					if (r == U) {
						tablePowerstatic[n][r][0] = p;
						tablePowerstatic[n][r][1] = k;
						tablePowerstatic[n][r][2] = m;
						r++;
					}
				}

				if (l == keyList.length - 1) {
					Rmax = U;
					p = BigInt;
				} else {
					l++;
					k = keyList[l][0];
					m = keyList[l][1];
					Rmax = this.input.getRate()[n][k][m];
					p = this.input.getPower()[n][k][m];
				}
			}
		}

		for (n = 1; n < this.N; n++) {
			map = mapList.get(n);
			myArray = new Integer[map.keySet().size()][];
			keyList = map.keySet().toArray(myArray);

			for (r = 0; r <= U; r++) {
				this.ratenp(tablePower, tablePowerstatic, keyList, n, r, U);
			}
		}

		r = U;
		while (tablePower[this.N - 1][r][0] >= this.P) {
			r--;
		}
		int R = r;
		for (n = this.N - 1; n >= 0; n--) {
			k = tablePower[n][r][1];
			m = tablePower[n][r][2];
			r = tablePower[n][r][3];

			this.X[n][k][m] = 1;
		}
		System.out.println("DPRate: Puissance utilsée = " + tablePower[this.N - 1][R][0] + "/"
				+ this.P + " rate max = " + R);

	}
//fonction auxilliaire au LPSolver pour l'implementation de la contrainte 1
	double[] contrainte1(int i, int size) {
		double[] list = new double[size];
		list[i] = 1;

		return list;
	}

	void online3(int pmax, int rmax) throws Exception {
		int n_activ = 0;
		int[] activated = new int[N];
		for (int k = 0; k < this.K; k++) {
			this.newUser(k, pmax, rmax);

			if (n_activ < this.N) {
				int[] indices = new int[2];
				double efficacite = -1;
				for (int n = 0; n < this.N; n++) {
					for (int m = 0; m < this.M; m++) {
						int r = this.input.getRate()[n][k][m];
						int p = this.input.getPower()[n][k][m];
						double candidate = ((double) r) / p;

						if (activated[n] == 0 && p <= this.P - p_used && efficacite <= candidate) {
							efficacite = candidate;
							indices[0] = n;
							indices[1] = m;
						}
					}
				}
				if (efficacite != -1) {
					int n = indices[0];
					int m = indices[1];
					activated[n] = 1;
					n_activ++;
					this.X[n][k][m] = 1;
					this.p_used += this.input.getPower()[n][k][m];
					this.rate_opt += this.input.getRate()[n][k][m];
				}
			}
		}
		if (p_used > P)
			System.out.println("P = " + P + " ALERT budget " + p_used);
		if (n_activ < N)
			this.nAlert++;
//		System.out.println("------------------------------------");
//		System.out.println("Nombre de canaux activ�s : " + n_activ);
//		System.out.println("Budget �puis� : " + p_used + "/" + P);
//		System.out.println("Rate atteint : " + rate);
	}

	void online2(int pmax, int rmax) throws Exception {
		int n_activ = 0;
		int[] activated = new int[N];
		for (int k = 0; k < this.K; k++) {
			this.newUser(k, pmax, rmax);

			if (n_activ < this.N && (k >= K - (N - n_activ) || worth(k, pmax, rmax))) {
				for (int t = 0; t < N - n_activ; t++) {
					Map<Integer[], Integer> map = new HashMap<>();
					for (int n = 0; n < this.N; n++) {
						for (int m = 0; m < this.M; m++) {
							Integer[] indices = new Integer[2];
							indices[0] = n;
							indices[1] = m;
							int p = this.input.getPower()[n][k][m];
							if (p <= this.P - p_used && activated[n] == 0)
								map.put(indices, -this.input.getRate()[n][k][m]);
						}
					}
					map = MapUtil.sortByValue(map);
					Integer[][] myArray = new Integer[map.keySet().size()][];
					Integer[][] keyList = map.keySet().toArray(myArray);

					if (keyList.length == 0)
						continue;
					Integer[] ind = new Integer[2]; // indices[0] = n, indices[1] = m
					ind = keyList[0];
					int n = ind[0];
					int m = ind[1];
					int p = this.input.getPower()[n][k][m];
					int r = this.input.getRate()[n][k][m];

					Random rand = new Random();
					int r_rand = rand.nextInt(rmax) + 1;
					int p_rand = rand.nextInt(pmax) + 1;

					int i = 0;
					while (i < keyList.length && ((double) r / p <= (double) r_rand / p_rand)) {
						ind = keyList[i++];
						n = ind[0];
						m = ind[1];
						p = this.input.getPower()[n][k][m];
						r = this.input.getRate()[n][k][m];
					}

					if (k >= K - (N - n_activ) || (double) r / p > (double) r_rand / p_rand) {
						activated[n] = 1;
						n_activ++;
						this.X[n][k][m] = 1;
						this.p_used += this.input.getPower()[n][k][m];
						this.rate_opt += this.input.getRate()[n][k][m];
					}
				}
			}
		}
		if (p_used > P + 0)
			this.nAlert++;
		if (n_activ < N)
			this.nAlert++;
	}

	void online(int alpha, int pmax, int rmax) throws Exception {
		int n_activ = 0;
		int[] activated = new int[N];
		for (int k = 0; k < this.K; k++) {
			this.newUser(k, pmax, rmax);

			if (n_activ < this.N && (k >= K - (N - n_activ) || worth(k, pmax, rmax))) {
				int[] indices = new int[2];
				double efficacite = -1;
				for (int n = 0; n < this.N; n++) {
					for (int m = 0; m < this.M; m++) {
						int r = this.input.getRate()[n][k][m];
						int p = this.input.getPower()[n][k][m];
						double candidate = ((double) r) / p;

						if (k >= K - (N - n_activ) && activated[n] != 1
								&& p <= this.input.getPower()[indices[0]][k][indices[1]]) {
							efficacite = candidate;
							indices[0] = n;
							indices[1] = m;
						} else if (activated[n] != 1 && p <= this.P - p_used && r >= (alpha / 100.0) * rmax
								&& efficacite <= candidate) {
							efficacite = candidate;
							indices[0] = n;
							indices[1] = m;
						}
					}
				}
				if (k >= K - (N - n_activ) || efficacite != -1) {
					n_activ++;
					int n = indices[0];
					int m = indices[1];

					this.X[n][k][m] = 1;
					activated[n] = 1;
					this.p_used += this.input.getPower()[n][k][m];
					this.rate_opt += this.input.getRate()[n][k][m];
				}
			}
		}
		if (p_used > P + 0) {
//			System.out.println("P = " + P + " ALERT budget " + p_used);
			this.nAlert++;
		}
		if (n_activ < N) {
//			System.out.println("ALERT canaux " + n_activ);
			this.nAlert++;
		}
//		System.out.println("------------------------------------");
//		System.out.println("Nombre de canaux activ�s : " + n_activ);
//		System.out.println("Budget �puis� : " + p_used + "/" + P);
//		System.out.println("Rate atteint : " + rate);
	}
//fonction appelée par l'algorithme online pour créer un nouvel utilisateur
	void newUser(int k, int pmax, int rmax) {
		for (int n = 0; n < this.N; n++) {
			for (int m = 0; m < this.M; m++) {
				Random rand = new Random();
				int p = rand.nextInt(pmax) + 1;
				int r = rand.nextInt(rmax) + 1;

				this.input.getPower()[n][k][m] = p;
				this.input.getRate()[n][k][m] = r;
				this.X[n][k][m] = 0;
			}
		}
	}

//	donne valeur maximale de puissance � ne pas d�passer pour garder une somme de puissance
//	similaire � celle du greedy algorithm
	double p_reg(int k, int n_activ, double p_used) throws Exception {
		if (k == 0)
			return this.P;

		Input input = new Input(n_activ + 1, M, k + 1, P);
		input.setPower(this.input.getPower());
		input.setRate(this.input.getRate());
		Solution s = new Solution(input);
		s.Greedy();

		return s.p_used - this.p_used;
	}

	double r_reg(int k, int n_activ, double p_used) throws Exception {
		if (k == 0)
			return 0;

		Input input = new Input(n_activ + 1, M, k + 1, P);
		input.setPower(this.input.getPower());
		input.setRate(this.input.getRate());
		Solution s = new Solution(input);
		s.LPSolver1();

		return s.rate_opt - this.rate_opt;
	}

//	vrai si l'efficacite maximale est sup�rieure � une efficacit� tir�e al�atoirement
	boolean worth(int k, int pmax, int rmax) {
		double efficacite = 0;
		for (int n = 0; n < this.N; n++) {
			for (int m = 0; m < this.M; m++) {
				double candidate = ((double) this.input.getRate()[n][k][m]) / this.input.getPower()[n][k][m];
				if (efficacite <= candidate)
					efficacite = candidate;
			}
		}

		double mean = 0;
		for (int i = 0; i < rmax; i++) {
			for (int j = 0; j < pmax; j++) {
				mean += (double) (i + 1) / (j + 1);
			}
		}
		mean = mean / (pmax * rmax);
		return efficacite >= (mean/1);
	}

// Indique si le triple [n][k][m] a �t� �limin�	
	boolean isElim(Integer[] indices, Integer n) {
		int k = indices[0];
		int m = indices[1];
		return (this.X[n][k][m] == ELIM);
	}
}

//Impl�mente une fonction qui permet d'ordonner une Map par valeurs d'entr�es
class MapUtil {
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
		List<Entry<K, V>> list = new ArrayList<>(map.entrySet());
		list.sort(Entry.comparingByValue());

		Map<K, V> result = new LinkedHashMap<>();
		for (Entry<K, V> entry : list) {
			result.put(entry.getKey(), entry.getValue());
		}

		return result;
	}
}
