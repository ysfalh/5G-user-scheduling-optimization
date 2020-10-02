import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class Input {
	private int N;
	private int M;
	private int K;
	private int P;
	private int[][][] power;
	private int[][][] rate;

	public Input() {
		this.N = 0;
		this.K = 0;
		this.M = 0;
		this.P = 0;
		this.power = new int[N][K][M];
		this.rate = new int[N][K][M];
	}
	
	public Input(int N, int M, int K, int P) {
		this.N = N;
		this.K = K;
		this.M = M;
		this.P = P;
		this.power = new int[N][K][M];
		this.rate = new int[N][K][M];
	}
	
	public Input(String path) throws FileNotFoundException {
		this.setInput(path);
	}

	void setCst(int[] Cst) {
		this.N = Cst[0];
		this.M = Cst[1];
		this.K = Cst[2];
		this.P = Cst[3];
	}

	int[] getCst() {
		int[] Cst = new int[4];
		Cst[0] = this.N;
		Cst[1] = this.M;
		Cst[2] = this.K;
		Cst[3] = this.P;
		return Cst;
		
	}
	
	void setPower(int[][][] power) {
		this.power = power;
	}
	
	int[][][] getPower(){
		return this.power;
	}
	
	void setRate(int[][][] rate) {
		this.rate = rate;
	}
	
	int[][][] getRate(){
		return this.rate;
	}
	
	void setInput(String path) throws FileNotFoundException {
		int[] Cst = new int[4];
		int i = 0;
		File file = new File(path);
		Scanner s = new Scanner(file);
		while (s.hasNext() && i < 4) {
			double token = Double.parseDouble(s.next());
			Cst[i] = (int) token;
			i++;
		}
		this.setCst(Cst);
		
		
//		set input Power
		
		int[][][] inputPower = new int[N][K][M];
		
		for(int n=0; n < this.N ; n++ ) {
			for (int k = 0; k < this.K ; k++) {
				for (int m = 0; m < this.M; m++) {
					double token = Double.parseDouble(s.next());
					inputPower[n][k][m] = (int) token;
				}
			}
		}
		
		this.setPower(inputPower);
		
		int[][][] inputRate = new int[N][K][M];
		
		for(int n=0; n < this.N ; n++ ) {
			for (int k = 0; k < this.K ; k++) {
				for (int m = 0; m < this.M; m++) {
					double token = Double.parseDouble(s.next());
					inputRate[n][k][m] = (int) token;
				}
			}
		}
		
		this.setRate(inputRate);
		s.close();
	}

	public static void main(String[] args) throws FileNotFoundException {
		Input in = new Input();
		String path = "C:\\Users\\Youssef\\Downloads\\testfiles\\test1.txt";
		in.setInput(path);
		System.out.println(in.getPower()[in.N-1][in.K-1][in.M-1]);
	}
	
}
