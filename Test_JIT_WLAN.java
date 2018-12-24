import java.sql.Date;
import java.text.SimpleDateFormat;
import java.util.Arrays;

public class Test_JIT_WLAN {

	public static void main(String [] args){
		int Num_user = 2;
		int Num_AP = 3;
		int Total_channel = 1;
		double H = 1;
		
		JITWLAN solver = new JITWLAN (Num_user,Num_AP,Total_channel,H);
		
		// overlap
		solver.setOverlap_G(0, 1, 1);
		solver.setOverlap_G(1, 0, 1);
		solver.setOverlap_G(1, 2, 1);
		solver.setOverlap_G(2, 1, 1);
		
		//weight
		solver.setWeight(0, 0, 16.2);
		solver.setWeight(1, 0, 4.5);
		solver.setWeight(1, 1, 6.8);
		solver.setWeight(2, 1, 10.2);
		
		//cost
		solver.setCost(0, 0.8);
		solver.setCost(1, 0.1);
		solver.setCost(2, 0.2);
		
		//solve LP
		System.out.println("------RoundMax------");
		double [] sol = solver.solveLP();
		int [][] solution_max = solver.Round_max(sol);
		for (int i = 0; i<Num_AP;i++){
			System.out.println();
			for (int j = 0; j<Num_user; j++){
				System.out.print(solution_max[i][j]);
			}
		}
		
		System.out.println();
		double mwis = solver.Cost_function(solution_max);		
		System.out.println();
		System.out.println(mwis);
		
		//another
		System.out.println("------RoundMin------");
		int [][] solution_min = solver.Round_min(Arrays.copyOfRange(sol, Num_user*Num_AP*Total_channel, sol.length));
		
		for (int i = 0; i<Num_AP;i++){
			System.out.println();
			for (int j = 0; j<Num_user; j++){
				System.out.print(solution_min[i][j]);
			}
		}
		System.out.println();
		
		double msc = solver.Cost_function(solution_min);	
		System.out.println(msc);
		int [][] solution_final = new int[Num_AP][Num_user];
		
		if(msc>mwis){
			solution_final = solution_min.clone();
		}else{
			solution_final = solution_max.clone();
		}
		
		System.out.println("------Final---------");
		for (int i = 0; i<Num_AP;i++){
			System.out.println();
			for (int j = 0; j<Num_user; j++){
				System.out.print(solution_final[i][j]);
			}
		}

		
	}
}
